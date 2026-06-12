"""
Batch Mixing Power Analysis — Single-Cell Proteomics
=====================================================
Tier 1: How many cells per plate and how many plates are needed for
        plate-effect correction to produce good batch mixing while
        preserving cell-type separation?

Tier 2: Rare cell simulation — using real K562 cells as a rare population
        embedded in a HEK293T-dominant background, under different
        distributions across plates.

Input : adata_protein.h5ad  (cells x proteins, log2 intensities in layer 'X_log2')

Output (Tier 1):
  power_results_batch_mixing.csv   — scalar metrics per iteration, with
                                     before/after correction columns and
                                     ARI/NMI from unsupervised clustering
  per_cell_scores.parquet          — per-cell silhouette + kBET scores
                                     (corrected only, to keep size manageable)
  iterations/<condition>.csv       — incremental crash-safe CSVs
  reference/reference_metrics.csv  — full-data reference run metrics
  reference/reference_per_cell.parquet

Output (Tier 2):
  rare_cell/rare_cell_results.csv  — scalar metrics per condition
  rare_cell/per_cell_scores.parquet

Design decisions:
- Before/after comparison: metrics computed on both raw (imputed, not corrected)
  and batch-corrected matrix. Delta columns capture the improvement.
  Per-cell parquet stores corrected scores only to keep storage manageable.
- Unsupervised clustering: k-means with k=2, ARI and NMI vs true cell type
  labels. Evaluates whether unsupervised methods recover known biology.
- Rare cell simulation: real K562 cells used as rare population (not synthetic)
  since the dataset provides real protein-level signal. Distribution across
  plates is varied to probe batch-biology confounding.
- Protein filter (50%) applied per iteration after subsampling so the protein
  set honestly reflects availability at that cell count.
- Ridge regression for batch correction: fits on observed values only,
  imputes residual NaNs with per-protein corrected mean.
"""

# ── Imports ───────────────────────────────────────────────────────────────────
import os
import warnings
from pathlib import Path
from itertools import product
from dataclasses import dataclass

import numpy as np
import pandas as pd
import anndata
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from joblib import Parallel, delayed
from tqdm import tqdm
import scib_metrics
from scib_metrics.nearest_neighbors import pynndescent, NeighborsResults

# ── Paths ─────────────────────────────────────────────────────────────────────
NOTEBOOK_DIR = Path().resolve()
ROOT_OUT     = Path(os.environ.get("SCI_DATA_DIR",
                    (NOTEBOOK_DIR / "../../Data").resolve()))
FIGURES_OUT  = Path(os.environ.get("SCI_FIGURES_DIR",
                    (NOTEBOOK_DIR / "../../Figures").resolve()))

ANNDATA_IN  = ROOT_OUT / "AnnData" / "out"
RESULTS_DIR = ROOT_OUT / "batch_mixing"

for d in [
    RESULTS_DIR / "iterations",
    RESULTS_DIR / "reference",
    RESULTS_DIR / "rare_cell",
]:
    d.mkdir(parents=True, exist_ok=True)

# ── Global parameters ─────────────────────────────────────────────────────────
LAYER        = "X_log2"
COL_PLATE    = "plate_idx"
COL_CELLTYPE = "source name"
COL_OUTLIER  = "completeness_outlier"
CELLTYPE_A   = "HEK293T"    # dominant population
CELLTYPE_B   = "K562"       # rare population in Tier 2
CONDITIONS   = [CELLTYPE_A, CELLTYPE_B]

MIN_CELL_PCT = 0.50    # protein missingness filter, applied per iteration
N_PCS        = 10      # PCA components
RIDGE_ALPHA  = 1.0     # Ridge regularisation strength
N_ITERATIONS = 50      # bootstrap iterations per condition
N_JOBS       = 60      # joblib parallel workers
RANDOM_STATE = 42
KNN_K        = 30      # fixed k for all neighbourhood metrics

# Tier 1 grids
N_PLATES_LIST     = [2, 3, 4, 5, 6, 7, 8, 9]
N_CELLS_PER_PLATE = [10, 20, 30, 50, 75, 90]

# Tier 2 grids
RARE_FRACTIONS    = [0.02, 0.05, 0.10, 0.15, 0.20]
N_TOTAL_LIST      = [50, 100, 200, 500]
RARE_DISTRIBUTIONS = ["uniform", "random", "clustered_2plates"]
N_PLATES_RARE     = 9     # fix plates for Tier 2 to isolate confounding question
N_ITERATIONS_RARE = 50


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_data(path: Path) -> anndata.AnnData:
    """
    Load adata_protein, remove outlier cells and non-biological conditions.
    No protein-level filtering here — that happens per iteration.
    """
    adata    = anndata.read_h5ad(path)
    n_before = adata.n_obs

    adata = adata[~adata.obs[COL_OUTLIER].astype(bool)].copy()
    print(f"Outlier filter   : {n_before} -> {adata.n_obs} cells "
          f"({n_before - adata.n_obs} removed)")

    adata = adata[adata.obs[COL_CELLTYPE].isin(CONDITIONS)].copy()
    print(f"Condition filter : {adata.n_obs} cells kept ({CONDITIONS})")

    adata.X = adata.layers[LAYER].copy()

    print(f"\nLoaded: {adata.n_obs} cells x {adata.n_vars} proteins")
    print(f"Cell types : {adata.obs[COL_CELLTYPE].value_counts().to_dict()}")
    print(f"Plates     : {adata.obs[COL_PLATE].value_counts().to_dict()}\n")

    return adata


# ══════════════════════════════════════════════════════════════════════════════
# SUBSAMPLING — Tier 1 (balanced, n_cells per cell type per plate)
# ══════════════════════════════════════════════════════════════════════════════

def subsample(adata: anndata.AnnData,
              plates: list,
              n_cells: int,
              rng: np.random.Generator) -> anndata.AnnData | None:
    """
    Subsample n_cells of each cell type from each of the given plates.
    Returns None with a warning if any (plate, cell_type) pair is under-sampled.
    """
    idx_keep = []
    for plate in plates:
        for ct in CONDITIONS:
            mask       = ((adata.obs[COL_PLATE] == plate) &
                          (adata.obs[COL_CELLTYPE] == ct))
            candidates = np.where(mask)[0]
            if len(candidates) < n_cells:
                warnings.warn(
                    f"SKIP: plate={plate} celltype={ct} has "
                    f"{len(candidates)} cells, need {n_cells}"
                )
                return None
            idx_keep.extend(rng.choice(candidates, n_cells, replace=False))
    return adata[idx_keep].copy()


# ══════════════════════════════════════════════════════════════════════════════
# SUBSAMPLING — Tier 2 (rare cell simulation)
# ══════════════════════════════════════════════════════════════════════════════

def subsample_rare(
    adata: anndata.AnnData,
    n_total: int,
    rare_fraction: float,
    distribution: str,
    all_plates: list,
    rng: np.random.Generator,
) -> anndata.AnnData | None:
    """
    Build a dataset with a rare cell population (K562) embedded in a
    dominant population (HEK293T), with controlled distribution across plates.

    Parameters
    ----------
    n_total       : total number of cells
    rare_fraction : fraction of n_total that are rare (K562)
    distribution  : "uniform"         — rare cells spread evenly across all plates
                    "random"          — rare cells distributed randomly across plates
                    "clustered_2plates" — rare cells restricted to 2 plates only
    all_plates    : list of available plate names
    rng           : random generator for reproducibility

    Returns None if insufficient cells are available for the requested design.
    """
    n_rare     = max(1, int(n_total * rare_fraction))
    n_dominant = n_total - n_rare
    n_plates   = len(all_plates)

    # ── Dominant cells: draw evenly from all plates ───────────────────────────
    n_dominant_per_plate = n_dominant // n_plates
    dominant_idx = []
    for plate in all_plates:
        mask       = ((adata.obs[COL_PLATE] == plate) &
                      (adata.obs[COL_CELLTYPE] == CELLTYPE_A))
        candidates = np.where(mask)[0]
        if len(candidates) < n_dominant_per_plate:
            warnings.warn(
                f"SKIP: not enough dominant cells in plate={plate} "
                f"(have {len(candidates)}, need {n_dominant_per_plate})"
            )
            return None
        dominant_idx.extend(
            rng.choice(candidates, n_dominant_per_plate, replace=False)
        )

    # ── Rare cells: distribution determines which plates they come from ────────
    if distribution == "uniform":
        # Spread evenly across all plates
        n_rare_per_plate = n_rare // n_plates
        rare_plates      = all_plates
        rare_per_plate   = {p: n_rare_per_plate for p in rare_plates}
        # Assign remainder to first plate
        remainder = n_rare - n_rare_per_plate * n_plates
        if remainder > 0:
            rare_per_plate[all_plates[0]] += remainder

    elif distribution == "random":
        # Assign each rare cell to a randomly chosen plate
        plate_assignments = rng.choice(all_plates, size=n_rare, replace=True)
        rare_per_plate    = {
            p: int((plate_assignments == p).sum()) for p in all_plates
        }

    elif distribution == "clustered_2plates":
        # All rare cells come from exactly 2 plates
        rare_plates    = list(rng.choice(all_plates, size=2, replace=False))
        n_per          = n_rare // 2
        rare_per_plate = {p: 0 for p in all_plates}
        rare_per_plate[rare_plates[0]] = n_per
        rare_per_plate[rare_plates[1]] = n_rare - n_per

    else:
        raise ValueError(f"Unknown distribution: '{distribution}'. "
                         f"Choose from: uniform, random, clustered_2plates")

    rare_idx = []
    for plate, n_req in rare_per_plate.items():
        if n_req == 0:
            continue
        mask       = ((adata.obs[COL_PLATE] == plate) &
                      (adata.obs[COL_CELLTYPE] == CELLTYPE_B))
        candidates = np.where(mask)[0]
        if len(candidates) < n_req:
            warnings.warn(
                f"SKIP: not enough rare cells in plate={plate} "
                f"(have {len(candidates)}, need {n_req})"
            )
            return None
        rare_idx.extend(rng.choice(candidates, n_req, replace=False))

    all_idx = dominant_idx + rare_idx
    return adata[all_idx].copy()


# ══════════════════════════════════════════════════════════════════════════════
# PROTEIN FILTER
# ══════════════════════════════════════════════════════════════════════════════

def filter_proteins(adata_sub: anndata.AnnData,
                    min_cell_pct: float = MIN_CELL_PCT) -> anndata.AnnData:
    """
    Keep proteins present in >= min_cell_pct of subsampled cells.
    Applied fresh per iteration — protein set reflects actual data availability.
    """
    X    = np.array(adata_sub.X, dtype=float)
    keep = np.mean(~np.isnan(X), axis=0) >= min_cell_pct
    return adata_sub[:, keep].copy()


# ══════════════════════════════════════════════════════════════════════════════
# IMPUTATION (used for uncorrected baseline only)
# ══════════════════════════════════════════════════════════════════════════════

def impute_rowmean(X: np.ndarray) -> np.ndarray:
    """
    Per-protein mean imputation on raw matrix.
    Used to produce an uncorrected but complete matrix for baseline PCA.
    No batch correction applied.
    """
    X_imp = X.copy().astype(float)
    for j in range(X_imp.shape[1]):
        col      = X_imp[:, j]
        nan_mask = np.isnan(col)
        if nan_mask.any() and (~nan_mask).any():
            X_imp[nan_mask, j] = np.nanmean(col)
        elif nan_mask.all():
            X_imp[:, j] = 0.0
    return X_imp


# ══════════════════════════════════════════════════════════════════════════════
# BATCH CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

def remove_plate_effect(X: np.ndarray,
                        plate_labels: np.ndarray,
                        alpha: float = RIDGE_ALPHA) -> np.ndarray:
    """
    Regress out plate as a fixed effect using Ridge regression.
    For each protein:
      1. Fit Ridge on observed (non-NaN) cells only.
      2. Subtract plate contribution from observed cells.
      3. Impute residual NaNs with per-protein corrected mean.
    Ridge shrinks sparse plate coefficients toward zero rather than
    producing unstable OLS estimates.
    """
    enc    = OneHotEncoder(sparse_output=False, drop=None)
    P      = enc.fit_transform(plate_labels.reshape(-1, 1))
    X_corr = X.copy().astype(float)

    for j in range(X.shape[1]):
        col      = X_corr[:, j]
        obs_mask = ~np.isnan(col)
        if obs_mask.sum() < 2:
            X_corr[:, j] = np.nan
            continue
        reg = Ridge(alpha=alpha, fit_intercept=True)
        reg.fit(P[obs_mask], col[obs_mask])
        X_corr[obs_mask, j] = col[obs_mask] - (reg.predict(P[obs_mask]) - reg.intercept_)

    for j in range(X_corr.shape[1]):
        col      = X_corr[:, j]
        nan_mask = np.isnan(col)
        if nan_mask.any():
            col_mean            = np.nanmean(col)
            X_corr[nan_mask, j] = col_mean if not np.isnan(col_mean) else 0.0

    return X_corr


# ══════════════════════════════════════════════════════════════════════════════
# EMBEDDING
# ══════════════════════════════════════════════════════════════════════════════

def run_pca(X: np.ndarray,
            n_pcs: int = N_PCS,
            random_state: int = RANDOM_STATE) -> tuple[np.ndarray, np.ndarray]:
    """PCA — returns (scores [n_cells x n_pcs], variance_explained [n_pcs])."""
    n_pcs = min(n_pcs, X.shape[0] - 1, X.shape[1])
    pca   = PCA(n_components=n_pcs, random_state=random_state)
    return pca.fit_transform(X), pca.explained_variance_ratio_


def run_tsne(X: np.ndarray,
             use_pca: bool = True,
             n_pcs: int = 10,
             random_state: int = RANDOM_STATE) -> np.ndarray:
    """
    tSNE for visual illustration only — not used for quantitative metrics.
    If use_pca=True (default), reduces to PCA first for stability.
    """
    if use_pca and X.shape[1] > n_pcs:
        X = PCA(n_components=n_pcs,
                random_state=random_state).fit_transform(X)
    perplexity = min(30, max(5, X.shape[0] // 4))
    return TSNE(n_components=2, random_state=random_state,
                perplexity=perplexity).fit_transform(X)


# ══════════════════════════════════════════════════════════════════════════════
# METRICS
# ══════════════════════════════════════════════════════════════════════════════

def compute_neighbors(pcs: np.ndarray,
                      k: int = KNN_K) -> NeighborsResults:
    """
    Compute k nearest neighbours from PCA coordinates.
    k fixed across all conditions for comparability — not scaled with N.
    """
    k = min(k, pcs.shape[0] - 1)
    return pynndescent(pcs, n_neighbors=k)


def compute_ilisi(nn: NeighborsResults,
                  batch_labels: np.ndarray) -> float:
    """iLISI — batch mixing. [0,1], higher = better mixing."""
    try:
        return float(scib_metrics.ilisi_knn(nn, batch_labels))
    except Exception as e:
        warnings.warn(f"iLISI failed: {e}")
        return np.nan


def compute_clisi(nn: NeighborsResults,
                  celltype_labels: np.ndarray) -> float:
    """cLISI — cell-type separation. [0,1], higher = better separation."""
    try:
        return float(scib_metrics.clisi_knn(nn, celltype_labels))
    except Exception as e:
        warnings.warn(f"cLISI failed: {e}")
        return np.nan


@dataclass
class KBetResult:
    acceptance_rate : float
    stat_per_cell   : np.ndarray
    pvalue_per_cell : np.ndarray


def compute_kbet(nn: NeighborsResults,
                 batch_labels: np.ndarray,
                 alpha: float = 0.05) -> KBetResult:
    """
    kBET acceptance rate + per-cell chi-squared statistics and p-values.
    Higher acceptance rate = better batch mixing.
    """
    try:
        acceptance_rate, test_stats, p_values = scib_metrics.kbet(
            nn, batch_labels, alpha=alpha
        )
        return KBetResult(
            acceptance_rate = float(acceptance_rate),
            stat_per_cell   = np.array(test_stats, dtype=float),
            pvalue_per_cell = np.array(p_values,   dtype=float),
        )
    except Exception as e:
        warnings.warn(f"kBET failed: {e}")
        n = nn.n_samples if hasattr(nn, "n_samples") else len(batch_labels)
        return KBetResult(
            acceptance_rate = np.nan,
            stat_per_cell   = np.full(n, np.nan),
            pvalue_per_cell = np.full(n, np.nan),
        )


def per_cell_silhouette(pcs: np.ndarray,
                         labels: np.ndarray) -> np.ndarray:
    """Per-cell silhouette scores. Range [-1, 1]. Returns NaN array on failure."""
    try:
        return np.array(silhouette_samples(pcs, labels), dtype=float)
    except Exception as e:
        warnings.warn(f"Silhouette failed: {e}")
        return np.full(pcs.shape[0], np.nan)


def compute_clustering_metrics(pcs: np.ndarray,
                                 celltype_labels: np.ndarray,
                                 n_clusters: int = 2,
                                 random_state: int = RANDOM_STATE) -> dict:
    """
    Unsupervised k-means clustering on PCA, evaluated against true cell type
    labels using ARI and NMI. Answers: if we didn't know the labels, would
    unsupervised clustering recover the cell types?

    ARI : adjusted Rand index [-1, 1], higher = better label recovery
    NMI : normalised mutual information [0, 1], higher = better label recovery
    """
    try:
        km          = KMeans(n_clusters=n_clusters, random_state=random_state,
                             n_init=10)
        pred_labels = km.fit_predict(pcs)
        ari         = adjusted_rand_score(celltype_labels, pred_labels)
        nmi         = normalized_mutual_info_score(celltype_labels, pred_labels)
        return {"ari": float(ari), "nmi": float(nmi)}
    except Exception as e:
        warnings.warn(f"Clustering metrics failed: {e}")
        return {"ari": np.nan, "nmi": np.nan}


def compute_knn_purity(pcs: np.ndarray,
                        labels: np.ndarray,
                        target_label: str,
                        k: int = KNN_K) -> float:
    """
    kNN purity for a specific label — fraction of a target cell's k nearest
    neighbours that share the same label. Used for rare cell recovery.
    Returns mean purity across all cells with target_label.
    """
    try:
        from sklearn.neighbors import NearestNeighbors
        k_eff    = min(k, pcs.shape[0] - 1)
        nbrs     = NearestNeighbors(n_neighbors=k_eff).fit(pcs)
        _, idx   = nbrs.kneighbors(pcs)
        mask     = np.array(labels) == target_label
        if mask.sum() == 0:
            return np.nan
        purity   = np.array([
            np.mean(np.array(labels)[idx[i]] == target_label)
            for i in np.where(mask)[0]
        ])
        return float(np.nanmean(purity))
    except Exception as e:
        warnings.warn(f"kNN purity failed: {e}")
        return np.nan


def _compute_metrics_for_matrix(
    X: np.ndarray,
    plate_labels: np.ndarray,
    celltype_labels: np.ndarray,
    n_proteins: int,
    obs_names: list,
    condition_meta: dict,
    corrected: bool,
    store_per_cell: bool = True,
) -> tuple[dict, pd.DataFrame | None]:
    """
    Core metric computation for a given (cells x proteins) matrix.
    Computes PCA, neighbourhood graph, all scalar metrics, and optionally
    per-cell scores.

    Parameters
    ----------
    corrected      : whether this matrix has been batch-corrected
                     (used as a label in output, not for computation)
    store_per_cell : whether to return a per-cell DataFrame.
                     Set False for uncorrected baseline to save storage.
    """
    pcs, var_exp = run_pca(X)
    nn           = compute_neighbors(pcs)

    ilisi  = compute_ilisi(nn, plate_labels)
    clisi  = compute_clisi(nn, celltype_labels)
    kbet_r = compute_kbet(nn, plate_labels)
    clust  = compute_clustering_metrics(pcs, celltype_labels)

    sil_batch    = per_cell_silhouette(pcs, plate_labels)
    sil_celltype = per_cell_silhouette(pcs, celltype_labels)

    scalar = {
        "corrected"           : corrected,
        "ilisi"               : ilisi,
        "clisi"               : clisi,
        "kbet_acceptance"     : kbet_r.acceptance_rate,
        "kbet_stat_mean"      : float(np.nanmean(kbet_r.stat_per_cell)),
        "kbet_pvalue_mean"    : float(np.nanmean(kbet_r.pvalue_per_cell)),
        "silhouette_batch"    : float(np.nanmean(sil_batch)),
        "silhouette_celltype" : float(np.nanmean(sil_celltype)),
        "ari"                 : clust["ari"],
        "nmi"                 : clust["nmi"],
        "var_explained_pc1"   : float(var_exp[0]),
        "var_explained_top10" : float(var_exp.sum()),
        "n_cells"             : int(pcs.shape[0]),
        "n_proteins"          : int(n_proteins),
        **condition_meta,
    }

    per_cell_df = None
    if store_per_cell:
        per_cell_df = pd.DataFrame({
            "cell_index"          : obs_names,
            "plate"               : plate_labels,
            "celltype"            : celltype_labels,
            "silhouette_batch"    : sil_batch,
            "silhouette_celltype" : sil_celltype,
            "kbet_statistic"      : kbet_r.stat_per_cell,
            "kbet_pvalue"         : kbet_r.pvalue_per_cell,
            "corrected"           : corrected,
            **{k: v for k, v in condition_meta.items()},
        })

    return scalar, per_cell_df


def compute_all_metrics(
    X_raw: np.ndarray,
    X_corr: np.ndarray,
    plate_labels: np.ndarray,
    celltype_labels: np.ndarray,
    n_proteins: int,
    obs_names: list,
    condition_meta: dict,
) -> tuple[dict, pd.DataFrame]:
    """
    Compute metrics on both uncorrected (raw imputed) and corrected matrices.
    Returns a single scalar dict with delta columns and a per-cell DataFrame
    for the corrected matrix only.

    The delta columns (delta_ilisi, delta_clisi, etc.) quantify how much
    batch correction improved each metric — the key signal for the
    before/after comparison.
    """
    # Uncorrected baseline — no per-cell storage to keep parquet lean
    scalar_raw, _ = _compute_metrics_for_matrix(
        X_raw, plate_labels, celltype_labels,
        n_proteins, obs_names, condition_meta,
        corrected=False, store_per_cell=False,
    )

    # Corrected — store per-cell scores
    scalar_corr, per_cell_df = _compute_metrics_for_matrix(
        X_corr, plate_labels, celltype_labels,
        n_proteins, obs_names, condition_meta,
        corrected=True, store_per_cell=True,
    )

    # Delta metrics: corrected - uncorrected
    # Positive delta = improvement after correction
    delta_metrics = [
        "ilisi", "clisi", "kbet_acceptance",
        "silhouette_batch", "silhouette_celltype",
        "ari", "nmi",
    ]
    scalar_combined = {**scalar_corr}
    for m in delta_metrics:
        raw_val  = scalar_raw.get(m,  np.nan)
        corr_val = scalar_corr.get(m, np.nan)
        scalar_combined[f"raw_{m}"]   = raw_val
        scalar_combined[f"delta_{m}"] = (corr_val - raw_val
                                          if not (np.isnan(raw_val) or
                                                  np.isnan(corr_val))
                                          else np.nan)

    return scalar_combined, per_cell_df


# ══════════════════════════════════════════════════════════════════════════════
# REFERENCE RUN — Tier 1
# ══════════════════════════════════════════════════════════════════════════════

def run_reference(adata: anndata.AnnData) -> dict:
    """
    Full-data reference run using all available cells.
    Establishes the baseline embedding quality ceiling.
    """
    print("[reference] Running full-data embedding...")

    adata_filt   = filter_proteins(adata)
    X_raw        = np.array(adata_filt.X, dtype=float)
    X_raw_imp    = impute_rowmean(X_raw)
    X_corr       = remove_plate_effect(X_raw, adata_filt.obs[COL_PLATE].values)

    condition_meta = {
        "n_plates"          : int(adata.obs[COL_PLATE].nunique()),
        "n_cells_per_plate" : "all",
        "iteration"         : 0,
        "plates_sampled"    : ",".join(
            str(p) for p in sorted(adata.obs[COL_PLATE].unique())
        ),
    }

    scalar, per_cell_df = compute_all_metrics(
        X_raw        = X_raw_imp,
        X_corr       = X_corr,
        plate_labels = adata_filt.obs[COL_PLATE].values,
        celltype_labels = adata_filt.obs[COL_CELLTYPE].values,
        n_proteins   = adata_filt.n_vars,
        obs_names    = adata_filt.obs_names.tolist(),
        condition_meta = condition_meta,
    )

    pd.DataFrame([scalar]).to_csv(
        RESULTS_DIR / "reference" / "reference_metrics.csv", index=False
    )
    per_cell_df.to_parquet(
        RESULTS_DIR / "reference" / "reference_per_cell.parquet", index=False
    )

    print(f"  [corrected]   iLISI={scalar['ilisi']:.3f}  "
          f"cLISI={scalar['clisi']:.3f}  "
          f"kBET={scalar['kbet_acceptance']:.3f}  "
          f"ARI={scalar['ari']:.3f}")
    print(f"  [raw]         iLISI={scalar['raw_ilisi']:.3f}  "
          f"cLISI={scalar['raw_clisi']:.3f}  "
          f"kBET={scalar['raw_kbet_acceptance']:.3f}  "
          f"ARI={scalar['raw_ari']:.3f}")
    print(f"  [delta]       iLISI={scalar['delta_ilisi']:+.3f}  "
          f"cLISI={scalar['delta_clisi']:+.3f}  "
          f"ARI={scalar['delta_ari']:+.3f}\n")
    return scalar


# ══════════════════════════════════════════════════════════════════════════════
# SINGLE ITERATION WORKER — Tier 1
# ══════════════════════════════════════════════════════════════════════════════

def run_single_iteration(
    iter_idx  : int,
    adata     : anndata.AnnData,
    n_plates  : int,
    n_cells   : int,
    all_plates: list,
) -> tuple[dict, pd.DataFrame] | None:
    """
    One bootstrap iteration for a given (n_plates, n_cells) condition.
    Computes metrics on both uncorrected and corrected matrices.
    Stateless and safe for parallel execution via joblib.
    """
    rng = np.random.default_rng(
        seed=iter_idx * 10_000 + n_plates * 1_000 + n_cells
    )

    plates_sample = list(rng.choice(all_plates, size=n_plates, replace=False))
    adata_sub     = subsample(adata, plates_sample, n_cells, rng)
    if adata_sub is None:
        return None

    try:
        adata_filt = filter_proteins(adata_sub)
        if adata_filt.n_vars < 10:
            warnings.warn(
                f"iter={iter_idx} n_plates={n_plates} n_cells={n_cells}: "
                f"only {adata_filt.n_vars} proteins pass filter — skipping"
            )
            return None

        X_raw     = np.array(adata_filt.X, dtype=float)
        X_raw_imp = impute_rowmean(X_raw)
        X_corr    = remove_plate_effect(X_raw, adata_filt.obs[COL_PLATE].values)

    except Exception as e:
        warnings.warn(f"iter={iter_idx} failed: {e}")
        return None

    condition_meta = {
        "n_plates"          : n_plates,
        "n_cells_per_plate" : n_cells,
        "iteration"         : iter_idx,
        "plates_sampled"    : ",".join(str(p) for p in sorted(plates_sample)),
    }

    scalar, per_cell_df = compute_all_metrics(
        X_raw           = X_raw_imp,
        X_corr          = X_corr,
        plate_labels    = adata_filt.obs[COL_PLATE].values,
        celltype_labels = adata_filt.obs[COL_CELLTYPE].values,
        n_proteins      = adata_filt.n_vars,
        obs_names       = adata_filt.obs_names.tolist(),
        condition_meta  = condition_meta,
    )

    return scalar, per_cell_df


# ══════════════════════════════════════════════════════════════════════════════
# OUTER LOOP — Tier 1
# ══════════════════════════════════════════════════════════════════════════════

def run_all(
    adata         : anndata.AnnData,
    n_plates_list : list = N_PLATES_LIST,
    n_cells_list  : list = N_CELLS_PER_PLATE,
    n_iterations  : int  = N_ITERATIONS,
    n_jobs        : int  = N_JOBS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run all (n_plates x n_cells) conditions in parallel.

    Returns
    -------
    power_results : pd.DataFrame — scalar metrics per iteration (corrected +
                                   raw + delta columns)
    per_cell_all  : pd.DataFrame — per-cell scores for corrected matrix only
    """
    all_plates = sorted(adata.obs[COL_PLATE].unique().tolist())
    conditions = list(product(n_plates_list, n_cells_list))

    print(f"Available plates  : {all_plates}")
    print(f"Conditions        : {len(conditions)}  "
          f"({len(n_plates_list)} n_plates x {len(n_cells_list)} n_cells)")
    print(f"Iterations        : {n_iterations} per condition")
    print(f"Total worker calls: {len(conditions) * n_iterations}")
    print(f"Parallel workers  : {n_jobs}\n")

    all_scalar_rows   = []
    all_per_cell_rows = []

    for n_plates, n_cells in tqdm(conditions, desc="Conditions"):
        if n_plates > len(all_plates):
            warnings.warn(
                f"n_plates={n_plates} > available ({len(all_plates)}) — skipping"
            )
            continue

        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(run_single_iteration)(
                iter_idx, adata, n_plates, n_cells, all_plates
            )
            for iter_idx in range(n_iterations)
        )
        results = [r for r in results if r is not None]

        if results:
            scalars      = [r[0] for r in results]
            per_cell_dfs = [r[1] for r in results]

            df_scalar   = pd.DataFrame(scalars)
            df_per_cell = pd.concat(per_cell_dfs, ignore_index=True)

            key = f"nPlates{n_plates}_nCells{n_cells}"
            df_scalar.to_csv(
                RESULTS_DIR / "iterations" / f"{key}.csv", index=False
            )

            all_scalar_rows.extend(scalars)
            all_per_cell_rows.append(df_per_cell)

            print(f"  n_plates={n_plates:2d}  n_cells={n_cells:3d} | "
                  f"n={len(df_scalar):3d} | "
                  f"iLISI={df_scalar['ilisi'].median():.3f} "
                  f"(Δ{df_scalar['delta_ilisi'].median():+.3f})  "
                  f"ARI={df_scalar['ari'].median():.3f} "
                  f"(Δ{df_scalar['delta_ari'].median():+.3f})")

    power_results = pd.DataFrame(all_scalar_rows)
    per_cell_all  = (pd.concat(all_per_cell_rows, ignore_index=True)
                     if all_per_cell_rows else pd.DataFrame())

    return power_results, per_cell_all


# ══════════════════════════════════════════════════════════════════════════════
# SINGLE ITERATION WORKER — Tier 2 (rare cell simulation)
# ══════════════════════════════════════════════════════════════════════════════

def run_single_iteration_rare(
    iter_idx     : int,
    adata        : anndata.AnnData,
    n_total      : int,
    rare_fraction: float,
    distribution : str,
    all_plates   : list,
) -> dict | None:
    """
    One bootstrap iteration for the rare cell simulation.
    Returns scalar metrics only (no per-cell parquet for speed).
    Includes rare-cell-specific metrics: kNN purity and silhouette
    of the rare cluster specifically.
    """
    rng = np.random.default_rng(
        seed=iter_idx * 100_000
            + int(rare_fraction * 1000) * 100
            + n_total
    )

    adata_sub = subsample_rare(
        adata, n_total, rare_fraction, distribution, all_plates, rng
    )
    if adata_sub is None:
        return None

    try:
        adata_filt = filter_proteins(adata_sub)
        if adata_filt.n_vars < 10:
            return None

        X_raw     = np.array(adata_filt.X, dtype=float)
        X_raw_imp = impute_rowmean(X_raw)
        X_corr    = remove_plate_effect(X_raw, adata_filt.obs[COL_PLATE].values)

    except Exception as e:
        warnings.warn(f"rare iter={iter_idx} failed: {e}")
        return None

    ct_labels    = adata_filt.obs[COL_CELLTYPE].values
    plate_labels = adata_filt.obs[COL_PLATE].values
    n_rare_actual = int((ct_labels == CELLTYPE_B).sum())

    condition_meta = {
        "n_total"      : n_total,
        "rare_fraction": rare_fraction,
        "n_rare_actual": n_rare_actual,
        "distribution" : distribution,
        "iteration"    : iter_idx,
        "plates_sampled": ",".join(str(p) for p in sorted(all_plates)),
    }

    # Corrected metrics
    pcs_corr, var_exp = run_pca(X_corr)
    nn_corr           = compute_neighbors(pcs_corr)
    clust_corr        = compute_clustering_metrics(pcs_corr, ct_labels)
    sil_rare_corr     = float(np.nanmean(
        per_cell_silhouette(pcs_corr, ct_labels)[ct_labels == CELLTYPE_B]
    ))
    knn_purity_corr   = compute_knn_purity(pcs_corr, ct_labels, CELLTYPE_B)

    # Uncorrected baseline metrics
    pcs_raw, _    = run_pca(X_raw_imp)
    clust_raw     = compute_clustering_metrics(pcs_raw, ct_labels)
    sil_rare_raw  = float(np.nanmean(
        per_cell_silhouette(pcs_raw, ct_labels)[ct_labels == CELLTYPE_B]
    ))
    knn_purity_raw = compute_knn_purity(pcs_raw, ct_labels, CELLTYPE_B)

    return {
        # Condition metadata
        **condition_meta,
        "n_proteins"              : int(adata_filt.n_vars),
        # Corrected metrics
        "ari"                     : clust_corr["ari"],
        "nmi"                     : clust_corr["nmi"],
        "ilisi"                   : compute_ilisi(nn_corr, plate_labels),
        "clisi"                   : compute_clisi(nn_corr, ct_labels),
        "kbet_acceptance"         : compute_kbet(nn_corr, plate_labels).acceptance_rate,
        "silhouette_rare_cluster" : sil_rare_corr,
        "knn_purity_rare"         : knn_purity_corr,
        # Uncorrected baseline
        "raw_ari"                 : clust_raw["ari"],
        "raw_nmi"                 : clust_raw["nmi"],
        "raw_silhouette_rare"     : sil_rare_raw,
        "raw_knn_purity_rare"     : knn_purity_raw,
        # Delta (improvement from correction)
        "delta_ari"               : clust_corr["ari"]  - clust_raw["ari"],
        "delta_silhouette_rare"   : sil_rare_corr       - sil_rare_raw,
        "delta_knn_purity_rare"   : knn_purity_corr     - knn_purity_raw,
    }


# ══════════════════════════════════════════════════════════════════════════════
# OUTER LOOP — Tier 2 (rare cell simulation)
# ══════════════════════════════════════════════════════════════════════════════

def run_rare_cell_simulation(
    adata         : anndata.AnnData,
    rare_fractions: list = RARE_FRACTIONS,
    n_total_list  : list = N_TOTAL_LIST,
    distributions : list = RARE_DISTRIBUTIONS,
    n_iterations  : int  = N_ITERATIONS_RARE,
    n_jobs        : int  = N_JOBS,
) -> pd.DataFrame:
    """
    Run the rare cell simulation across all (rare_fraction x n_total x
    distribution) conditions.

    Returns a DataFrame with one row per iteration, including corrected,
    uncorrected, and delta metrics for rare cell recovery.
    """
    all_plates = sorted(adata.obs[COL_PLATE].unique().tolist())
    conditions = list(product(rare_fractions, n_total_list, distributions))

    print(f"Rare cell simulation")
    print(f"  Conditions        : {len(conditions)}  "
          f"({len(rare_fractions)} fractions x "
          f"{len(n_total_list)} n_total x "
          f"{len(distributions)} distributions)")
    print(f"  Iterations        : {n_iterations} per condition")
    print(f"  Total worker calls: {len(conditions) * n_iterations}")
    print(f"  Parallel workers  : {n_jobs}\n")

    all_rows = []

    for rare_frac, n_total, dist in tqdm(conditions, desc="Rare cell conditions"):
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(run_single_iteration_rare)(
                iter_idx, adata, n_total, rare_frac, dist, all_plates
            )
            for iter_idx in range(n_iterations)
        )
        results = [r for r in results if r is not None]

        if results:
            df_cond = pd.DataFrame(results)
            key     = (f"rare{int(rare_frac*100)}pct_"
                       f"n{n_total}_{dist}")
            df_cond.to_csv(
                RESULTS_DIR / "rare_cell" / f"{key}.csv", index=False
            )
            all_rows.extend(results)

            print(f"  rare={rare_frac:.0%}  n={n_total:4d}  dist={dist:20s} | "
                  f"n={len(df_cond):3d} | "
                  f"ARI={df_cond['ari'].median():.3f} "
                  f"(Δ{df_cond['delta_ari'].median():+.3f})  "
                  f"knn_purity={df_cond['knn_purity_rare'].median():.3f}")

    return pd.DataFrame(all_rows)


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    adata = load_data(ANNDATA_IN / "adata_protein.h5ad")

    # ── Tier 1: batch mixing power analysis ───────────────────────────────────
    print("=" * 60)
    print("TIER 1: Batch mixing power analysis")
    print("=" * 60)

    run_reference(adata)

    power_results, per_cell_all = run_all(adata)

    scalar_path  = RESULTS_DIR / "power_results_batch_mixing.csv"
    parquet_path = RESULTS_DIR / "per_cell_scores.parquet"

    power_results.to_csv(scalar_path, index=False)
    per_cell_all.to_parquet(parquet_path, index=False)

    print(f"\nTier 1 done.")
    print(f"  Scalar  : {len(power_results)} rows -> {scalar_path}")
    print(f"  Per-cell: {len(per_cell_all)} rows -> {parquet_path}")

    # ── Tier 2: rare cell simulation ──────────────────────────────────────────
    print("\n" + "=" * 60)
    print("TIER 2: Rare cell simulation")
    print("=" * 60)

    rare_results = run_rare_cell_simulation(adata)

    rare_path = RESULTS_DIR / "rare_cell" / "rare_cell_results.csv"
    rare_results.to_csv(rare_path, index=False)

    print(f"\nTier 2 done.")
    print(f"  Results : {len(rare_results)} rows -> {rare_path}")
    print("\nAll done.")