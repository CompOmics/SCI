"""Sample-to-sample correlations with missingness equalised across labs.

Builds a reference feature set with equal missingness per lab, computes all
pairwise Spearman/Pearson correlations within each condition and compares
correlation distributions. Plots for these results live in
:mod:`SCI_py.plot.correlations`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from scipy.stats import mannwhitneyu
from tqdm import tqdm

from ..constants.labs import LABS

# ---------------------------------------------------------------------------
# Reference set construction
# ---------------------------------------------------------------------------

def create_reference_set(
    adata,
    conditions:              list[str],
    layer:                   str = '',
    min_sample_completeness: float = 0.30,
    min_protein_completeness: float = 0.50,
    random_seed:             int   = 42,
) -> tuple[dict, pd.DataFrame]:
    """
    Build a reference protein set with equalised missingness across labs.

    Three steps are applied per condition:

    1. **Sample filtering** — remove samples whose protein count falls
       below ``min_sample_completeness × median`` for their lab.
    2. **Reference protein selection** — keep proteins detected in at least
       ``min_protein_completeness`` of samples in *every* lab.
    3. **Missingness equalisation** — randomly mask observed values so that
       each lab's per-protein missingness matches the worst-case lab,
       ensuring fair cross-lab comparisons.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx`` and ``source name`` in ``obs``.
    conditions:
        Source names to process.
    layer:
        Expression matrix to use in ``.layers``. Defaults to ``X`` matrix.
    min_sample_completeness:
        Fraction of median protein count below which a sample is excluded.
    min_protein_completeness:
        Minimum per-lab detection fraction for a protein to be retained.
    random_seed:
        Seed for reproducible random masking.

    Returns
    -------
    masked_data:
        ``{condition: {lab: X_masked}}`` — reference matrices per
        condition and lab with equalised missingness.
    stats:
        DataFrame with masking statistics per condition and lab.
    """
    rng         = np.random.default_rng(random_seed)
    masked_data = {}
    stats       = []

    for cond in conditions:
        obs_cond = adata.obs[adata.obs['source name'] == cond].copy()

        # Step 1 — sample filter
        kept_samples = _filter_samples(adata, obs_cond, min_sample_completeness)
        obs_filtered = obs_cond.loc[kept_samples]

        # Step 2 — reference proteins
        lab_completeness = _compute_lab_completeness(adata, obs_filtered)
        if not lab_completeness:
            continue

        completeness_matrix = np.vstack(list(lab_completeness.values()))
        reference_idx       = np.where(
            np.all(completeness_matrix >= min_protein_completeness, axis=0)
        )[0]
        print(f'{cond}: {len(reference_idx)} reference proteins')

        # Step 3 — equalise missingness
        valid_labs          = [l for l in LABS if l in lab_completeness]
        missingness_matrix  = 1 - completeness_matrix[:, reference_idx]
        target_missingness  = np.max(missingness_matrix, axis=0)

        masked_data[cond] = {}
        for lab_i, lab in enumerate(valid_labs):
            lab_obs   = obs_filtered[obs_filtered['plate_idx'] == lab]
            idx       = [adata.obs.index.get_loc(s) for s in lab_obs.index]
            
            if layer == '':
                X_lab     = adata.X[np.ix_(idx, reference_idx)].copy()
            elif layer in adata.layers.keys():
                X_lab = adata.layers[layer][np.ix_(idx, reference_idx)].copy()
            else:
                raise Exception(f'{layer} does not exist in {adata.X.layers}. Leave empty to select X.')
                
            X_lab, n_masked = _equalise_missingness(
                X_lab, missingness_matrix[lab_i], target_missingness, rng
            )
            masked_data[cond][lab] = X_lab

            total_observed = np.sum(~np.isnan(X_lab))
            stats.append({
                'condition'              : cond,
                'lab'                    : lab,
                'n_reference_proteins'   : len(reference_idx),
                'n_samples'              : len(idx),
                'total_masked'           : sum(n_masked),
                'pct_masked'             : sum(n_masked) / (total_observed + sum(n_masked)) * 100,
                'mean_masked_per_protein': np.mean(n_masked),
            })

    return masked_data, pd.DataFrame(stats)


def equalise_plate_missingness(
    matrix:                    np.ndarray,
    plate_labels:              np.ndarray | pd.Series,
    min_sample_completeness:  float = 0.30,
    min_protein_completeness: float = 0.50,
    random_seed:               int   = 42,
) -> dict:
    """
    Equalise per-feature missingness across plates within a single sample
    group. This is the same three-step procedure as `create_reference_set`
    (sample filtering, reference-feature selection, missingness
    equalisation), but applied directly to a matrix / plate-label pair
    instead of stratifying by condition first — every sample passed in is
    treated as one group.

    Row order is always preserved: with `min_sample_completeness = 0` (no
    samples dropped), output rows line up 1:1 with `matrix`'s rows. When
    samples are dropped, `sample_mask` records which original rows survived,
    in their original order.

    Parameters
    ----------
    matrix:
        (n_samples x n_features) array. NaN marks a missing value.
    plate_labels:
        Length-n_samples plate/batch label per row, aligned with `matrix`.
    min_sample_completeness:
        Fraction of a plate's detectable features below which a sample is
        dropped. 0 disables sample filtering.
    min_protein_completeness:
        Minimum per-plate detection fraction required, in every plate, for a
        feature to be kept as a reference feature.
    random_seed:
        Seed for reproducible random masking.

    Returns
    -------
    dict
        - ``'masked'``: (n_kept_samples x n_reference_features) matrix, with
          extra values set to NaN so every plate has the same per-feature
          missingness (matching the plate with the most missing data for
          that feature)
        - ``'sample_mask'``: boolean mask into `matrix`'s rows (kept samples)
        - ``'feature_mask'``: boolean mask into `matrix`'s columns (reference features)
        - ``'stats'``: DataFrame of per-plate masking statistics
    """
    plate_labels = np.asarray(plate_labels)
    n_samples, n_features = matrix.shape
    plates = pd.unique(plate_labels)
    rng    = np.random.default_rng(random_seed)

    # --- Step 1: sample filtering -------------------------------------------
    sample_mask = np.ones(n_samples, dtype=bool)
    for plate in plates:
        plate_idx = np.where(plate_labels == plate)[0]
        X_plate   = matrix[plate_idx]

        detected_in_plate = np.sum(~np.isnan(X_plate), axis=0) > 0
        n_detectable       = detected_in_plate.sum()
        if n_detectable == 0:
            continue

        n_observed = np.sum(~np.isnan(X_plate[:, detected_in_plate]), axis=1)
        threshold  = min_sample_completeness * n_detectable
        sample_mask[plate_idx] = n_observed >= threshold

    kept_idx          = np.where(sample_mask)[0]
    kept_plate_labels = plate_labels[kept_idx]
    X_kept             = matrix[kept_idx]

    # --- Step 2: reference feature selection --------------------------------
    valid_plates = [p for p in plates if np.sum(kept_plate_labels == p) > 0]
    completeness = np.vstack([
        np.sum(~np.isnan(X_kept[kept_plate_labels == plate]), axis=0)
        / np.sum(kept_plate_labels == plate)
        for plate in valid_plates
    ])
    feature_mask           = np.all(completeness >= min_protein_completeness, axis=0)
    reference_completeness = completeness[:, feature_mask]

    # --- Step 3: missingness equalisation -----------------------------------
    missingness        = 1 - reference_completeness
    target_missingness = np.max(missingness, axis=0)

    X_ref = X_kept[:, feature_mask].copy()
    stats = []
    for i, plate in enumerate(valid_plates):
        plate_idx = np.where(kept_plate_labels == plate)[0]
        X_plate, n_masked = _equalise_missingness(
            X_ref[plate_idx], missingness[i], target_missingness, rng
        )
        X_ref[plate_idx] = X_plate

        total_masked   = int(sum(n_masked))
        total_observed = int(np.sum(~np.isnan(X_plate)))
        stats.append({
            'plate':                   plate,
            'n_reference_features':    int(feature_mask.sum()),
            'n_samples':               len(plate_idx),
            'total_masked':            total_masked,
            'pct_masked':              total_masked / (total_observed + total_masked) * 100
                                       if (total_observed + total_masked) > 0 else 0.0,
            'mean_masked_per_feature': np.mean(n_masked),
        })

    return {
        'masked':       X_ref,
        'sample_mask':  sample_mask,
        'feature_mask': feature_mask,
        'stats':        pd.DataFrame(stats),
    }


def compute_pairwise_correlations(
    masked_data:    dict,
    conditions:     list[str],
    adata,
    min_coobserved: int       = 50,
    methods:        list[str] = ('spearman', 'pearson'),
) -> dict[str, pd.DataFrame]:
    """
    Compute all pairwise sample correlations within each condition.

    For each condition, all lab matrices from ``masked_data`` are
    concatenated and every sample pair ``(i, j)`` is correlated over
    their co-observed proteins.  Pairs with fewer than ``min_coobserved``
    co-observed proteins are set to ``NaN``.

    Each pair is classified as:

    * ``intra_lab``                  — both samples from the same lab.
    * ``inter_lab_same_instrument``  — different labs, same instrument type.
    * ``inter_lab_diff_instrument``  — different labs, different instruments.

    Parameters
    ----------
    masked_data:
        Output of ``create_reference_set`` —
        ``{condition: {lab: X_masked}}``.
    conditions:
        Conditions to process.
    adata:
        AnnData with ``plate_idx`` and ``instrument_type`` in ``obs``.
    min_coobserved:
        Minimum number of co-observed proteins required to compute a
        correlation.  Pairs below this threshold are excluded.
    methods:
        Correlation methods to compute.  Any subset of
        ``['spearman', 'pearson']``.

    Returns
    -------
    dict mapping condition → DataFrame of pairwise correlation records.
    """
    results = {}

    for cond in conditions:
        cond_data = masked_data.get(cond, {})
        if not cond_data:
            continue

        X_all, meta_df = _build_condition_matrix(cond_data, adata)
        n_total        = X_all.shape[0]
        print(f'{cond}: {n_total} samples × {X_all.shape[1]} proteins '
              f'→ {n_total * (n_total - 1) // 2:,} pairs')

        X_ranked = _precompute_transforms(X_all, methods)
        records = _compute_pairs(
            X_all,
            X_ranked,
            meta_df,
            cond,
            methods,
            min_coobserved
        )

        results[cond] = pd.DataFrame(records)
        print(f'  Done: {len(records):,} valid pairs')

    return results


def compute_pairwise_stats(
    group_a: pd.Series,
    group_b: pd.Series,
    label_a: str,
    label_b: str,
    meaningful_delta: float = 0.02,  # minimum correlation difference worth reporting
) -> dict:
    """
    Compare two correlation distributions.
    Uses median difference as the primary effect size.
    Mann-Whitney p-value is reported but should not be used alone
    given large sample sizes — focus on delta and practical_sig instead.
    """
    a, b   = group_a.dropna(), group_b.dropna()
    _, p   = mannwhitneyu(a, b, alternative='two-sided')
    delta  = group_a.median() - group_b.median()

    # Common language effect size (probability that a random value from A > B)
    # Equivalent to the rank-biserial correlation
    n_a, n_b = len(a), len(b)
    u_stat, _ = mannwhitneyu(a, b, alternative='greater')
    cles = u_stat / (n_a * n_b)   # ranges 0-1, 0.5 = no effect

    return {
        'group_a'        : label_a,
        'group_b'        : label_b,
        'n_a'            : n_a,
        'n_b'            : n_b,
        'median_a'       : group_a.median(),
        'median_b'       : group_b.median(),
        'delta'          : delta,
        'cles'           : cles,           # effect size: >0.56 small, >0.64 medium, >0.71 large
        'p_value'        : p,
        'practical_sig'  : abs(delta) >= meaningful_delta,  # is the difference meaningful?
    }

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# For create_reference_set
def _filter_samples(adata, obs_cond, min_completeness):
    """Return samples with >= `min_completeness` x the proteins detectable in their lab."""
    kept = []
    for lab in LABS:
        lab_obs = obs_cond[obs_cond['plate_idx'] == lab]
        if lab_obs.empty:
            continue

        idx    = [adata.obs.index.get_loc(s) for s in lab_obs.index]
        X_lab  = adata.X[idx, :]

        # Only consider proteins observed in at least one sample in this lab
        detected_in_lab = np.sum(~np.isnan(X_lab), axis=0) > 0
        X_lab_filtered  = X_lab[:, detected_in_lab]

        # Threshold = min_completeness × number of proteins detectable in this lab
        n_detectable = detected_in_lab.sum()
        n_proteins   = np.sum(~np.isnan(X_lab_filtered), axis=1)
        threshold    = min_completeness * n_detectable

        kept.extend(lab_obs.index[n_proteins >= threshold])
    return kept


def _compute_lab_completeness(adata, obs_filtered):
    """Return per-lab arrays of the fraction of samples in which each protein is observed."""
    completeness = {}
    for lab in LABS:
        lab_obs = obs_filtered[obs_filtered['plate_idx'] == lab]
        if lab_obs.empty:
            continue
        idx = [adata.obs.index.get_loc(s) for s in lab_obs.index]
        X   = adata.X[idx, :]
        completeness[lab] = np.sum(~np.isnan(X), axis=0) / len(idx)
    return completeness


def _equalise_missingness(X, current_miss, target_miss, rng):
    """Randomly set observed values to NaN per column until it reaches `target_miss`; return X and counts masked."""
    n_masked = []
    for p in range(X.shape[1]):
        n_to_mask = int(round((target_miss[p] - current_miss[p]) * X.shape[0]))
        if n_to_mask > 0:
            observed  = np.where(~np.isnan(X[:, p]))[0]
            mask_idx  = rng.choice(observed, size=n_to_mask, replace=False)
            X[mask_idx, p] = np.nan
        n_masked.append(n_to_mask)
    return X, n_masked

# For compute_pairwise_correlations
def _build_condition_matrix(
    cond_data: dict,
    adata,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Concatenate per-lab matrices and build sample metadata."""
    X_list, meta = [], []

    for lab in LABS:
        if lab not in cond_data:
            continue
        X_lab      = cond_data[lab]
        instrument = (adata.obs[adata.obs['plate_idx'] == lab]['instrument_type'].iloc[0]
                      if 'instrument_type' in adata.obs.columns else 'unknown')

        meta.extend([{'lab': lab, 'instrument_type': instrument}] * X_lab.shape[0])
        X_list.append(X_lab)

    return np.concatenate(X_list, axis=0), pd.DataFrame(meta).reset_index(drop=True)


def _precompute_transforms(
    X: np.ndarray,
    methods: list[str],
) -> np.ndarray | None:
    """Precompute ranked matrix for Spearman correlation."""
    if 'spearman' not in methods:
        return None

    X_ranked = np.full_like(X, np.nan, dtype=float)
    for p in range(X.shape[1]):
        col      = X[:, p]
        observed = ~np.isnan(col)
        if observed.sum() > 1:
            X_ranked[observed, p] = rankdata(col[observed])

    return X_ranked


def _compute_pairs(
    X_all:    np.ndarray,
    X_ranked: np.ndarray | None,
    meta_df:  pd.DataFrame,
    cond:     str,
    methods:  list[str],
    min_coobserved: int,
) -> list[dict]:
    """Iterate over all sample pairs and compute correlations."""
    n_total = X_all.shape[0]
    records = []

    for i in tqdm(range(n_total), desc=cond):
        lab_i  = meta_df.loc[i, 'lab']
        inst_i = meta_df.loc[i, 'instrument_type']

        for j in range(i + 1, n_total):
            lab_j  = meta_df.loc[j, 'lab']
            inst_j = meta_df.loc[j, 'instrument_type']

            record = {
                'sample_i'       : i,
                'sample_j'       : j,
                'lab_i'          : lab_i,
                'lab_j'          : lab_j,
                'instrument_i'   : inst_i,
                'instrument_j'   : inst_j,
                'condition'      : cond,
                'same_lab'       : lab_i == lab_j,
                'same_instrument': inst_i == inst_j,
                'lab_pair'       : '__'.join(sorted([lab_i, lab_j])),
                'comparison_type': (
                    'intra_lab'                 if lab_i == lab_j  else
                    'inter_lab_same_instrument' if inst_i == inst_j else
                    'inter_lab_diff_instrument'
                ),
            }

            if 'spearman' in methods:
                obs       = ~np.isnan(X_ranked[i]) & ~np.isnan(X_ranked[j])
                n_coobs   = obs.sum()
                record['n_coobserved'] = n_coobs
                record['spearman']     = (
                    spearmanr(X_ranked[i, obs], X_ranked[j, obs])[0]
                    if n_coobs >= min_coobserved else np.nan
                )

            if 'pearson' in methods:
                obs_log          = ~np.isnan(X_all[i]) & ~np.isnan(X_all[j])
                n_coobs_log      = obs_log.sum()
                record['n_coobserved_pearson'] = n_coobs_log
                record['pearson'] = (
                    np.corrcoef(X_all[i, obs_log], X_all[j, obs_log])[0, 1]
                    if n_coobs_log >= min_coobserved else np.nan
                )

            if any(not np.isnan(record.get(m, np.nan)) for m in methods):
                records.append(record)

    return records