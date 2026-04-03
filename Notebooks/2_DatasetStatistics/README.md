# Dataset Statistics

This directory contains notebooks that characterise the quantitative quality and measurement properties of the SCI dataset across all nine participating laboratories. The analyses span identification depth, measurement reproducibility, spatial plate uniformity, and chromatographic stability. Results contribute to **Figure 2** and **Figure 3** of the manuscript.

**Prerequisite:** `0_Create_data_objects/02_Results_to_AnnData.ipynb` must be run before either notebook. `02_QC_Metrics.ipynb` writes completeness metrics back to the AnnData objects and should be run before any downstream analysis notebook.

---

## Overview

```
Data/AnnData/
  adata_precursor.h5ad  ──┬──► 01_Counts_CVs.ipynb ──► Figures/Figure_2/ (panels A, B, C)
  adata_protein.h5ad    ──┘
  adata_precursor_qc.h5ad  ──┐
  adata_protein_qc.h5ad    ──┤
                              └──► 02_QC_Metrics.ipynb ──► Figures/Figure_3/ (panels A–F + Supp.)
Data/Parquet/                                          ──► Data/AnnData/ (updated with QC columns)
  search_results.parquet  ──┘
```

---

## Notebooks

### `01_Counts_CVs.ipynb` — Identification Counts & Coefficients of Variation

**Figure 2.** Establishes the baseline analytical performance of the dataset by quantifying how many features are detected per cell and how reproducible those measurements are within each biological condition.

#### Analysis steps

**Identification counts (Panel A)**

The number of non-missing entries in `adata.X` is computed per sample and stored as `n_precursors` / `n_proteins` in `adata.obs`. Box plots are generated in two layouts for both precursor and protein level:

- *Combined* — all labs in a single panel, conditions on the x-axis. Provides a dataset-wide identification depth summary.
- *Split* — one panel per lab, x-labels suppressed (colours remain informative). Reveals inter-site variability and platform differences (Astral vs. timsTOF).

**Coefficients of variation (Panel B)**

Per-feature CVs (`std / mean × 100`) are computed independently for K562, HEK293T, QC - HeLa, and site-specific QC conditions and stored in `adata.uns['CVs']`. A minimum of 10 samples per condition is required for a feature to contribute to the distribution (`min_n = 10`). Box plots are generated for two proteome subsets:

- *Common proteome* — features passing `min_n` in all conditions simultaneously; enables a fair cross-condition CV comparison on the same feature set.
- *All features* — features passing `min_n` within each condition independently; maximises the number of features per condition but makes direct cross-condition comparisons less straightforward.

**Binary missingness heatmap (Panel C)**

Per-lab binary identification heatmaps (detected / not detected) across all precursors and proteins. Features are ordered by global detection frequency computed from K562 and HEK293T only (`sort_conditions`), fixing the feature axis across all lab figures to allow direct visual comparison. Samples are grouped by condition and sorted by descending identification count within each condition. K562 and HEK293T are capped at 50 samples per lab (`max_samples_per_condition = 50`); control and blank conditions are shown without subsampling.

#### Outputs (`Figures/Figure_2/`)

| Panel | File | Description |
|---|---|---|
| A | `A/counts_overview_precursors.png` | Precursor ID counts, all labs combined |
| A | `A/counts_overview_precursors_SPLIT.png` | Precursor ID counts, one panel per lab |
| A | `A/counts_overview_proteins.png` | Protein ID counts, all labs combined |
| A | `A/counts_overview_proteins_SPLIT.png` | Protein ID counts, one panel per lab |
| B | `B/CV_overview_proteins_common_box.png` | Protein CVs, common proteome |
| B | `B/CV_overview_proteins_all_box.png` | Protein CVs, per-condition proteome |
| B | `B/CV_overview_precursors_common_box.png` | Precursor CVs, common proteome |
| B | `B/CV_overview_precursors_all_box.png` | Precursor CVs, per-condition proteome |
| C | `C/{lab}_identification_heatmap.png` | Per-lab binary missingness heatmap (precursor & protein) |
| C | `C/legend_heatmap_conditions.png` | Condition colour legend |
| C | `C/legend_heatmap_identification.png` | Identification/missing colour legend |

---

### `02_QC_Metrics.ipynb` — QC Metrics

**Figure 3.** Evaluates three complementary quality dimensions spanning the quantitative, spatial, and chromatographic properties of the data. Completeness columns computed here are written back to the AnnData objects and used as sample-level QC filters in all downstream notebooks.

#### Analysis steps

**Missingness / completeness curves (Panels E, F)**

Completeness curves visualise the trade-off between missing-value filter stringency and data retention. Two perspectives are computed per lab over thresholds 0–100%:

- *Feature-level completeness* — fraction of features detected in ≥ *t*% of samples. Reveals whether missingness is diffuse or concentrated in rarely detected features.
- *Sample-level completeness* — fraction of samples with ≥ *t*% of the lab's detectable features. Identifies near-empty injections from failed electrospray, blocked columns, or dispensing errors.

Analysis is run separately for **single-cell conditions** (K562, HEK293T, suffix `_SC`) and **QC injection conditions** (QC - HeLa, site-specific QC, suffix `_QC`). Because QC injections use bulk material, their completeness is expected to be substantially higher; a low-completeness QC injection is therefore a stronger indicator of an instrument or LC problem than the same value in a single cell.

Per-sample completeness (`completeness_sample_SC/QC`) and per-feature detection rates (`completeness_feature_{lab}_SC/QC`) are stored in `adata.obs` and `adata.var` respectively via `store_completeness`.

**Completeness outlier flagging**

Samples are labelled `completeness_outlier = True` if they fail both the SC and QC completeness thresholds (< 25%). The AND logic makes the filter robust to the inherently lower completeness of single-cell data: a sample must be low in both contexts to be excluded. The 25% threshold sits in the gap between the main population of well-measured samples and the cluster of clearly failed injections, as read directly from the completeness curves. The flag is propagated from the precursor to the protein AnnData by index alignment.

**Plate heatmaps (Panels C, D)**

Per-well identification counts are shown on the physical 16 × 24 plate grid. Each well is coloured by condition with opacity scaled linearly from 0.3 (minimum count) to 1.0 (maximum count within the plate). White boundary lines mark the distinct SCI plate sections (cell wells, QC wells, blank wells). QC injection wells are overlaid on the same figure via `extra_adata`, enabling direct comparison of cell and QC well coverage within a single plate — a pattern confined to cell wells points to a dispensing problem, while one spanning QC wells indicates an instrument or LC issue. One figure is produced per lab, at both protein and precursor level.

**Retention time stability (Panels A, B + Supplementary)**

The elution profile of each sample is represented as the proportion of identified precursors per RT bin (25 bins), with RT normalised to [0, 1] per sample to remove the effect of absolute gradient length differences between labs. Only samples passing the completeness outlier filter are included. Two figure types are generated:

- *Per-lab supplementary figures* — individual sample curves (transparent) overlaid with per-condition group means. Reveal whether instability is condition-specific or lab-wide.
- *Overview figures* — group mean curves only, coloured by lab, generated separately for single-cell (K562, HEK293T) and QC conditions.

#### Outputs (`Figures/Figure_3/`)

| Panel | File | Description |
|---|---|---|
| A | `A/RT_stability_overview_SC.png` | RT stability overview, single-cell conditions, coloured by lab |
| B | `B/RT_stability_overview_QC.png` | RT stability overview, QC conditions, coloured by lab |
| C | `C/plate_heatmap_protein_{lab}.png` | Per-lab protein ID plate heatmap |
| D | `D/plate_heatmap_precursor_{lab}.png` | Per-lab precursor ID plate heatmap |
| E | `E/completeness_precursors_SC.png` | Precursor completeness curves, single cells |
| E | `E/completeness_precursors_QC.png` | Precursor completeness curves, QC injections |
| F | `F/completeness_proteins_SC.png` | Protein completeness curves, single cells |
| F | `F/completeness_proteins_QC.png` | Protein completeness curves, QC injections |
| Supp. | `Supplementary/RT_stability_{lab}.png` | Per-lab RT stability coloured by condition |

#### AnnData updates

The following columns are added and persisted before downstream use:

| Object | Slot | Column | Description |
|---|---|---|---|
| precursor / protein | `obs` | `completeness_sample_SC` | Fraction of lab-detected features found in each single cell |
| precursor / protein | `obs` | `completeness_sample_QC` | Fraction of lab-detected features found in each QC injection |
| precursor / protein | `var` | `completeness_feature_{lab}_SC` | Per-lab single-cell detection rate for each feature |
| precursor / protein | `var` | `completeness_feature_{lab}_QC` | Per-lab QC detection rate for each feature |
| precursor / protein | `obs` | `completeness_outlier` | `True` if sample fails both SC and QC completeness thresholds (< 25%) |

---

## Execution Order

```
01_Counts_CVs.ipynb        ← independent; no AnnData modifications
       ↓
02_QC_Metrics.ipynb        ← writes completeness columns back to AnnData
                              → required before downstream analysis notebooks
```

> `01_Counts_CVs.ipynb` can be run independently of `02_QC_Metrics.ipynb`. However, `02_QC_Metrics.ipynb` must be run before proceeding to the analysis notebooks in `3_Correlations`, `4_VariancePartitioning`, and `5_DifferentialExpression`, as those notebooks filter samples on `completeness_outlier`.
