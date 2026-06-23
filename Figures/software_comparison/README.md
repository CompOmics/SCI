# Software Comparison Figures

Cross-search-engine comparison figures. Files use an `{section}{number}_{description}` naming
convention (e.g. `A1_precursor_counts_overview`).

All files are saved in **png + svg + pdf**.

---

## Overview of sections

| Section | Comparison | Generating notebook |
|---------|-----------|---------------------|
| [A — DIA-NN 2.2 vs 2.6](#section-a--diann-22-vs-diann-26) | Version improvement | `2_DatasetStatistics/00_DIANN_version_comparison.ipynb` |
| [B — DIA-NN 2.6 vs Spectronaut](#section-b--diann-26-vs-spectronaut) | Software comparison | `2_DatasetStatistics/00_DIANN_version_comparison.ipynb` |
| [C — Completeness and CVs](#section-c--completeness-and-cvs) | Quality metrics across datasets | `2_DatasetStatistics/00_DIANN_version_comparison.ipynb` |
| [D — DE result overlap](#section-d--differential-expression-overlap) | MSqRob agreement | `5_DifferentialExpression/05_DIANN_Spectronaut_overlap.ipynb` |

---

## Section A — DIA-NN 2.2 vs DIA-NN 2.6

**Notebook:** [`2_DatasetStatistics/00_DIANN_version_comparison.ipynb`](../../Notebooks/2_DatasetStatistics/00_DIANN_version_comparison.ipynb)

Datasets compared: DIA-NN 2.2, DIA-NN 2.6 first pass, DIA-NN 2.6 second pass.

| File | Description |
|------|-------------|
| `A1_precursor_counts_overview` | Grouped boxplot of precursor counts per condition across all DIA-NN dataset versions |
| `A1_precursor_counts_per_lab` | Lab-by-lab faceted boxplot of precursor counts |
| `A2_protein_counts` | Two panels: protein IDs per condition (left) and per lab (right) for 2.2 vs 2.6 |
| `A3_version_improvement` | Paired per-sample scatter + per-lab gain boxplots for 2.2 → 2.6 (precursor 1st pass, precursor 2nd pass, protein 2nd pass) |
| `A4_precursors_overlap` | Feature overlap bars + per-sample coverage histograms for DIA-NN 2.2 vs 2.6 |
| `A5_correlation_precursor` | Per-sample Pearson r histograms of log2 precursor intensities (2.2 vs 2.6) |
| `A5_correlation_protein` | Per-sample Pearson r histograms of log2 protein intensities (2.2 vs 2.6) |

---

## Section B — DIA-NN 2.6 vs Spectronaut

**Notebook:** [`2_DatasetStatistics/00_DIANN_version_comparison.ipynb`](../../Notebooks/2_DatasetStatistics/00_DIANN_version_comparison.ipynb)

| File | Description |
|------|-------------|
| `B1_precursors_counts_per_lab` | Per-lab faceted boxplot of precursor ID counts (DIA-NN 2.6 vs Spectronaut) |
| `B1_proteins_counts_per_lab` | Same for protein ID counts |
| `B2_precursors_overlap` | Feature overlap bars + directional coverage histograms at precursor level |
| `B2_proteins_overlap` | Same at protein level |
| `B3_correlation_precursor` | Per-sample Pearson r histograms for precursor intensities |
| `B3_correlation_protein` | Per-sample Pearson r histograms for protein intensities |
| `B4_paired_per_lab` | Paired per-sample ID scatter per lab (precursor and protein rows) |
| `B5_overlap_per_sample` | Per-sample stacked precursor overlap bars (only DIA-NN / shared / only Spectronaut), faceted by lab |
| `B5_overlap_overview` | Mean per-sample precursor overlap: all samples and per-lab summary |
| `B6_binned_overlap_boxplot_precursors` | Precursor overlap % per abundance rank bin (boxplot) |
| `B6_binned_overlap_boxplot_proteins` | Same for proteins |
| `B6_binned_overlap_profile_precursors` | Per-rank abundance profile: example stacked bar + all-sample mean±SD (precursors) |
| `B6_binned_overlap_profile_proteins` | Same for proteins |

---

## Section C — Completeness and CVs

**Notebook:** [`2_DatasetStatistics/00_DIANN_version_comparison.ipynb`](../../Notebooks/2_DatasetStatistics/00_DIANN_version_comparison.ipynb)

| File | Description |
|------|-------------|
| `C1_completeness_per_lab` | Completeness curves per lab, one line per dataset (all five: DIA-NN 2.2, 2.6 1st/2nd pass, Spectronaut) |
| `C2_cv_per_lab` | CV per lab across datasets |

---

## Section D — Differential Expression Overlap

**Notebook:** [`5_DifferentialExpression/05_DIANN_Spectronaut_overlap.ipynb`](../../Notebooks/5_DifferentialExpression/05_DIANN_Spectronaut_overlap.ipynb)

Compares MSqRob DE results (K562 vs HEK293T) from DIA-NN and Spectronaut. Models: integrated
(all labs) and per-plate (nine plates: Astral_1–6, timsTOF_1–3).

| File | Description |
|------|-------------|
| `D1_proteins_per_plate` | Grouped bar chart of proteins quantified per plate (DIA-NN vs Spectronaut) |
| `D2_protein_overlap_upset` | UpSet plot of protein detection overlap across four result sets (DIA-NN integrated/per-plate, Spectronaut integrated/per-plate) |
| `D3_logFC_scatter_integrated` | logFC concordance scatter (integrated models), coloured by significance group |
| `D4_significance_overlap_bar` | Bar chart of significant protein counts by overlap group (Both / DIA-NN only / Spectronaut only / Neither) |
| `D4_significance_overlap_upset` | UpSet plot of significant protein overlap (FDR < 0.01, \|log2FC\| > 0.5) |
| `D5_logFC_per_plate_grid` | 3×3 grid of per-plate logFC scatter plots, coloured by significance group |
| `D5_pearson_r_per_plate` | Bar chart of per-plate Pearson r, coloured by instrument type (Astral vs timsTOF) |
| `D6_outlier_scatter` | Integrated model logFC scatter with opposite-direction outlier proteins highlighted |
| `D7_outlier_perplate_heatmap` | Side-by-side heatmaps of logFC for outlier proteins across all plates |
| `D9_volcano_per_plate` | Per-plate panel of volcano plots: DIA-NN / Spectronaut / significant protein count bars (27 subplots total) |

---

[← Back to Figures overview](../README.md)
