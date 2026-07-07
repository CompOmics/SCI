# 2_DatasetStatistics

Descriptive statistics and quality control for the single-cell SCI dataset across all participating laboratories. These notebooks operate on the AnnData objects produced in `0_Create_data_objects` and cover identification depth, reproducibility, QC metrics, and identification drift over time. Figures produced here correspond to Figures 3 and 4 of the manuscript. All notebooks accept the `SCI_DATA_DIR` and `SCI_FIGURES_DIR` environment variables to switch between search engines (DIA-NN by default, Spectronaut when overridden).

## Notebooks

### [00_SearchEngineComparison.ipynb](00_SearchEngineComparison.ipynb) (Informative)
Five-way comparison of DIA-NN 2.2, DIA-NN 2.6 (both first and second pass), and Spectronaut on the same samples. Covers precursor and protein identification counts, version-to-version improvement, feature overlap and subset coverage, per-sample quantitative correlation, overlap as a function of abundance rank, completeness curves, and CVs per lab. Marked informative but not included in the manuscript.

### [01_Counts_CVs.ipynb](01_Counts_CVs.ipynb)
Characterises three core quality dimensions per single cell: identification counts (precursors and proteins per sample, by lab and condition), coefficients of variation (within-condition reproducibility, `CV = std / mean x 100`), and the binary missingness landscape across the full feature space. Also reports the percentage of contaminant identifications per sample.

### [02_QC_Metrics.ipynb](02_QC_Metrics.ipynb)
Computes QC-metric figures: completeness curves for single-cell and QC injection conditions, per-lab plate heatmaps of per-well identification counts, and retention time stability. Also writes updated AnnData objects annotated with completeness.

### [03_Drift_over_time.ipynb](03_Drift_over_time.ipynb)
Characterises identification drift at three time scales: inter-day variability, long-term stability across the acquisition campaign, and intra-day drift by injection order (per-lab OLS regression with BH-corrected significance). Also covers retention time displacement of precursors over the run sequence.
