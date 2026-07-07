# 3_Correlations

Sample-to-sample reproducibility analysis based on pairwise correlations of protein-level quantification. Figures produced here correspond to Figure 5 of the manuscript. The notebook accepts the `SCI_DATA_DIR` and `SCI_FIGURES_DIR` environment variables to switch between search engines.

## Notebooks

### [01_sample_correlations.ipynb](01_sample_correlations.ipynb)
Computes pairwise sample correlations (Pearson and Spearman) from protein-level intensities, with a missingness-equalisation step so that correlations are compared on a consistent set of co-observed features. Includes diagnostics on the reference-set construction, then reports intra-lab reproducibility (violin plots per lab), inter-lab similarity (clustered heatmaps of median correlations), instrument-platform effects, example sample-pair scatter plots, and statistical testing of intra-lab reproducibility.
