# SCI Notebooks

Analysis notebooks for the **Single-Cell Initiative (SCI)** inter-laboratory proteomics study. The project benchmarks single-cell proteomics across seven participating laboratories, nine 384-well plates, and two instrument platforms (Orbitrap Astral, timsTOF), using human cell lines (K562 and HEK293T as single cells, HeLa as a QC standard) as common reference material.

Notebooks are organised into six numbered directories that follow the logical progression of the analysis, from raw data ingestion to biological signal detection. Each directory contains its own `README.md` describing the individual notebooks in more detail.

---

## Repository structure

```
Notebooks/
├── 0_Create_data_objects/    Data ingestion and object construction
├── 1_QCQuant/                Search software benchmarking (QC mixtures)
├── 2_DatasetStatistics/      Identification depth, reproducibility, QC metrics, drift
├── 3_Correlations/           Sample correlations within and across labs
├── 4_VariancePartitioning/   Variance decomposition across experimental covariates
└── 5_DifferentialExpression/ Differential expression between cell lines
```

All notebooks read from and write to `Data/` and `Figures/` directories that sit alongside `Notebooks/` at the repository root. Both paths can be overridden via environment variables (`SCI_DATA_DIR`, `SCI_FIGURES_DIR`). The figures directory override is also how the Spectronaut figure set is produced from the same notebooks that default to DIA-NN.

---

## Analysis flow

The diagram below shows the data dependencies between notebook groups. Arrows indicate that one group produces outputs consumed by another.

```
SDRF metadata files
cellenONE isolation reports   ──►  0_Create_data_objects  ──►  Data/metadata/metadata.parquet
DIA-NN / Spectronaut results                               ──►  Data/AnnData/*.h5ad
                                          │
                                          ▼
                               2_DatasetStatistics
                               (identification counts, CVs, QC
                                metrics, drift over time;
                                writes completeness / outlier flags)
                                          │
                        ┌─────────────────┼─────────────────┐
                        ▼                 ▼                 ▼
              3_Correlations   4_VariancePartitioning   5_DifferentialExpression
              (sample-level    (variance decomposition, (K562 vs HEK293T,
               correlations)    batch correction)        cross-engine comparison)

   1_QCQuant runs on the quantitative QC mixtures, independent of the main
   single-cell chain; its results inform which search engine and settings are used.
```

Note that all results (AnnData objects and other outputs) are avaible on [Zenodo](10.5281/zenodo.21181080) as well. Search results for Spectronaut and the mixture datasets are available on PRIDE with accession number **PXD080610**.

---

## Directory summaries

### [`0_Create_data_objects/`](0_Create_data_objects/README.md) — Data object construction

Transforms raw search output and sample annotations into the analysis-ready data objects used throughout the project. Must be run first; all other directories depend on its outputs.

| Notebook | Description |
|---|---|
| `00_FASTA_overview.ipynb` | Summarises the mixture and human FASTA databases (per-species and contaminant counts) |
| `01_metadata.ipynb` | Merges SDRF files and cellenONE isolation reports into a single sample-level metadata table |
| `02_Results_to_AnnData_DIANN.ipynb` | Converts DIA-NN reports to precursor- and protein-level AnnData objects |
| `02_Results_to_AnnData_Spectronaut.ipynb` | Same conversion for Spectronaut output |

**Key outputs:** `Data/metadata/metadata.parquet`; `Data/AnnData/.../adata_precursor.h5ad`, `adata_protein.h5ad`, and their `_qc` variants, for each search engine.

---

### [`1_QCQuant/`](1_QCQuant/README.md) — Search software benchmarking

Benchmarks six DIA search engines on the quantitative QC mixtures (three-species samples with known ratios), covering identification depth, quantification accuracy, dynamic range, feature overlap, and precursor characteristics. Corresponds to Figure 2. Independent of the main single-cell chain; results inform the search engine and settings used elsewhere.

| Notebook | Description |
|---|---|
| `01_parse_results.ipynb` | Parses six search engines into ProteoBench format (`datapoints.pkl`, `intermediates.pkl`) |
| `02_Settings_comparison.ipynb` | Count-versus-error benchmark scatter figures across tools, settings, and instruments |
| `03_Performance_comparison.ipynb` | Dynamic range and quantification accuracy at the best setting per tool |
| `04_Overlap_between_search_engines.ipynb` | Precursor overlap (UpSet) and quantification-error panels per lab |
| `05_PrecursorCharacteristics.ipynb` | Peptide-level characteristics (length, charge, missed cleavages, modifications) |

> Notebooks `02` to `05` consume the parsed pickle objects, which can be regenerated by `01` or downloaded from Zenodo, and can therefore be run independently once those files are present. Additional runs with alternative settings/quantification workflows can be added as well.

> In the future, this analysis will be converted to a specific [ProteoBench](https://proteobench.cubimed.rub.de/) module to stimulate community contributions.

---

### [`2_DatasetStatistics/`](2_DatasetStatistics/README.md) — Dataset statistics and QC

Characterises identification depth, measurement reproducibility, QC metrics, and identification drift over time. Corresponds to Figures 3 and 4.

| Notebook | Description |
|---|---|
| `00_SearchEngineComparison.ipynb` | Five-way comparison of DIA-NN 2.2, DIA-NN 2.6, and Spectronaut (informative, not in manuscript) |
| `01_Counts_CVs.ipynb` | Identification counts, coefficients of variation, and missingness heatmaps |
| `02_QC_Metrics.ipynb` | Completeness curves, plate heatmaps, RT stability; writes completeness and outlier flags to AnnData |
| `03_Drift_over_time.ipynb` | Inter-day, long-term, and intra-day identification drift and RT displacement |

> `02_QC_Metrics.ipynb` writes the `completeness_outlier` sample flag back to the AnnData objects, which is used as a sample filter in the downstream differential expression analysis.

---

### [`3_Correlations/`](3_Correlations/README.md) — Sample correlations

Quantifies pairwise sample correlations from protein-level quantification, with missingness equalisation, to characterise intra-lab reproducibility, inter-lab similarity, and instrument-platform effects. Corresponds to Figure 5.

| Notebook | Description |
|---|---|
| `01_sample_correlations.ipynb` | Pairwise (Pearson and Spearman) correlations, diagnostics, and reproducibility testing |

---

### [`4_VariancePartitioning/`](4_VariancePartitioning/README.md) — Variance decomposition

Decomposes the variance in single-cell quantification into technical and biological contributions across four successive correction stages, at the precursor and protein level. Corresponds to Figure 6. See the folder README for the conda environment and Jupyter kernel setup.

| Notebook | Description |
|---|---|
| `01_scplainer.ipynb` | scplainer (`scp`) linear models on precursors: variance partitioning, GMF embeddings, silhouette scores |
| `02_VariancePartitioning.ipynb` | `variancePartition` mixed-effects models on proteins (R notebook) |

---

### [`5_DifferentialExpression/`](5_DifferentialExpression/README.md) — Differential expression

Tests for differential protein abundance between K562 and HEK293T single cells, and compares the results between search engines. Corresponds to Figure 6.

| Notebook | Description |
|---|---|
| `01_MSQROB.ipynb` | `msqrob2` mixed-effects DE models (integrated and per-plate), volcano and overlap figures (R notebook) |
| `02_DIANN_Spectronaut_overlap.ipynb` | DIA-NN versus Spectronaut DE comparison: overlap, logFC concordance, discordant proteins |

---

## Recommended execution order

```
0_Create_data_objects/01_metadata.ipynb
         ↓
0_Create_data_objects/02_Results_to_AnnData_{DIANN,Spectronaut}.ipynb
         ↓
2_DatasetStatistics/          ← run 02_QC_Metrics before downstream work
                                 (writes the outlier sample filter)
         ↓
3_Correlations/   4_VariancePartitioning/   5_DifferentialExpression/
```

> `1_QCQuant` is independent of this chain and only needs to be run when evaluating search software. It reads `metadata.parquet` from `0_Create_data_objects` and, for the benchmark figures, the parsed ProteoBench objects (regenerated locally or downloaded from Zenodo).
