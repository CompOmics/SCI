# SCI Notebooks

Analysis notebooks for the **Single-Cell Interlaboratory (SCI)** proteomics study. The project benchmarks single-cell proteomics across nine participating laboratories and two instrument platforms (Orbitrap Astral, timsTOF), using three human cell lines (K562, HEK293T, HeLa) as a common reference material.

Notebooks are organised into six numbered directories that follow the logical progression of the analysis, from raw data ingestion to biological signal detection. Each directory contains its own `README.md` with detailed descriptions of individual notebooks, input/output tables, and figure references.

---

## Repository structure

```
Notebooks/
├── 0_Create_data_objects/    Data ingestion and object construction
├── 1_QCQuant/                Search software benchmarking and selection
├── 2_DatasetStatistics/      Identification depth, reproducibility, and QC metrics
├── 3_Correlations/           Sample correlations within and across labs
├── 4_VariancePartitioning/   Variance decomposition across experimental covariates
└── 5_DifferentialExpression/ Biological signal detection across conditions
```

All notebooks read from and write to `Data/` and `Figures/` directories that sit alongside `Notebooks/` at the repository root. Both paths can be overridden via environment variables (`SCI_DATA_DIR`, `SCI_FIGURES_DIR`).

---

## Analysis flow

The diagram below shows the data dependencies between notebook groups. Arrows indicate that one group produces outputs consumed by another.

```
Raw Spectronaut TSVs
SDRF metadata files          ──► 0_Create_data_objects  ──► Data/AnnData/*.h5ad
cellenONE isolation reports                              ──► Data/Parquet/*.parquet
                                                         ──► Data/metadata/metadata.parquet
                                          │
                                          ▼
                               1_QCQuant
                               (search software selection —
                                informs which Spectronaut
                                settings are used in 0_*)
                                          │
                                          ▼
                               2_DatasetStatistics
                               (identification counts, CVs,
                                plate heatmaps, RT stability,
                                completeness outlier flagging)
                                          │
                              ┌───────────┴────────────┐
                              ▼                        ▼
                    4_VariancePartitioning      3_Correlations
                    (variance decomposition,    (sample-level
                     batch correction)          correlations)
                              │                        ▲
                              │   batch-corrected data │
                              └────────────────────────┘
                                          │
                                          ▼
                               5_DifferentialExpression
                               (biological signal detection,
                                condition comparisons)
```

> **Note on the 3 ↔ 4 dependency:** The `3_Correlations` directory contains two classes of notebooks. Sample correlation analyses on the raw/normalized data can be run after `2_DatasetStatistics`. However, notebooks that assess correlations on **batch-corrected** data require output from `4_VariancePartitioning` first. Check the individual notebook dependencies in `3_Correlations/README.md` before running.

---

## Directory summaries

### `0_Create_data_objects/` — Data Object Construction

Transforms raw instrument output and metadata into the analysis-ready data objects used throughout the project. Must be run first; all other directories depend on its outputs.

| Notebook | Description |
|---|---|
| `01_metadata.ipynb` | Merges SDRF files and cellenONE isolation reports; applies an 8-step metadata parsing and imputation pipeline |
| `02_Results_to_AnnData.ipynb` | Converts Spectronaut quantification tables to AnnData objects at precursor and protein level |
| `03_AnnData_to_QFeatures.ipynb` | Transfers AnnData objects to R-side QFeatures format for Scplainer and VariancePartition *(to be written)* |
| `04_Assembly_AnnData.ipynb` | Consolidates analysis results from all downstream notebooks into a single annotated AnnData object *(to be written)* |

**Key outputs:** `Data/AnnData/adata_precursor.h5ad`, `adata_protein.h5ad`, `adata_precursor_qc.h5ad`, `adata_protein_qc.h5ad`; `Data/Parquet/search_results.parquet`; `Data/metadata/metadata.parquet`.

---

### `1_QCQuant/` — Search Software Benchmarking

Compares quantification software outputs to select the optimal search strategy for the main SCI dataset. Results from this directory inform the Spectronaut settings and report format used when producing the raw inputs for `0_Create_data_objects`.

> Run before (or alongside) `0_Create_data_objects` when re-processing the raw data from scratch. Not required when working from pre-computed AnnData objects.

---

### `2_DatasetStatistics/` — Dataset Statistics & QC

Characterises the identification depth, measurement reproducibility, spatial plate quality, and chromatographic stability of the dataset. Produces the QC metrics that gate sample inclusion in all downstream analyses.

| Notebook | Figure | Description |
|---|---|---|
| `01_Counts_CVs.ipynb` | Figure 2 | Identification counts (precursor & protein), CVs per condition, binary missingness heatmaps |
| `02_QC_Metrics.ipynb` | Figure 3 | Completeness curves, plate heatmaps, RT stability |

**Key outputs written back to AnnData:** `completeness_sample_SC/QC` (obs), `completeness_feature_{lab}_SC/QC` (var), `completeness_outlier` (obs).

> `02_QC_Metrics.ipynb` must be run before any analysis notebook, as `completeness_outlier` is used as a sample filter throughout.

---

### `3_Correlations/` — Sample Correlations

Quantifies pairwise and group-level sample correlations to characterise the sources of similarity and divergence in the dataset. Analyses are structured around three questions:

- How correlated are samples within the same lab?
- How correlated are samples across labs for the same condition?
- How does batch correction (from `4_VariancePartitioning`) improve cross-lab correlation?

**Dependency:** notebooks using batch-corrected data require output from `4_VariancePartitioning`. See `3_Correlations/README.md` for the per-notebook dependency map.

---

### `4_VariancePartitioning/` — Variance Decomposition

Decomposes the total variance in protein/precursor intensities into contributions from experimental covariates — lab, instrument type, cell line, gradient time, morphological features, and technical noise — using linear mixed models (Scplainer / VariancePartition). Produces batch-corrected expression matrices used by `3_Correlations` and `5_DifferentialExpression`.

**Dependency:** requires the QFeatures objects from `0_Create_data_objects/03_AnnData_to_QFeatures.ipynb`.

---

### `5_DifferentialExpression/` — Differential Expression

Tests for biological signal between conditions (K562 vs. HEK293T) before and after batch correction, quantifying the ability of the dataset to detect true biological differences in the presence of inter-laboratory technical variation. Results are placed in context of the variance components identified in `4_VariancePartitioning`.

**Dependency:** requires batch-corrected data from `4_VariancePartitioning`.

---

## Recommended execution order

```
0_Create_data_objects/01_metadata.ipynb
         ↓
0_Create_data_objects/02_Results_to_AnnData.ipynb
         ↓
2_DatasetStatistics/01_Counts_CVs.ipynb    (can run in parallel with 02)
2_DatasetStatistics/02_QC_Metrics.ipynb    ← gates all downstream work
         ↓
0_Create_data_objects/03_AnnData_to_QFeatures.ipynb
         ↓
4_VariancePartitioning/                    ← produces batch-corrected data
         ↓
3_Correlations/                            ← run raw-data notebooks first;
                                              batch-corrected notebooks after VP
         ↓
5_DifferentialExpression/
         ↓
0_Create_data_objects/04_Assembly_AnnData.ipynb
```

> `1_QCQuant` is independent of this chain and only needs to be re-run when evaluating a new search software or Spectronaut configuration.
