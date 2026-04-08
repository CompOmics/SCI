# Single-Cell Interlaboratory (SCI) Proteomics

Analysis code for the SCI benchmarking study, which evaluates single-cell proteomics across **nine participating laboratories** and **two instrument platforms** (Orbitrap Astral, timsTOF) using two human cell lines (K562, HEK293T) as common reference materials.

The repository contains:

- **`SCI_py/`** — Python package with all metrics, plotting, and utility functions used by the analysis notebooks.
- **`SCI_R/`** — R scripts for statistical modelling (Scplainer, VariancePartition).
- **`Notebooks/`** — Jupyter notebooks organised by analysis stage, from raw data ingestion to biological signal detection.

---

**NOTE: This repo is still under construction and will thus still undergo a lot of changes.**

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/)
- Git

The conda environment installs **Python 3.13** and **R 4.5** together with all required packages (see `environment.yml`). No separate R installation is needed.

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/CompOmics/SCI.git
cd SCI
```

### 2. Create the conda environment
```bash
conda env create -f environment.yml
conda activate sci_environment
```

### 3. Install the SCI_py package
```bash
pip install -e ./SCI_py
```

The `-e` flag installs the package in editable mode — any changes to the source files under `SCI_py/` are immediately reflected without reinstalling.

---

## Verifying the installation
```python
from SCI_py.metrics import correlations, drift
from SCI_py.plot import identifications
import anndata
print("Installation successful.")
```

---

## Data

Pre-processed AnnData objects and supporting files are distributed as `Data.zip` alongside this repository. Extract the archive at the repository root:

```bash
unzip Data.zip
```

This creates a `Data/` directory with the following layout expected by the notebooks:

```
Data/
├── AnnData/
│   ├── adata_precursor.h5ad      # Precursor-level quantification, single cells
│   ├── adata_protein.h5ad        # Protein-level quantification, single cells
│   ├── adata_precursor_qc.h5ad   # Precursor-level quantification, QC injections
│   └── adata_protein_qc.h5ad     # Protein-level quantification, QC injections
├── Parquet/
│   ├── search_results.parquet    # Long-format search results, single cells
│   └── search_results_qc.parquet # Long-format search results, QC injections
└── metadata/
    └── metadata.parquet          # Merged sample-level metadata
```

To use a different data location (e.g. when running on a cluster or in CI), set environment variables before launching JupyterLab:

```bash
export SCI_DATA_DIR=/path/to/your/Data
export SCI_FIGURES_DIR=/path/to/your/Figures
```

All notebooks read these variables at startup and fall back to the default `Data/` and `Figures/` paths if they are not set.

---

## Repository structure

```
SCI/
├── environment.yml               Conda environment specification
├── Data/                         Analysis data (see above; not tracked by git)
├── Figures/                      Figure outputs (not tracked by git)
├── SCI_py/                       Python package
│   ├── setup.py
│   └── SCI_py/
│       ├── constants/            Lab maps, condition colours, plotting constants
│       ├── metrics/              Computation modules
│       │   ├── completeness.py   Completeness curve computation and storage
│       │   ├── correlations.py   Reference set construction, pairwise correlations
│       │   ├── cvs.py            Coefficient of variation matrices
│       │   ├── drift.py          Intra-/inter-day drift quantification
│       │   └── rt_stability.py   Retention time elution profile computation
│       └── plot/                 Visualisation modules (one per metrics module)
│           ├── correlations.py
│           ├── drift.py
│           ├── identifications.py
│           ├── plate_heatmap.py
│           └── rt_stability.py
├── SCI_R/                        R scripts for Scplainer / VariancePartition
└── Notebooks/
    ├── 0_Create_data_objects/    Data ingestion and AnnData / QFeatures construction
    ├── 1_QCQuant/                Search software benchmarking and selection
    ├── 2_DatasetStatistics/      Identification depth, reproducibility, QC metrics
    ├── 3_Correlations/           Sample correlations within and across labs
    ├── 4_VariancePartitioning/   Variance decomposition and batch correction
    └── 5_DifferentialExpression/ Biological signal detection across conditions
```

---

## Notebooks overview

Notebooks are numbered to reflect the recommended execution order. Each directory contains its own `README.md` with input/output tables and figure references.

### `0_Create_data_objects/` — Data Object Construction

Transforms raw Spectronaut output and sample metadata into the AnnData and QFeatures objects consumed by all downstream notebooks. **Must be run first.**

| Notebook | Description |
|---|---|
| `01_metadata.ipynb` | Merges SDRF files and cellenONE isolation reports into a unified metadata table |
| `02_Results_to_AnnData.ipynb` | Converts Spectronaut quantification tables to AnnData objects at precursor and protein level |
| `03_AnnData_to_QFeatures.ipynb` | Transfers AnnData objects to R-side QFeatures for Scplainer and VariancePartition |
| `04_Assembly_AnnData.ipynb` | Consolidates all analysis results into a single annotated AnnData object for release |

### `1_QCQuant/` — Search Software Benchmarking

Compares quantification software outputs to select the optimal Spectronaut configuration. Only needs re-running when evaluating a new search strategy.

### `2_DatasetStatistics/` — Dataset Statistics & QC

Characterises identification depth, CV, spatial plate quality, chromatographic stability, and sample completeness. The completeness outlier flags written by `02_QC_Metrics.ipynb` are required by all downstream notebooks.

| Notebook | Figure | Description |
|---|---|---|
| `01_Counts_CVs.ipynb` | Figure 3 | Identification counts, CVs, binary missingness heatmaps |
| `02_QC_Metrics.ipynb` | Figure 4 | Completeness curves, plate heatmaps, RT stability |
| `03_Drift_over_time.ipynb` | Figure 4 | Inter-day variability, long-term stability, intra-day drift |

### `3_Correlations/` — Sample Correlations

Pairwise and group-level sample correlations to characterise the sources of similarity and divergence in the dataset, before and after batch correction.

| Notebook | Description |
|---|---|
| `01_sample_correlations.ipynb` | Intra- and inter-lab correlations on the reference proteome with equalised missingness |
| `02_correlation_after_bec.ipynb` | Repeat correlation analysis on batch-corrected data (requires `4_VariancePartitioning`) |

### `4_VariancePartitioning/` — Variance Decomposition

Decomposes protein/precursor intensity variance into contributions from lab, instrument type, cell line, gradient time, morphological features, and noise using linear mixed models. Produces batch-corrected expression matrices for downstream use.

### `5_DifferentialExpression/` — Differential Expression

Tests for biological signal between K562 and HEK293T before and after batch correction, quantifying the dataset's ability to detect true biological differences in the presence of inter-laboratory technical variation.

---

## Recommended execution order

```
0_Create_data_objects/01_metadata.ipynb
        ↓
0_Create_data_objects/02_Results_to_AnnData.ipynb
        ↓
2_DatasetStatistics/01_Counts_CVs.ipynb   ─┐
2_DatasetStatistics/02_QC_Metrics.ipynb   ─┤ (can run in parallel)
2_DatasetStatistics/03_Drift_over_time.ipynb ─┘
        ↓
0_Create_data_objects/03_AnnData_to_QFeatures.ipynb
        ↓
4_VariancePartitioning/
        ↓
3_Correlations/                           ← raw-data notebooks run after step 2;
                                             batch-corrected notebooks run after VP
        ↓
5_DifferentialExpression/
        ↓
0_Create_data_objects/04_Assembly_AnnData.ipynb
```

> `1_QCQuant` is independent of this chain and only needs to be re-run when evaluating a new Spectronaut configuration.
>
> `02_QC_Metrics.ipynb` must complete before any analysis notebook, as it writes the `completeness_outlier` sample filter used throughout.
