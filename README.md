# Single-Cell Interlaboratory (SCI) Proteomics

Analysis code for the SCI benchmarking study, which evaluates single-cell proteomics across **seven participating laboratories**, **nine 384-well plates**, and **two instrument platforms** (Orbitrap Astral, timsTOF). Two human cell lines (K562, HEK293T) are used as single-cell reference material and HeLa as a QC standard. Search results are analysed for two DIA search engines (DIA-NN and Spectronaut).

The repository contains:

- **`SCI_py/`** — Python package with the metrics, plotting, parsing, and utility functions used by the analysis notebooks.
- **`SCI_R/`** — R plotting and helper scripts used by the R notebooks (MSqRob, scplainer, variancePartition).
- **`Notebooks/`** — Jupyter notebooks organised by analysis stage, from raw data ingestion to differential expression. See [`Notebooks/README.md`](Notebooks/README.md) for the full walkthrough and recommended execution order.

The code is released under the Apache License 2.0 (see [`LICENSE`](LICENSE)).

---

## 1. System requirements

### Software dependencies

All dependencies are managed through conda. The repository provides four environment specifications, each covering a distinct part of the analysis:

| Environment file | Env name | Interpreter | Used by |
|---|---|---|---|
| [`environment.yml`](environment.yml) | `sci_environment` | Python 3.13, R 4.5 (base) | All Python notebooks (folders 0, 1, 2, 3, and the comparison notebook in 5) |
| [`scplainer_env.yaml`](scplainer_env.yaml) | `r45_scp` | R 4.5, Python 3.14 | `4_VariancePartitioning/01_scplainer.ipynb` |
| [`varpart_env.yaml`](varpart_env.yaml) | `r44_varpart` | R 4.4, Python 3.14 | `4_VariancePartitioning/02_VariancePartitioning.ipynb` |
| [`msqrob_env.yaml`](msqrob_env.yaml) | `r45_msqrob` | R 4.5, Python 3.14 | `5_DifferentialExpression/01_MSQROB.ipynb` |

The full package lists are in the environment files linked above. The R and Python interpreters are pinned to major versions; individual packages are resolved by conda and pip at install time.

### Operating systems

Analysis was performed on Debian GNU/Linux 13.4 ("trixie"), kernel 6.12.85-amd64 (x86_64 architecture). No part of the analysis pipeline was run or tested on macOS or Windows. Raw-file metadata extraction and format conversion (ThermoRawFileParser, tdf2mzml) were executed via pinned Docker container images and are therefore OS-independent of the host system; all custom Python and R analysis code was developed and tested natively on this Linux environment.

### Hardware

No non-standard hardware (no GPU) is required. The R mixed-effects models (MSqRob, variancePartition) and the full precursor-level AnnData objects are memory intensive.

Recommended RAM 32GB minimum

---

## 2. Installation guide

### 1. Clone the repository
```bash
git clone https://github.com/CompOmics/SCI.git
cd SCI
```

### 2. Create the conda environment(s)

For the Python notebooks:
```bash
conda env create -f environment.yml
conda activate sci_environment
pip install -e ./SCI_py
```

The `-e` flag installs the `SCI_py` package in editable mode. For the R notebooks, create the corresponding environment(s) as needed:
```bash
conda env create -f scplainer_env.yaml   # 4_VariancePartitioning/01_scplainer.ipynb
conda env create -f varpart_env.yaml     # 4_VariancePartitioning/02_VariancePartitioning.ipynb
conda env create -f msqrob_env.yaml      # 5_DifferentialExpression/01_MSQROB.ipynb
```

### 3. Register the Jupyter kernels

Each notebook expects a kernel matching its environment (`sci_env`, `R (r45_scp)`, `R (r44_varpart)`, `R (r45_msqrob)`). Step-by-step kernel registration (including the R `IRkernel` setup and a `reticulate` troubleshooting section) is documented in [`Notebooks/4_VariancePartitioning/README.md`](Notebooks/4_VariancePartitioning/README.md) and [`Notebooks/5_DifferentialExpression/README.md`](Notebooks/5_DifferentialExpression/README.md).

### 4. Verify the Python installation
```python
from SCI_py.metrics import correlations, drift
from SCI_py.plot import identifications
import anndata
print("Installation successful.")
```

---

## 3. Data

Pre-processed data objects (AnnData files, metadata, and parsed benchmark results) are distributed via Zenodo: **[https://doi.org/10.5281/zenodo.21181080](https://doi.org/10.5281/zenodo.21181080)**.

Download and extract the archive so that a `Data/` directory sits at the repository root. The notebooks expect, at minimum:

```
Data/
├── metadata/metadata.parquet          # sample-level metadata
├── AnnData/DIANN/...                   # DIA-NN precursor/protein AnnData objects
├── AnnData/Spectronaut/...             # Spectronaut precursor/protein AnnData objects
└── ProteoBench/                        # parsed QC-mixture benchmark objects (1_QCQuant)
```

The exact subfolder layout is set in each notebook's path-configuration cell and matches the archive on Zenodo. To use a different data location, set environment variables before launching JupyterLab:

```bash
export SCI_DATA_DIR=/path/to/your/Data
export SCI_FIGURES_DIR=/path/to/your/Figures
```

All notebooks read these variables at startup and fall back to the default `Data/` and `Figures/` paths (relative to the repository root) if they are not set. Alternatively, you can change the default paths to your specified folders within each notebook.

---

## 4. Demo and expected output

The analysis is delivered as notebooks rather than a command-line tool. Each notebook is self-contained: opening it, selecting the matching kernel, and running all cells reproduces its results with no manual intervention.

**Suggested demo:** with the Zenodo `Data/` extracted and the `sci_environment` kernel selected, run [`Notebooks/2_DatasetStatistics/01_Counts_CVs.ipynb`](Notebooks/2_DatasetStatistics/01_Counts_CVs.ipynb) end to end.

**Expected output:** the notebook prints per-lab identification and CV tables inline and writes the corresponding figures (count boxplots, CV boxplots, missingness heatmaps) to `Figures/<search_engine>/Figure_3/`. Because the notebooks are deterministic given the same input data, re-running a notebook reproduces the same figures. A reference copy of every figure is version-controlled under [`Figures/`](Figures/), so the demo output can be compared directly against the committed files.

---

## 5. Instructions for use

### Reproducing the study
Run the notebooks in the numbered order described in [`Notebooks/README.md`](Notebooks/README.md). In short: build the data objects (folder 0), then run the downstream analyses (folders 1 to 5). Folder-level READMEs describe each notebook, its inputs, and the figures it produces.

### Running on your own data
The notebooks operate on the AnnData objects produced in `0_Create_data_objects`. To analyse a new dataset, generate `metadata.parquet` and the AnnData objects with the folder-0 notebooks (adjusting the input paths and, if needed, the site-specific filename converters in `SCI_py`), then run the downstream notebooks unchanged. Point `SCI_DATA_DIR` and `SCI_FIGURES_DIR` at your own locations to keep inputs and outputs separate from the repository.

---

## Repository structure

```
SCI/
├── environment.yml            Python analysis environment (sci_environment)
├── scplainer_env.yaml         R environment for scplainer (r45_scp)
├── varpart_env.yaml           R environment for variancePartition (r44_varpart)
├── msqrob_env.yaml            R environment for MSqRob (r45_msqrob)
├── LICENSE                    Apache License 2.0
├── Data/                      Analysis data (from Zenodo; not tracked by git)
├── Figures/                   Figure outputs (version-controlled reference copies)
├── SCI_py/                    Python package
│   ├── setup.py
│   └── SCI_py/
│       ├── constants/           Lab maps, condition colours, plotting constants
│       ├── data/                AnnData construction and I/O
│       ├── metadata/            SDRF / cellenONE parsing
│       ├── metrics/             Completeness, correlations, CVs, drift, RT stability
│       ├── plot/                Visualisation modules
│       └── software_comparison/ QC-mixture parsing and benchmark comparison
├── SCI_R/                     R plotting/helper scripts for the R notebooks
└── Notebooks/                 Analysis notebooks (see Notebooks/README.md)
    ├── 0_Create_data_objects/    Data ingestion and AnnData construction
    ├── 1_QCQuant/                Search software benchmarking (QC mixtures)
    ├── 2_DatasetStatistics/      Identification depth, reproducibility, QC, drift
    ├── 3_Correlations/           Sample correlations within and across labs
    ├── 4_VariancePartitioning/   Variance decomposition across covariates
    └── 5_DifferentialExpression/ Differential expression between cell lines
```

---

## Citation

> **TODO :** Citation and DOI.
