# Data Object Construction

This directory contains the notebooks that transform raw instrument output and metadata into the analysis-ready data objects used throughout the SCI project. **These notebooks must be run in order before any downstream analysis notebook.**

---

## Overview

The SCI dataset integrates data from **9 participating laboratories** across two instrument platforms (Orbitrap Astral, timsTOF). Raw acquisitions were processed centrally through Spectronaut and the resulting quantification tables, together with sample-level metadata from two independent sources (SDRF files and cellenONE dispensing reports), are assembled here into a unified set of data objects.

Two object types are produced and maintained in parallel throughout the project:

| Object type | Format | Primary use |
|---|---|---|
| **AnnData** | `.h5ad` (Python / `anndata`) | Exploratory analysis, visualisation, differential expression |
| **QFeatures** | `.rds` (R / `QFeatures`) | Scplainer modelling, variance partitioning |

---

## Data Flow

```
Raw Spectronaut TSV                    SDRF.xlsx (per lab)
        │                              cellenONE reports (per plate)
        │                                         │
        ▼                                         ▼
┌───────────────────┐             ┌───────────────────────────┐
│  02_Results_to_   │             │     01_metadata.ipynb     │
│  AnnData.ipynb    │◄────────────│                           │
│                   │  metadata   │  Merges SDRF + cellenONE, │
│  Pivots long-     │  .parquet   │  applies 8-step parsing   │
│  format results   │             │  and imputation pipeline  │
│  to AnnData       │             └───────────────────────────┘
└────────┬──────────┘                          │
         │                                     │ metadata.parquet
         │ adata_precursor.h5ad                ▼
         │ adata_protein.h5ad        Data/metadata/
         │ adata_precursor_qc.h5ad
         │ adata_protein_qc.h5ad
         ▼
    Data/AnnData/
         │
         ├──────────────────────────────────────┐
         │                                      │
         ▼                                      ▼
┌─────────────────────┐             ┌───────────────────────┐
│ 03_AnnData_to_      │             │ Analysis notebooks    │
│ QFeatures.ipynb     │             │ (1_QCQuant,           │
│                     │             │  2_DatasetStatistics, │
│ Converts AnnData    │             │  3_Correlations,      │
│ objects to R-side   │             │  4_VariancePartition, │
│ QFeatures for       │             │  5_DiffExpression)    │
│ Scplainer and       │             └──────────┬────────────┘
│ VariancePartition   │                        │ results written
└─────────────────────┘                        │ back to AnnData
         │                                     │
         │ QFeatures .rds                      ▼
         ▼                          ┌───────────────────────┐
    Data/QFeatures/                 │ 04_Assembly_AnnData   │
                                    │ .ipynb                │
                                    │                       │
                                    │ Consolidates analysis │
                                    │ results from all      │
                                    │ notebooks into a      │
                                    │ single, annotated     │
                                    │ AnnData object        │
                                    └───────────────────────┘
                                             │
                                             ▼
                                    Data/AnnData/
                                    adata_assembled.h5ad
```

---

## Notebooks

### `01_metadata.ipynb` — Metadata Construction

**Run first. Required by all other notebooks.**

Consolidates sample-level metadata from two independent sources into a single parquet file.

**Inputs**
- `metadata/` — one `SDRF.xlsx` per lab directory
- `cellenONE/cellenONE_Reports/` — one `Reordered_*_isolated.xls` per plate × cell-type

**Processing steps**
1. **SDRF ingestion** — per-lab filename normalisation to canonical sample names, plate number mapping, well coordinate conversion (letter → integer, XY → `R{row}-C{col}`), reinjection detection (Wyss: supplementary CSV; Broad: `"redo"` in filename), and a BYU source-name label correction (dispensing error swapping QC and blank well identities).
2. **cellenONE ingestion** — merges isolation reports from all plates and cell types, attaching morphological measurements (diameter, elongation, circularity) and per-cell image paths.
3. **Merge** — SDRF and cellenONE joined on plate number + XY position. Manual well-level QC annotations (`Dying`, `Dividing`, `Doublet`, `CellInBlank`) from plate preparation documentation are applied as boolean columns.
4. **Parsing pipeline** — sample name lowercasing; temporal column decomposition (plate arrival, gradient time, thawing timepoint → `_date`, `_time`, `_datetime`); derivation of `instrument_type`, `sample_id`, and anonymised `plate_idx`; source name normalisation (`Cell A/B` → `K562/HEK293T`); normalised injection order per plate; scalar and lookup-based missing value imputation; type enforcement.

**Outputs**
- `Data/metadata/metadata.parquet`

---

### `02_Results_to_AnnData.ipynb` — Search Results → AnnData

**Depends on:** `01_metadata.ipynb`

Converts raw Spectronaut quantification tables into AnnData objects at precursor and protein level.

**Inputs**
- Spectronaut BGS Factory Report (Normal) TSV — main SCI run (all labs merged)
- Spectronaut Factory Report TSVs — site-specific quantitative QC mixture runs
- `Data/metadata/metadata.parquet`

**Processing steps**
1. **Spectronaut parsing** — ProteoBench loader with modification notation standardisation to UNIMOD accessions (`[Oxidation (M)]` → `[UNIMOD:35]`, etc.) and construction of a canonical precursor ID (`sequence + charge`).
2. **Metadata merge** — joins search results to sample metadata on `plate_id` + `sample_name`.
3. **Column partitioning** — columns are grouped into feature annotations (`var`), sample annotations (`obs`), and quantitative values (`X`) before pivoting.
4. **AnnData construction** — `df_to_anndata` pivots the long-format table into a sparse `obs × var` intensity matrix. `align_obs_across_anndata` enforces a shared observation index between the precursor and protein objects.
5. **Log transformation** — `log2(X)` stored as the `X_log2` layer alongside the raw intensities in `X`.

**Outputs** (`Data/AnnData/`, `Data/Parquet/`)

| File | Description |
|---|---|
| `adata_precursor.h5ad` | Precursor-level, human samples |
| `adata_protein.h5ad` | Protein-level, human samples |
| `adata_precursor_qc.h5ad` | Precursor-level, QC mixtures |
| `adata_protein_qc.h5ad` | Protein-level, QC mixtures |
| `search_results.parquet` | Filtered long-format results, human samples |
| `search_results_qc.parquet` | Filtered long-format results, QC mixtures |

---

### `03_AnnData_to_QFeatures.ipynb` — AnnData → QFeatures *(to be written)*

**Depends on:** `02_Results_to_AnnData.ipynb`

Transfers the AnnData objects to the R ecosystem for statistical modelling. The QFeatures format (Bioconductor) is the expected input for **Scplainer** (single-cell proteomics linear modelling) and **VariancePartition** (variance component decomposition). This notebook handles the cross-language handoff, ensuring that observation metadata, feature annotations, and quantification matrices are faithfully preserved in the R-side representation.

**Inputs**
- `Data/AnnData/adata_precursor.h5ad`
- `Data/AnnData/adata_protein.h5ad`

**Outputs** (`Data/QFeatures/`)
- `qfeatures_precursor.rds`
- `qfeatures_protein.rds`

---

### `04_Assembly_AnnData.ipynb` — Consolidation *(to be written)*

**Depends on:** all analysis notebooks (`1_QCQuant` through `5_DifferentialExpression`)

Assembles a single, fully annotated AnnData object by collecting results written back to disk by the individual analysis notebooks. Most analysis notebooks update the AnnData objects in place (adding embeddings, cluster labels, statistical results to `obs`, `var`, or `uns`); this notebook performs a final consolidation into a canonical object suitable for public release and figure reproduction.

**Inputs**
- AnnData objects as updated by the analysis notebooks

**Outputs**
- `Data/AnnData/adata_assembled.h5ad` — complete, annotated dataset

---

## Execution Order

```
01_metadata.ipynb
       ↓
02_Results_to_AnnData.ipynb
       ↓
03_AnnData_to_QFeatures.ipynb   ←── required before running 4_VariancePartitioning
       ↓
   [analysis notebooks: 1_QCQuant → 5_DifferentialExpression]
       ↓
04_Assembly_AnnData.ipynb
```

> Notebooks `03` and `04` only need to be re-run when the upstream AnnData objects change (e.g. after reprocessing the Spectronaut results or adding new plates). Analysis notebooks can generally be re-run independently once the AnnData objects exist.

---

## Environment

Python notebooks use the `sci_env` conda environment. R notebooks (`03`) require the `QFeatures`, `scpdata`, and `VariancePartition` Bioconductor packages. See the top-level `README.md` for installation instructions.
