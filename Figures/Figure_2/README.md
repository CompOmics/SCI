# Figure 2 — Multi-engine Benchmarking (ProteoBench)

These figures compare DIA search engine performance (DIA-NN, Spectronaut, AlphaDIA, PEAKS,
Chimerys/MSAID, FragPipe) across nine labs (six Astral, three timsTOF) using a three-species
benchmark mix (human, yeast, *E. coli*). Metrics: precursor identification count and absolute
quantification error (mean |log2 FC error|).

**All files are saved in png + svg + pdf unless noted otherwise.**

---

## Subfolders

| Subfolder | Topic | Generating notebook |
|-----------|-------|---------------------|
| [overview/](#overview) | Benchmark scatter plots | `1_QCQuant/02_Settings_comparison.ipynb` |
| [dynamic_range/](#dynamic-range) | MA-plots and pairwise error | `1_QCQuant/03_Performance_comparison.ipynb` |
| [overlap_error/](#overlap-and-quantification-error) | UpSet + boxplot panels | `1_QCQuant/04_Overlap_between_search_engines.ipynb` |
| [precursor_characteristics/](#precursor-characteristics) | Peptide feature distributions | `1_QCQuant/05_PrecursorCharacteristics.ipynb` |

---

## Overview

**Notebook:** [`1_QCQuant/02_Settings_comparison.ipynb`](../../Notebooks/1_QCQuant/02_Settings_comparison.ipynb)

Scatter plots of precursor count vs mean absolute quantification error, one point per
lab × search engine × settings combination.

| File | Description |
|------|-------------|
| `scatter_all_tools` | All tools and all settings, colour by search engine |
| `scatter_per_tool` | One subplot per search engine, colour by lab |
| `scatter_instrument_split` | Astral (left) vs timsTOF (right), all settings, colour by search engine |
| `scatter_tool_instrument_grid` | Grid: rows = search engines, columns = Astral / timsTOF |
| `scatter_recommended` | Recommended (no-MBR) and MBR settings side by side |
| `scatter_recommended_split_no_mbr` | Recommended settings, split by instrument type |
| `scatter_recommended_split_mbr` | MBR settings, split by instrument type |
| `scatter_averaged` | Instrument-averaged mean performance per search engine |

---

## Dynamic range

**Notebook:** [`1_QCQuant/03_Performance_comparison.ipynb`](../../Notebooks/1_QCQuant/03_Performance_comparison.ipynb)

Compares measured log2 fold changes against the expected species ratios. Labs: Astral_1–6,
timsTOF_1–3.

| File | Description |
|------|-------------|
| `species_distributions` | KDE of measured log2 FC per lab and species (HUMAN / YEAST / ECOLI) |
| `MAplots_{lab}` (×9) | MA-plots for one lab: mean log intensity vs measured log2 FC, one panel per search engine |
| `pairwise_error_DIANN_vs_Spectronaut_{lab}` (×9) | Top: mean absolute error vs abundance rank per species; bottom: delta error (which tool wins) |

---

## Overlap and quantification error

**Notebook:** [`1_QCQuant/04_Overlap_between_search_engines.ipynb`](../../Notebooks/1_QCQuant/04_Overlap_between_search_engines.ipynb)

Panels arranged as 3×3 grids (one cell per lab).

| File | Description |
|------|-------------|
| `upsetplot_panel_3x3` | Precursor overlap across six search engines (unique / all-but-one / all-tools subsets) |
| `boxplot_panel_3x3` | `|log2 FC error|` per search engine, split by unique-to-tool vs all-tools precursors |

---

## Precursor characteristics

**Notebook:** [`1_QCQuant/05_PrecursorCharacteristics.ipynb`](../../Notebooks/1_QCQuant/05_PrecursorCharacteristics.ipynb)

Saved for three instrument groups: `All/`, `Astral/`, `timsTOF/`.

| File | Description |
|------|-------------|
| `pep_length` | Peptide length distribution (two panels: all precursors vs excluding all-tools intersection) |
| `charge_state` | Charge state distribution |
| `missed_cleavages` | Missed cleavage distribution |
| `modifications` | Modification class distribution |

---

[← Back to Figures overview](../README.md)
