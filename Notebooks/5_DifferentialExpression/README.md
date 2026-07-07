# 5_DifferentialExpression

Differential expression analysis between the two cell lines (K562 versus HEK293T) and comparison of the results across search engines. Figures produced here correspond to Figure 6 of the manuscript. The notebooks use the R kernel (`msqrob2`) and Python for the cross-engine comparison.

## Notebooks

### [01_MSQROB.ipynb](01_MSQROB.ipynb)  (R notebook)
Fits `msqrob2` ridge-penalised robust mixed-effects models for differential protein abundance between K562 and HEK293T. Precursor intensities are aggregated to proteins, then four integrated models are fitted that add one layer of technical covariates at a time, and each plate is additionally fitted independently. Produces volcano-plot grids, UpSet plots of significant-protein overlap across models, logFC and standard-error comparisons, batch-corrected PCA sanity checks, and integrated-versus-per-plate overlap analysis.

### [02_DIANN_Spectronaut_overlap.ipynb](02_DIANN_Spectronaut_overlap.ipynb)
Compares the MSqRob DE results obtained from DIA-NN and Spectronaut on the same samples. Covers protein detection and overlap across the four result sets, logFC concordance on shared proteins, significance overlap, per-plate logFC correlations, direction-of-effect agreement, and identification of opposite-direction discordant proteins (significant in both tools but with opposite sign). Figures feed the software-comparison section of the supplement.

> `01_MSQROB.ipynb` is an R notebook and runs in the `r45_msqrob` environment (setup below). `02_DIANN_Spectronaut_overlap.ipynb` is a Python notebook and runs in the main `sci_environment` kernel described in the repository root [`README.md`](../../README.md).

---

# MSqRob Notebook: Setup Instructions

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
- [Jupyter](https://jupyter.org/install) (or JupyterLab)
- Git (to clone this repository)

---

## 1. Create the conda environment

```bash
conda env create -f msqrob_env.yaml
```

This installs all required Python and R packages. Depending on your machine this may take 10–20 minutes.

To verify the environment was created:

```bash
conda env list
```

You should see `r45_msqrob` in the list.

---

## 2. Activate the environment

```bash
conda activate r45_msqrob
```

---

## 3. Register the environment as a Jupyter kernel

With the environment activated, register it so Jupyter can find it:

```bash
python -m ipykernel install --user --name r45_msqrob --display-name "Python (r45_msqrob)"
```

---

## 4a. Register the R kernel (IRkernel)

Still with the environment activated, launch R and install the kernel:

```bash
R -e "install.packages('IRkernel', repos='https://cloud.r-project.org')"
R -e "IRkernel::installspec(name='r45_msqrob_r', displayname='R (r45_msqrob)')"
```

---

## 5. Verify the kernels are registered

```bash
jupyter kernelspec list
```

You should see both `r45_msqrob` and `r45_msqrob_r` in the output.

---

## 6. Launch Jupyter

```bash
jupyter lab
```

Or if you prefer classic Jupyter:

```bash
jupyter notebook
```

---

## 7. Select the correct kernel

When opening a notebook:

- For **R notebooks** (`01_MSQROB.ipynb`): select `R (r45_msqrob)`

To change the kernel in an already open notebook: **Kernel → Change Kernel** from the menu.

---

## Updating the environment

If the `msqrob_env.yaml` file is updated:

```bash
conda env update -f msqrob_env.yaml --prune
```

The `--prune` flag removes packages that are no longer in the file.

---

## Removing the environment

```bash
conda deactivate
conda env remove -n r45_msqrob
```

To also remove the Jupyter kernels:

```bash
jupyter kernelspec remove r45_msqrob
jupyter kernelspec remove r45_msqrob_r
```

---

## Troubleshooting

**`reticulate` cannot find the conda environment in R**

Make sure the Python path is set before loading `reticulate`. At the top of your R notebook:

```r
Sys.setenv(RETICULATE_PYTHON = path.expand("~/miniconda3/envs/r45_msqrob/bin/python"))
library(reticulate)
use_condaenv("r45_msqrob", required = TRUE)
```

Adjust the path to match your Miniconda installation directory.

**Kernel not appearing in Jupyter**

Re-run the kernel registration step (step 3 or 4) with the environment activated, then restart Jupyter.

**Package missing after environment creation**

Some packages may need to be installed manually if they are not available via conda:

```bash
conda activate r45_msqrob
pip install <package>
# or
R -e "install.packages('<package>')"
```
