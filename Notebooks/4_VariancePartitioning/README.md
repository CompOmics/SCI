# 4_VariancePartitioning

Variance decomposition of the single-cell quantification, quantifying how much variability is explained by technical factors (normalisation, plate effects, injection position) versus biological factors (cell type, cell size). The two notebooks apply complementary frameworks at the precursor and protein level and each proceed through the same four successive correction stages. Figures produced here correspond to Figure 6 of the manuscript. Both accept the `SCI_DATA_DIR` and `SCI_FIGURES_DIR` environment variables to switch between search engines.

Environment setup instructions for the required conda environment and Jupyter kernels follow below the notebook overview. Both are R-notebooks run in jupyter notebooks.

## Notebooks

### [01_scplainer.ipynb](01_scplainer.ipynb)
Fits scplainer (`scp`) linear models to precursor-level log-intensities across four correction stages. Decomposes per-precursor variance by covariate (average bar chart and per-feature violin distributions), runs GMF probabilistic matrix factorisation for dimensionality reduction with rank selection, evaluates batch mixing and biological separation via silhouette scores, and runs a moderated t-test for differential abundance per stage. Exports the batch-corrected object to AnnData for downstream Python analyses.

### [02_VariancePartitioning.ipynb](02_VariancePartitioning.ipynb) (informative)
Fits `variancePartition` mixed-effects models on protein-level log2 intensities. Applies protein filtering (validated with a retention heatmap) and row-mean imputation, then fits three models corresponding to the correction stages, using `limma::lmFit` residuals to regress out technical effects. Produces per-covariate variance-fraction violin plots and PCA plots per stage, and exports the variance-fraction matrices as CSV and RDS.

---

# SCP Analysis Notebook: Setup Instructions

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
- [Jupyter](https://jupyter.org/install) (or JupyterLab)
- Git (to clone this repository)

---

## 1. Create the conda environment

```bash
conda env create -f scplainer_env.yaml
```

This installs all required Python and R packages. Depending on your machine this may take 10–20 minutes.

To verify the environment was created:

```bash
conda env list
```

You should see `r45_scp` in the list.

---

## 2. Activate the environment

```bash
conda activate r45_scp
```

---

## 3. Register the environment as a Jupyter kernel

With the environment activated, register it so Jupyter can find it:

```bash
python -m ipykernel install --user --name r45_scp --display-name "Python (r45_scp)"
```

---

## 4a. Register the R kernel (IRkernel)

Still with the environment activated, launch R and install the kernel:

```bash
R -e "install.packages('IRkernel', repos='https://cloud.r-project.org')"
R -e "IRkernel::installspec(name='r45_scp_r', displayname='R (r45_scp)')"
```

---

## 5. Verify the kernels are registered

```bash
jupyter kernelspec list
```

You should see both `r45_scp` and `r45_scp_r` in the output.

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

- For **R notebooks**: select `R (r45_scp)`

To change the kernel in an already open notebook: **Kernel → Change Kernel** from the menu.

---

## Updating the environment

If the `scplainer_env.yaml` file is updated:

```bash
conda env update -f scplainer_env.yaml --prune
```

The `--prune` flag removes packages that are no longer in the file.

---

## Removing the environment

```bash
conda deactivate
conda env remove -n r45_scp
```

To also remove the Jupyter kernels:

```bash
jupyter kernelspec remove r45_scp
jupyter kernelspec remove r45_scp_r
```

---

## Troubleshooting

**`reticulate` cannot find the conda environment in R**

Make sure the Python path is set before loading `reticulate`. At the top of your R notebook:

```r
Sys.setenv(RETICULATE_PYTHON = path.expand("~/miniconda3/envs/r45_scp/bin/python"))
library(reticulate)
use_condaenv("r45_scp", required = TRUE)
```

Adjust the path to match your Miniconda installation directory.

**Kernel not appearing in Jupyter**

Re-run the kernel registration step (step 3 or 4) with the environment activated, then restart Jupyter.

**Package missing after environment creation**

Some packages may need to be installed manually if they are not available via conda:

```bash
conda activate r45_scp
pip install <package>
# or
R -e "install.packages('<package>')"
```