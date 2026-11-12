"""SCI_py: analysis code for the SCI single-cell proteomics multi-lab study.

Subpackages, roughly in the order the notebooks use them:

- :mod:`SCI_py.constants` -- fixed study facts: labs, samples, colours
- :mod:`SCI_py.metadata` -- build the per-sample metadata table
- :mod:`SCI_py.data` -- load search-engine reports into AnnData
- :mod:`SCI_py.benchmark` -- search-engine comparison on the HYE QC/Quant samples
- :mod:`SCI_py.metrics` -- computations (no plotting)
- :mod:`SCI_py.plot` -- figures; same module names as :mod:`SCI_py.metrics`

Import from the subpackage, e.g. ``from SCI_py.plot import plot_drift``.
"""
