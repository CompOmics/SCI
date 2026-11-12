SCI_py
======

``SCI_py`` holds the analysis code for the HUPO Single-Cell Initiative (SCI)
inter-laboratory benchmarking study: nine 384-well plates from seven labs,
measured on Orbitrap Astral and timsTOF instruments, searched with DIA-NN and
Spectronaut.

The notebooks in ``Notebooks/`` do the analyses; this package holds the code
they share.

.. code-block:: bash

   conda activate sci_environment
   pip install -e ./SCI_py

.. code-block:: python

   from SCI_py.constants import LAB_COLORS, LAB_ORDER
   from SCI_py.metrics import compute_cv_matrix
   from SCI_py.plot import plot_cvs_overview

.. toctree::
   :maxdepth: 2

   layout
   notebooks
   api
