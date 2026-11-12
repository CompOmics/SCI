Where to find things
====================

The package has six subpackages. Always import from the subpackage, for
example ``from SCI_py.plot import plot_drift``, not from the module file.

.. code-block:: text

   SCI_py/
   ├── constants/   fixed study facts
   │   ├── labs.py        lab ids, anonymised names (Astral_1 … timsTOF_3), conditions
   │   ├── samples.py     plate layout, QC/benchmark run names, SDRF value maps
   │   └── plotting.py    colours and display order
   ├── metadata/    build the per-sample metadata table
   │   ├── read.py        read SDRF, cellenONE and docx annotations
   │   ├── columns.py     clean, impute and type the merged table
   │   └── filenames.py   raw-file name → lab and sample name
   ├── data/        search-engine reports → AnnData
   │   ├── reports.py     load DIA-NN / Spectronaut reports
   │   └── anndata_utils.py
   ├── benchmark/   search-engine comparison on the HYE QC/Quant samples
   │   ├── readers.py, entries.py, variants.py, cache.py, fasta.py, filters.py
   │   ├── parameters/    parse each engine's settings files
   │   └── proteobench/   ProteoBench conversion (copied from ProteoBench)
   ├── metrics/     computations, no plotting
   └── plot/        figures; same module names as metrics/

Rules of thumb
--------------

- **metrics computes, plot draws.** A module in :mod:`SCI_py.metrics` returns
  numbers or tables; the module with the same name in :mod:`SCI_py.plot`
  draws them. ``metrics.drift`` ↔ ``plot.drift``, ``metrics.cvs`` ↔
  ``plot.cvs``, and so on.
- **software_*** modules belong to the search-engine comparison:
  ``plot.software_settings``, ``plot.software_performance``,
  ``plot.software_overlap`` and ``plot.software_versions``.
- **Saving figures.** Plot functions take ``save_dir`` and ``filename``
  and write PNG, SVG and PDF through :func:`~SCI_py.plot.export.save_figure`.
- **Lab names.** Internal lab ids (``IMP_Astral_1``) are mapped to the
  anonymised names used in figures (``Astral_5``) by
  ``LAB_NAME_MAPPING`` in :mod:`SCI_py.constants.labs`.

I want to …
-----------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Task
     - Use
   * - build the metadata table
     - :func:`~SCI_py.metadata.read.read_SDRF`, :func:`~SCI_py.metadata.read.read_cellenONE`,
       :func:`~SCI_py.metadata.read.add_docx_annotations`, :func:`~SCI_py.metadata.columns.parse_metadata_columns`
   * - load a DIA-NN or Spectronaut report
     - :func:`~SCI_py.data.reports.load_diann`, :func:`~SCI_py.data.reports.load_spectronaut`
   * - turn a report into AnnData
     - :func:`~SCI_py.data.anndata_utils.df_to_anndata`, :func:`~SCI_py.data.anndata_utils.align_obs_across_anndata`
   * - compare search engines on the QC/Quant samples
     - :mod:`SCI_py.benchmark` (pipeline), then the ``plot.software_*`` modules
   * - compute completeness / CVs / correlations / drift
     - :mod:`SCI_py.metrics`
   * - plot them
     - :mod:`SCI_py.plot`
   * - get lab colours and order
     - ``LAB_COLORS``, ``LAB_ORDER`` in :mod:`SCI_py.constants.plotting`
