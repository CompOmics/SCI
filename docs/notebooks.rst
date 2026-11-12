Functions used per notebook
============================

Which ``SCI_py`` names each notebook imports, linked to their documentation.
Notebooks in folders 4 and 5 are R notebooks and do not use the package.

0_Create_data_objects
---------------------

``01_metadata.ipynb``

- :mod:`SCI_py.metadata`: :func:`~SCI_py.metadata.read.read_SDRF`, :func:`~SCI_py.metadata.read.read_cellenONE`, :func:`~SCI_py.metadata.read.add_docx_annotations`, :func:`~SCI_py.metadata.columns.parse_metadata_columns`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.metadata.plot_metadata_overview`, :func:`~SCI_py.plot.metadata.plot_cell_qc_overview`

``02_Results_to_AnnData_DIANN.ipynb``

- :mod:`SCI_py.data`: :func:`~SCI_py.data.anndata_utils.df_to_anndata`, :func:`~SCI_py.data.anndata_utils.align_obs_across_anndata`, :func:`~SCI_py.data.anndata_utils.replace_zeros_with_nan`, :func:`~SCI_py.data.reports.load_spectronaut`, :func:`~SCI_py.data.reports.load_diann`

``02_Results_to_AnnData_Spectronaut.ipynb``

- :mod:`SCI_py.data`: :func:`~SCI_py.data.anndata_utils.df_to_anndata`, :func:`~SCI_py.data.anndata_utils.align_obs_across_anndata`, :func:`~SCI_py.data.reports.load_spectronaut`

1_QCQuant
---------

``01_parse_results.ipynb``

- :mod:`SCI_py.constants`: ``LAB_NAME_MAPPING``
- :mod:`SCI_py.benchmark`: :func:`~SCI_py.benchmark.proteobench.convert.parse_tables_to_proteobench`, :func:`~SCI_py.benchmark.fasta.parse_fasta_headers`, :func:`~SCI_py.benchmark.readers.read_result_folder`, :func:`~SCI_py.benchmark.readers.read_diann_firstpass`, :func:`~SCI_py.benchmark.variants.apply_tic_norm`, :func:`~SCI_py.benchmark.entries.enumerate_entries`, :class:`~SCI_py.benchmark.cache.EntryCache`

``02_Settings_comparison.ipynb``

- :mod:`SCI_py.constants`: ``SOFTWARE_COLORS``, ``LAB_COLORS``, ``LAB_ORDER``
- :mod:`SCI_py.benchmark`: :func:`~SCI_py.benchmark.filters.filter_recommended`, :func:`~SCI_py.benchmark.filters.apply_version_filter`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.software_settings.plot_scatter_all_tools`, :func:`~SCI_py.plot.software_settings.plot_scatter_per_tool`, :func:`~SCI_py.plot.software_settings.plot_scatter_instrument_split`, :func:`~SCI_py.plot.software_settings.plot_scatter_tool_instrument_grid`, :func:`~SCI_py.plot.software_settings.plot_scatter_recommended`, :func:`~SCI_py.plot.software_settings.plot_scatter_recommended_split`, :func:`~SCI_py.plot.software_settings.plot_scatter_averaged`

``03_Performance_comparison.ipynb``

- :mod:`SCI_py.constants`: ``SOFTWARE_COLORS``
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.software_performance.plot_species_distributions`, :func:`~SCI_py.plot.software_performance.plot_ma_per_tool`, :func:`~SCI_py.plot.software_performance.plot_ma_multilab`, :func:`~SCI_py.plot.software_performance.plot_dynamic_range_multi`, :func:`~SCI_py.plot.software_performance.compare_tools_error`

``04_Overlap_between_search_engines.ipynb``

- :mod:`SCI_py.constants`: ``LAB_COLORS``, ``LAB_ORDER``, ``SOFTWARE_COLORS``
- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.software_overlap.intermediate_to_upsetplotdf`, :func:`~SCI_py.metrics.software_overlap.get_boxplot_data`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.software_overlap.plot_upsetplot_panel`, :func:`~SCI_py.plot.software_overlap.plot_boxplot_panel`, :func:`~SCI_py.plot.software_overlap.plot_upsetplot_panel_supplemental`, :func:`~SCI_py.plot.software_overlap.plot_boxplot_panel_supplemental`, :func:`~SCI_py.plot.software_overlap.plot_full_upsetplot_panel`

``05_PrecursorCharacteristics.ipynb``

- :mod:`SCI_py.constants`: ``SOFTWARE_COLORS``, ``LAB_ORDER``, ``SAMPLE_TYPE_A``, ``SAMPLE_TYPE_B``
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.export.save_figure`

2_DatasetStatistics
-------------------

``00_SearchEngineComparison.ipynb``

- :mod:`SCI_py.constants`: ``LAB_COLORS``, ``LAB_ORDER``, ``CONDITION_COLORS``
- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.completeness.compute_completeness`, :func:`~SCI_py.metrics.cvs.compute_cv_matrix`, :func:`~SCI_py.metrics.software_versions.shared_obs`, :func:`~SCI_py.metrics.software_versions.compute_overlap_decomposition`, :func:`~SCI_py.metrics.software_versions.compute_binned_overlap`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.software_versions.plot_id_counts_overview`, :func:`~SCI_py.plot.software_versions.plot_id_counts_per_lab`, :func:`~SCI_py.plot.software_versions.plot_protein_counts`, :func:`~SCI_py.plot.software_versions.plot_version_improvement`, :func:`~SCI_py.plot.software_versions.plot_feature_overlap`, :func:`~SCI_py.plot.software_versions.plot_quantitative_correlation`, :func:`~SCI_py.plot.software_versions.plot_tool_counts_per_lab`, :func:`~SCI_py.plot.software_versions.plot_tool_overlap`, :func:`~SCI_py.plot.software_versions.plot_paired_per_lab`, :func:`~SCI_py.plot.software_versions.plot_overlap_decomposition_per_lab`, :func:`~SCI_py.plot.software_versions.plot_overlap_overview`, :func:`~SCI_py.plot.software_versions.plot_binned_overlap_boxplot`, :func:`~SCI_py.plot.software_versions.plot_binned_overlap_profile`, :func:`~SCI_py.plot.software_versions.plot_completeness_multi_dataset`, :func:`~SCI_py.plot.software_versions.plot_cv_multi_dataset`

``01_Counts_CVs.ipynb``

- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.cvs.compute_cv_matrix`, :func:`~SCI_py.metrics.cvs.filter_common_proteome`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.identifications.plot_counts_overview`, :func:`~SCI_py.plot.cvs.plot_cvs_overview`, :func:`~SCI_py.plot.identifications.plot_identification_heatmap`
- :mod:`SCI_py.constants`: ``CONDITION_ORDER``
- :mod:`SCI_py.data`: :func:`~SCI_py.data.anndata_utils.inspect_anndata_size`

``02_QC_Metrics.ipynb``

- :mod:`SCI_py.constants`: ``LABS``, ``CONDITION_COLORS``, ``LAB_ORDER``
- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.completeness.compute_completeness`, :func:`~SCI_py.metrics.completeness.compute_combined_curve`, :func:`~SCI_py.metrics.completeness.store_completeness`, :func:`~SCI_py.metrics.rt_stability.compute_rt_proportions`, :func:`~SCI_py.metrics.rt_stability.compute_rt_norm_layer`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.completeness.plot_completeness`, :func:`~SCI_py.plot.plate_heatmap.plot_plate_heatmap`, :func:`~SCI_py.plot.rt_stability.plot_rt_stability`, :func:`~SCI_py.plot.plate_heatmap.get_well_color`, :func:`~SCI_py.plot.plate_heatmap.prepare_plate_obs`, :func:`~SCI_py.plot.export.save_figure`

``03_Drift_over_time.ipynb``

- :mod:`SCI_py.constants`: ``LABS``, ``CONDITION_COLORS``
- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.drift.compute_drift`, :func:`~SCI_py.metrics.drift.fit_drift_stats`, :func:`~SCI_py.metrics.drift.compute_rt_drift_correlations`, :func:`~SCI_py.metrics.drift.compute_time_norm`, :func:`~SCI_py.metrics.rt_stability.compute_rt_proportions`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.drift.plot_drift_significance`, :func:`~SCI_py.plot.drift.plot_drift`, :func:`~SCI_py.plot.drift.plot_features_by_day`, :func:`~SCI_py.plot.drift.plot_feature_stability`, :func:`~SCI_py.plot.drift.plot_rt_drift_precursor`, :func:`~SCI_py.plot.drift.plot_rt_drift_correlations`, :func:`~SCI_py.plot.drift.plot_drift_injection`, :func:`~SCI_py.plot.rt_stability.plot_rt_stability`, :func:`~SCI_py.plot.export.save_figure`

3_Correlations
--------------

``01_sample_correlations.ipynb``

- :mod:`SCI_py.metrics`: :func:`~SCI_py.metrics.correlations.create_reference_set`, :func:`~SCI_py.metrics.correlations.compute_pairwise_correlations`, :func:`~SCI_py.metrics.correlations.compute_pairwise_stats`
- :mod:`SCI_py.plot`: :func:`~SCI_py.plot.correlations.plot_masking_stats`, :func:`~SCI_py.plot.correlations.plot_masking_heatmap`, :func:`~SCI_py.plot.correlations.plot_missingness_distributions`, :func:`~SCI_py.plot.correlations.plot_coobserved_diagnostic`, :func:`~SCI_py.plot.correlations.plot_intralab_correlation`, :func:`~SCI_py.plot.correlations.plot_interlab_correlation_heatmap`, :func:`~SCI_py.plot.correlations.plot_instrument_correlation`, :func:`~SCI_py.plot.correlations.plot_sample_pair_correlations`, :func:`~SCI_py.plot.correlations.plot_intralab_pairwise_stats`
