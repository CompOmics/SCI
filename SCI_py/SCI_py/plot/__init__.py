"""Figures. Every plot function can save its figure as PNG, SVG and PDF.

Modules with the same name as one in :mod:`SCI_py.metrics` plot its results.

- :mod:`.export` -- :func:`save_figure` and legend helpers used by every plot
- :mod:`.metadata` -- metadata overview (plates, cell types, acquisition timeline)
- :mod:`.plate_heatmap` -- values laid out on the 384-well plate
- :mod:`.identifications` -- identification counts and heatmaps
- :mod:`.completeness`, :mod:`.cvs`, :mod:`.correlations`, :mod:`.drift`,
  :mod:`.rt_stability`, :mod:`.dimensionality_reduction` -- dataset statistics
- :mod:`.software_settings` -- benchmark: search settings compared per tool
- :mod:`.software_performance` -- benchmark: quantitative accuracy per tool
- :mod:`.software_overlap` -- benchmark: overlap between search engines
- :mod:`.software_versions` -- DIA-NN versions and DIA-NN vs Spectronaut
"""

from .export import (
    save_figure,
    save_legend,
    save_legend_from_mapping,
)
from .metadata import (
    plot_cell_types_per_plate,
    plot_qc_flags_per_plate,
    plot_acquisition_timeline,
    plot_plates_per_instrument,
    plot_plate_arrival,
    plot_metadata_overview,
    plot_diameter_per_plate_and_condition,
    plot_injection_order_per_plate,
    plot_plate_position,
    plot_acquisition_date_heatmap,
    plot_cell_qc_overview,
)
from .plate_heatmap import (
    get_well_color,
    plot_plate_heatmap,
    prepare_plate_obs,
)
from .identifications import (
    plot_counts_overview,
    plot_identification_heatmap,
)
from .completeness import plot_completeness
from .cvs import plot_cvs_overview
from .correlations import (
    plot_masking_stats,
    plot_masking_heatmap,
    plot_missingness_distributions,
    plot_coobserved_diagnostic,
    plot_intralab_correlation,
    plot_interlab_correlation_heatmap,
    plot_instrument_correlation,
    plot_sample_pair_correlations,
    plot_intralab_pairwise_stats,
)
from .drift import (
    plot_features_by_day,
    plot_feature_stability,
    plot_drift,
    plot_drift_significance,
    plot_rt_drift_correlations,
    plot_rt_drift_precursor,
    plot_drift_injection,
)
from .rt_stability import plot_rt_stability
from .dimensionality_reduction import (
    plot_pca_comparison,
    plot_embedding_comparison,
)
from .software_settings import (
    plot_scatter_all_tools,
    plot_scatter_per_tool,
    plot_scatter_instrument_split,
    plot_scatter_tool_instrument_grid,
    plot_scatter_recommended,
    plot_scatter_recommended_split,
    plot_scatter_averaged,
)
from .software_performance import (
    plot_species_distributions,
    plot_ma_per_tool,
    plot_ma_multilab,
    plot_dynamic_range_multi,
    compare_tools_error,
)
from .software_overlap import (
    plot_upsetplot_panel,
    plot_full_upsetplot_panel,
    plot_boxplot_panel,
    plot_upsetplot_panel_supplemental,
    plot_boxplot_panel_supplemental,
)
from .software_versions import (
    plot_id_counts_overview,
    plot_id_counts_per_lab,
    plot_protein_counts,
    plot_version_improvement,
    plot_feature_overlap,
    plot_quantitative_correlation,
    plot_tool_counts_per_lab,
    plot_tool_overlap,
    plot_paired_per_lab,
    plot_overlap_decomposition_per_lab,
    plot_overlap_overview,
    plot_binned_overlap_boxplot,
    plot_binned_overlap_profile,
    plot_completeness_multi_dataset,
    plot_cv_multi_dataset,
)

__all__ = [
    "save_figure",
    "save_legend",
    "save_legend_from_mapping",
    "plot_cell_types_per_plate",
    "plot_qc_flags_per_plate",
    "plot_acquisition_timeline",
    "plot_plates_per_instrument",
    "plot_plate_arrival",
    "plot_metadata_overview",
    "plot_diameter_per_plate_and_condition",
    "plot_injection_order_per_plate",
    "plot_plate_position",
    "plot_acquisition_date_heatmap",
    "plot_cell_qc_overview",
    "get_well_color",
    "plot_plate_heatmap",
    "prepare_plate_obs",
    "plot_counts_overview",
    "plot_identification_heatmap",
    "plot_completeness",
    "plot_cvs_overview",
    "plot_masking_stats",
    "plot_masking_heatmap",
    "plot_missingness_distributions",
    "plot_coobserved_diagnostic",
    "plot_intralab_correlation",
    "plot_interlab_correlation_heatmap",
    "plot_instrument_correlation",
    "plot_sample_pair_correlations",
    "plot_intralab_pairwise_stats",
    "plot_features_by_day",
    "plot_feature_stability",
    "plot_drift",
    "plot_drift_significance",
    "plot_rt_drift_correlations",
    "plot_rt_drift_precursor",
    "plot_drift_injection",
    "plot_rt_stability",
    "plot_pca_comparison",
    "plot_embedding_comparison",
    "plot_scatter_all_tools",
    "plot_scatter_per_tool",
    "plot_scatter_instrument_split",
    "plot_scatter_tool_instrument_grid",
    "plot_scatter_recommended",
    "plot_scatter_recommended_split",
    "plot_scatter_averaged",
    "plot_species_distributions",
    "plot_ma_per_tool",
    "plot_ma_multilab",
    "plot_dynamic_range_multi",
    "compare_tools_error",
    "plot_upsetplot_panel",
    "plot_full_upsetplot_panel",
    "plot_boxplot_panel",
    "plot_upsetplot_panel_supplemental",
    "plot_boxplot_panel_supplemental",
    "plot_id_counts_overview",
    "plot_id_counts_per_lab",
    "plot_protein_counts",
    "plot_version_improvement",
    "plot_feature_overlap",
    "plot_quantitative_correlation",
    "plot_tool_counts_per_lab",
    "plot_tool_overlap",
    "plot_paired_per_lab",
    "plot_overlap_decomposition_per_lab",
    "plot_overlap_overview",
    "plot_binned_overlap_boxplot",
    "plot_binned_overlap_profile",
    "plot_completeness_multi_dataset",
    "plot_cv_multi_dataset",
]
