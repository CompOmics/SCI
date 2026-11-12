"""Computations on the AnnData objects and benchmark tables. No plotting here.

Each module has a plotting counterpart with the same name in :mod:`SCI_py.plot`
(e.g. :mod:`SCI_py.metrics.drift` and :mod:`SCI_py.plot.drift`).

- :mod:`.completeness` -- data completeness per lab / condition
- :mod:`.cvs` -- coefficients of variation
- :mod:`.correlations` -- sample-to-sample correlations, missingness equalisation
- :mod:`.drift` -- signal drift over acquisition time
- :mod:`.rt_stability` -- retention-time stability
- :mod:`.dimensionality_reduction` -- PCA and other embeddings
- :mod:`.batch_correction` -- regress out covariates
- :mod:`.software_overlap` -- overlap of identifications between search engines
- :mod:`.software_versions` -- DIA-NN version and DIA-NN vs Spectronaut comparisons
"""

from .completeness import (
    compute_completeness,
    compute_combined_curve,
    store_completeness,
)
from .cvs import (
    compute_cv_matrix,
    filter_common_proteome,
)
from .correlations import (
    create_reference_set,
    equalise_plate_missingness,
    compute_pairwise_correlations,
    compute_pairwise_stats,
)
from .drift import (
    compute_time_norm,
    compute_drift,
    compute_rt_drift_correlations,
    fit_drift_stats,
)
from .rt_stability import (
    compute_rt_norm_layer,
    compute_rt_proportions,
)
from .dimensionality_reduction import (
    select_pca_features,
    compute_pca_embeddings,
    compute_embeddings,
)
from .batch_correction import (
    build_design_matrix,
    regress_out_covariates,
)
from .software_overlap import (
    intermediate_to_upsetplotdf,
    get_boxplot_data,
)
from .software_versions import (
    shared_obs,
    compute_overlap_decomposition,
    compute_binned_overlap,
)

__all__ = [
    "compute_completeness",
    "compute_combined_curve",
    "store_completeness",
    "compute_cv_matrix",
    "filter_common_proteome",
    "create_reference_set",
    "equalise_plate_missingness",
    "compute_pairwise_correlations",
    "compute_pairwise_stats",
    "compute_time_norm",
    "compute_drift",
    "compute_rt_drift_correlations",
    "fit_drift_stats",
    "compute_rt_norm_layer",
    "compute_rt_proportions",
    "select_pca_features",
    "compute_pca_embeddings",
    "compute_embeddings",
    "build_design_matrix",
    "regress_out_covariates",
    "intermediate_to_upsetplotdf",
    "get_boxplot_data",
    "shared_obs",
    "compute_overlap_decomposition",
    "compute_binned_overlap",
]
