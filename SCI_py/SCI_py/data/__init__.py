"""Load search-engine reports and turn them into AnnData objects.

- :mod:`.reports` -- load DIA-NN and Spectronaut reports, with lab and sample name added
- :mod:`.anndata_utils` -- build, align and inspect AnnData objects
"""

from .reports import (
    read_spectronaut,
    load_spectronaut,
    load_diann,
)
from .anndata_utils import (
    df_to_anndata,
    align_obs_across_anndata,
    replace_zeros_with_nan,
    inspect_anndata_size,
    get_proteins_by_condition,
    simplify_protein_group_ids,
)

__all__ = [
    "read_spectronaut",
    "load_spectronaut",
    "load_diann",
    "df_to_anndata",
    "align_obs_across_anndata",
    "replace_zeros_with_nan",
    "inspect_anndata_size",
    "get_proteins_by_condition",
    "simplify_protein_group_ids",
]
