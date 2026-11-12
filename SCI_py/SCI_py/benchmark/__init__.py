"""Search-engine comparison on the HYE QC/Quant samples (the "QC/Quant benchmark").

Pipeline, as run in ``Notebooks/1_QCQuant/01_parse_results.ipynb``:

1. :func:`enumerate_entries` -- list every lab / engine / settings result folder
2. :func:`read_result_folder` -- read one folder into a *table* dict (see :mod:`.readers`)
3. :func:`parse_tables_to_proteobench` -- turn tables into ProteoBench datapoints
4. :class:`EntryCache` -- cache each entry on disk so reruns are fast

Modules
-------
- :mod:`.readers` -- per-engine readers and result-folder readers
- :mod:`.entries`, :mod:`.variants` -- which result folders exist and their variants
  (``main`` / ``ticnorm`` / ``firstpass``)
- :mod:`.fasta` -- FASTA header parsing
- :mod:`.cache` -- per-entry pickle cache
- :mod:`.filters` -- row filters for the ProteoBench results table
- :mod:`.constants` -- settings labels per engine
- :mod:`.parameters` -- parse each engine's settings/parameter files
- :mod:`.proteobench` -- ProteoBench conversion (code copied from ProteoBench)
"""

from .entries import enumerate_entries
from .readers import (
    read_result_folder,
    read_diann_firstpass,
)
from .variants import (
    firstpass_comment,
    needs_ticnorm,
    has_firstpass,
    apply_tic_norm,
)
from .fasta import parse_fasta_headers
from .cache import EntryCache
from .filters import (
    filter_recommended,
    apply_version_filter,
)
from .constants import SEARCH_SETTINGS_NAMES
from .proteobench import parse_tables_to_proteobench

__all__ = [
    "enumerate_entries",
    "read_result_folder",
    "read_diann_firstpass",
    "firstpass_comment",
    "needs_ticnorm",
    "has_firstpass",
    "apply_tic_norm",
    "parse_fasta_headers",
    "EntryCache",
    "filter_recommended",
    "apply_version_filter",
    "SEARCH_SETTINGS_NAMES",
    "parse_tables_to_proteobench",
]
