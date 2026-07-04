"""Helpers for parsing QC/Quant search-engine results into ProteoBench format.

This subpackage holds the logic previously defined inline in
``Notebooks/1_QCQuant/01_parse_results_clean.ipynb``:

- :mod:`fasta`    -- FASTA header parsing (:func:`parse_fasta_headers`)
- :mod:`readers`  -- per-engine result-folder readers
- :mod:`variants` -- ``main`` / ``ticnorm`` / ``firstpass`` variant logic
- :mod:`entries`  -- enumeration of all entries under the QC/Quant root
- :mod:`cache`    -- per-entry pickle cache (:class:`EntryCache`)
"""

from .fasta import parse_fasta_headers
from .readers import read_result_folder, read_diann_firstpass
from .variants import (
    apply_tic_norm,
    needs_ticnorm,
    has_firstpass,
    firstpass_comment,
)
from .entries import enumerate_entries
from .cache import EntryCache

__all__ = [
    'parse_fasta_headers',
    'read_result_folder',
    'read_diann_firstpass',
    'apply_tic_norm',
    'needs_ticnorm',
    'has_firstpass',
    'firstpass_comment',
    'enumerate_entries',
    'EntryCache',
]
