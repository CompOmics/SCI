"""Build the per-sample metadata table (SDRF + cellenONE + manual annotations).

- :mod:`.read` -- read the SDRF, cellenONE and docx-annotation sources
- :mod:`.columns` -- clean, impute and type the merged table (:func:`parse_metadata_columns`)
- :mod:`.filenames` -- recognise the lab from a raw-file name and convert it to a sample name

Typical usage
-------------
>>> from SCI_py.metadata import read_SDRF, read_cellenONE, add_docx_annotations, parse_metadata_columns
>>> sdrf     = read_SDRF(ROOT_SDRF_IN)
>>> metadata = sdrf.merge(read_cellenONE(ROOT_CELLENONE_IN),
...                       on=['comment[X position]', 'comment[Y position]', 'comment[Plate Number]'],
...                       how='left')
>>> metadata = parse_metadata_columns(add_docx_annotations(metadata))
"""

from .read import (
    read_SDRF,
    read_cellenONE,
    add_docx_annotations,
)
from .columns import (
    xy_to_rc,
    parse_temporal_columns,
    compute_injection_order,
    impute_missing_values,
    convert_column_types,
    parse_metadata_columns,
)
from .filenames import (
    get_institute,
    FILENAME_CONVERTERS,
)

__all__ = [
    "read_SDRF",
    "read_cellenONE",
    "add_docx_annotations",
    "xy_to_rc",
    "parse_temporal_columns",
    "compute_injection_order",
    "impute_missing_values",
    "convert_column_types",
    "parse_metadata_columns",
    "get_institute",
    "FILENAME_CONVERTERS",
]
