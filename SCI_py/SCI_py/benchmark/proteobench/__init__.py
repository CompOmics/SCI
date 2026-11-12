"""Convert search-engine results into ProteoBench datapoints.

Most of this code is copied from ProteoBench
(https://github.com/Proteobench/ProteoBench) so the package does not depend on it.

- :mod:`.convert` -- :func:`parse_tables_to_proteobench`, the entry point
- :mod:`.input_formats` -- load each tool's output in ProteoBench's standard format
- :mod:`.settings` -- per-tool parse settings (column names, modifications)
- :mod:`.quantscores` -- quantification scores per precursor
- :mod:`.datapoint` -- summary metrics for one run (one point in the ProteoBench plot)
"""

from .convert import (
    parse_tables_to_proteobench,
    get_datapoint,
    PARSE_SETTINGS_DIR,
)
from .input_formats import load_input_file
from .settings import ParseSettingsBuilder
from .quantscores import QuantScores
from .datapoint import QuantDatapoint

__all__ = [
    "parse_tables_to_proteobench",
    "get_datapoint",
    "PARSE_SETTINGS_DIR",
    "load_input_file",
    "ParseSettingsBuilder",
    "QuantScores",
    "QuantDatapoint",
]
