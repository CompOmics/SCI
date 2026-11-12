"""Fixed facts about the SCI study, shared by all other subpackages.

- :mod:`.labs` -- lab ids, their anonymised names (``Astral_1`` ... ``timsTOF_3``) and conditions
- :mod:`.samples` -- plate layout, QC/benchmark run names, SDRF value maps and imputations
- :mod:`.plotting` -- colours and display order for labs, conditions and software
"""

from .labs import (
    LAB_NAME_MAPPING,
    LAB_NAMES,
    CONDITIONS,
    LABS,
)
from .samples import (
    MAPPER_PATH,
    LAB_TO_PLATENR,
    QC_RUNS,
    N_ROWS,
    PLATE_ANNOTATIONS,
    PLATE_NOTES,
    INSTITUTE_PATTERNS,
    ORGANISM_MAP,
    ORGANISM_PART_MAP,
    CELL_TYPE_MAP,
    DISEASE_MAP,
    CELL_LINE_MAP,
    SOURCE_NAME_MAPS,
    SCALAR_IMPUTATIONS,
    SAMPLE_TYPE_A,
    SAMPLE_TYPE_B,
)
from .plotting import (
    CONDITION_COLORS,
    LAB_COLORS,
    CONDITION_ORDER,
    LAB_ORDER,
    SOFTWARE_COLORS,
)

__all__ = [
    "LAB_NAME_MAPPING",
    "LAB_NAMES",
    "CONDITIONS",
    "LABS",
    "MAPPER_PATH",
    "LAB_TO_PLATENR",
    "QC_RUNS",
    "N_ROWS",
    "PLATE_ANNOTATIONS",
    "PLATE_NOTES",
    "INSTITUTE_PATTERNS",
    "ORGANISM_MAP",
    "ORGANISM_PART_MAP",
    "CELL_TYPE_MAP",
    "DISEASE_MAP",
    "CELL_LINE_MAP",
    "SOURCE_NAME_MAPS",
    "SCALAR_IMPUTATIONS",
    "SAMPLE_TYPE_A",
    "SAMPLE_TYPE_B",
    "CONDITION_COLORS",
    "LAB_COLORS",
    "CONDITION_ORDER",
    "LAB_ORDER",
    "SOFTWARE_COLORS",
]
