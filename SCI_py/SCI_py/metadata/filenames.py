"""Map raw-file names to labs and to sample names.

Each lab named its raw files differently. :func:`get_institute` recognises the
lab (plate id, e.g. ``'IMP_Astral_1'``) from a filename, and
:data:`FILENAME_CONVERTERS` maps each plate id to a function that turns that
lab's raw-file name into the shared sample name used in the SDRF tables.

Example
-------
>>> plate_id = get_institute(run_name)
>>> sample_name = FILENAME_CONVERTERS[plate_id](run_name + '.RAW')
"""

from __future__ import annotations

from typing import Callable

from ..constants.samples import INSTITUTE_PATTERNS

__all__ = [
    "get_institute",
    "FILENAME_CONVERTERS",
]


# ---------------------------------------------------------------------------
# Filename utilities
# ---------------------------------------------------------------------------

def get_institute(filename: str) -> str | None:
    """
    Infer the canonical plate ID from a raw data filename.

    Parameters
    ----------
    filename:
        Raw filename as stored in the SDRF ``comment[data file]`` column.

    Returns
    -------
    str or None
        Canonical plate ID (e.g. ``"Wyss_Astral_1"``) or ``None`` if the
        filename cannot be matched.
    """
    for pattern, plate_id in INSTITUTE_PATTERNS.items():
        if pattern in filename:
            return plate_id
    return None


def _filename_default(filename: str) -> str:
    """Sample name for the default scheme (DTU_Astral_1, LSMBO_Astral_1): drop the first 6 '_' fields and the last one."""
    return "_".join(filename.split("_")[6:-1])


def _filename_Wyss(filename: str) -> str:
    """Sample name for Wyss_Astral_1 files: drop the first 6 '_' fields (and the last one unless the name ends in .RAW)."""
    if filename.upper().endswith(".RAW"):
        return "_".join(filename.split("_")[6:]).split(".")[0]
    result = "_".join(filename.split("_")[6:-1]).split(".")[0]
    return result.rstrip("_")


def _filename_KI(filename: str) -> str:
    """Sample name for KI_Astral_1/KI_Astral_2 files: rename '7_8min' to '140SPD', then drop the first 6 and last '_' fields."""
    return "_".join(filename.replace("7_8min", "140SPD").split("_")[6:-1])


def _filename_IMP_Astral(filename: str) -> str:
    """Sample name for IMP_Astral_1 files: drop the first 6 and the last 2 '_' fields."""
    return "_".join(filename.split("_")[6:-2])


def _filename_IMP_timsTOF(filename: str) -> str:
    """Sample name for IMP_timsTOF_1 files: drop the first 9 and the last 3 '_' fields."""
    return "_".join(filename.split("_")[9:-3])


def _filename_Broad(filename: str) -> str:
    """Sample name for Broad_timsTOF_1 files: drop the first 6 '_' fields and a trailing block that depends on '.RAW' and 'redo'."""
    has_redo = "redo" in filename
    if filename.upper().endswith(".RAW"):
        return "_".join(filename.split("_")[6: -5 if has_redo else -4])
    return "_".join(filename.split("_")[6: -2 if has_redo else -1])


def _filename_BYU(filename: str) -> str:
    """Sample name for BYU_timsTOF_1 files: drop the first 6 and the last 4 '_' fields."""
    return "_".join(filename.split("_")[6:-4])


#: Maps plate ID → filename normalisation function.
FILENAME_CONVERTERS: dict[str, Callable] = {
    "DTU_Astral_1"   : _filename_default,
    "LSMBO_Astral_1" : _filename_default,
    "KI_Astral_1"    : _filename_KI,
    "KI_Astral_2"    : _filename_KI,
    "Wyss_Astral_1"  : _filename_Wyss,
    "IMP_Astral_1"   : _filename_IMP_Astral,
    "IMP_timsTOF_1"  : _filename_IMP_timsTOF,
    "Broad_timsTOF_1": _filename_Broad,
    "BYU_timsTOF_1"  : _filename_BYU,
}
