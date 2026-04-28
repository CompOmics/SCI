import os

MAPPER_PATH = os.path.join(os.path.dirname(__file__), "mapper.csv")

LAB_TO_PLATENR = {
    'BYU_timsTOF_1'  : -1,  # Unsure
    'Broad_timsTOF_1': 8,
    'KI_Astral_1'    : 10,
    'LSMBO_Astral_1' : 13,
    'KI_Astral_2'    : 11,
    'DTU_Astral_1'   : 14,
    'IMP_Astral_1'   : 5,
    'Wyss_Astral_1'  : 4,
    'IMP_timsTOF_1'  : 6,
}

QC_RUNS = [
    'HYEA_QC_B13',
    'HYEA_QC_B14',
    'HYEA_QC_B15',
    'HYEA_QC_B16',
    'HYEA_QC_B17',
    'HYEA_QC_B18',
    'HYEA_QC_B19',
    'HYEA_QC_B20',
    'HYEA_QC_B21',
    'HYEA_QC_B22',
    'HYEA_QC_B23',
    'HYEB_QC_O10',
    'HYEB_QC_O11',
    'HYEB_QC_O12',
    'HYEB_QC_O2',
    'HYEB_QC_O3',
    'HYEB_QC_O4',
    'HYEB_QC_O5',
    'HYEB_QC_O6',
    'HYEB_QC_O7',
    'HYEB_QC_O8',
    'HYEB_QC_O9',
]

_LETTERS = "ABCDEFGHIJKLMNOP"          # 16 columns
_X_TO_COL = {letter: i + 1 for i, letter in enumerate(_LETTERS)}   # A→1 … P→16
_COL_TO_X = {v: k for k, v in _X_TO_COL.items()}                   # 1→A … 16→P
N_ROWS = 24   # rows per plate (Y axis)


# ---------------------------------------------------------------------------
# Manual QC annotations
# ---------------------------------------------------------------------------

# Well-level annotations extracted from plate preparation documentation.
# Keys are plate numbers (int); values are dicts mapping QC flag → list of
# well labels in R{row}-C{col} notation.

PLATE_ANNOTATIONS: dict[int, dict[str, list[str]]] = {
    2: {
        "Dying"   : ["R19-C10", "R9-C8"],
        "Dividing": ["R9-C4",   "R10-C13"],
    },
    3: {
        "Dying"   : ["R16-C3", "R19-C5", "R6-C8",  "R10-C10",
                     "R11-C10", "R12-C11", "R9-C14"],
        "Dividing": ["R16-C5", "R14-C10", "R6-C11"],
    },
    4: {
        "Dying"   : ["R5-C3", "R4-C14", "R10-C7", "R10-C6", "R4-C7", "R9-C9"],
        "Dividing": ["R14-C8", "R18-C10", "R6-C14"],
    },
    5: {
        "Dying"      : ["R9-C5", "R7-C7", "R5-C8"],
        "CellInBlank": ["R2-C11"],
    },
    6: {
        "Dying"   : ["R19-C5", "R8-C12"],
        "Dividing": ["R13-C6", "R5-C13", "R12-C14"],
    },
    7: {
        "Dying"   : ["R13-C3", "R18-C6"],
        "Dividing": ["R20-C10", "R10-C11", "R12-C12"],
    },
    8: {
        "Dying"   : ["R17-C4", "R9-C3"],
        "Dividing": ["R15-C12", "R9-C9"],
        "Doublet" : ["R20-C5"],
    },
    9: {
        "Dying"   : ["R7-C4", "R6-C8",  "R11-C9"],
        "Dividing": ["R13-C3", "R8-C12", "R12-C13"],
    },
    10: {
        "Dying"   : ["R11-C7"],
        "Dividing": ["R20-C11", "R6-C6", "R4-C9"],
    },
    13: {
        "Dying"   : ["R17-C3", "R14-C6", "R21-C8", "R13-C9"],
        "Dividing": ["R16-C12"],
    },
    14: {
        "Dying"   : ["R7-C6", "R6-C12"],
        "Dividing": ["R4-C8"],
    },
}

# Free-text notes per plate (not encoded as per-well flags).
PLATE_NOTES: dict[int, str] = {
    1 : "HEK dispensed instead of K562 in wells C13-N13",
    11: "K562 less viable on that day (81.4%)",
    12: "K562 less viable on that day (81.4%)",
}


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

#: Maps institute identifier strings in raw filenames → canonical plate IDs.
INSTITUTE_PATTERNS: dict[str, str] = {
    "Wyss"          : "Wyss_Astral_1",
    "Broad"         : "Broad_timsTOF_1",
    "Schoof"        : "DTU_Astral_1",
    "LSMBO"         : "LSMBO_Astral_1",
    "RM_Astral"     : "IMP_Astral_1",
    "RM_timsTOF"    : "IMP_timsTOF_1",
    "BYU"           : "BYU_timsTOF_1",
    "20251007_KI_AV": "KI_Astral_1",
    "20260113_KI_AV": "KI_Astral_2",
}


# ---------------------------------------------------------------------------
# Source-name lookup tables
# ---------------------------------------------------------------------------

#: Maps source name → characteristics[organism].
ORGANISM_MAP: dict[str, str] = {
    "QC - HYE A"       : "Mixed",
    "QC - HYE B"       : "Mixed",
    "QC - Hela"        : "Homo sapiens",
    "Site specific QC" : "Homo sapiens",
    "K562"             : "Homo sapiens",
    "HEK293T"          : "Homo sapiens",
}

#: Maps source name → characteristics[organism part].
ORGANISM_PART_MAP: dict[str, str] = {
    "QC - Hela"        : "cervix",
    "Site specific QC" : "cervix",
    "K562"             : "bone marrow",
    "HEK293T"          : "kidney",
}

#: Maps source name → characteristics[cell type].
CELL_TYPE_MAP: dict[str, str] = {
    "QC - Hela"        : "malignant epithelial cell",
    "Site specific QC" : "malignant epithelial cell",
    "K562"             : "lymphoblast cells",
    "HEK293T"          : "embryonic kidney cells",
}

#: Maps source name → characteristics[disease].
DISEASE_MAP: dict[str, str] = {
    "QC - Hela"        : "epidermoid carcinoma",
    "Site specific QC" : "epidermoid carcinoma",
    "K562"             : "CML",
    "HEK293T"          : "healthy",
}

#: Maps source name → characteristics[cell line].
CELL_LINE_MAP: dict[str, str] = {
    "QC - Hela"        : "Hela",
    "Site specific QC" : "Hela",
    "K562"             : "K562",
    "HEK293T"          : "HEK293",
}

#: All columns derived from ``source name`` via lookup.
SOURCE_NAME_MAPS: dict[str, dict[str, str]] = {
    "characteristics[organism]"      : ORGANISM_MAP,
    "characteristics[organism part]" : ORGANISM_PART_MAP,
    "characteristics[cell type]"     : CELL_TYPE_MAP,
    "characteristics[disease]"       : DISEASE_MAP,
    "characteristics[cell line]"     : CELL_LINE_MAP,
}


# ---------------------------------------------------------------------------
# Missing value imputation
# ---------------------------------------------------------------------------

#: Scalar fill values for missing metadata.
#: Columns absent from the table are silently skipped.
SCALAR_IMPUTATIONS: dict[str, str | int | float] = {
    # SDRF controlled vocabulary
    "characteristics[pooled sample]"        : "not applicable",
    "comment[LC column]"                    : "not available",
    "comment[Sample Type]"                  : "not applicable",
    "image_path"                            : "not applicable",
    "Background"                            : "not applicable",
    "ImageFile"                             : "not applicable",
    "Teg"                                   : "not applicable",
    "Well"                                  : "Not applicable",
    # Numeric sentinels
    "characteristics[injection]"            : -1,
    "comment[technical replicate]"          : -1,
    "comment[fraction identifier]"          : -1,
    "characteristics[biological replicate]" : -1,
    "DropNo"                                : -1.0,
    "X"                                     : -1.0,
    "Y"                                     : -1.0,
    "Diameter"                              : -1.0,
    "Elongation"                            : -1.0,
    "Circularity"                           : -1.0,
    "Intensity"                             : -1.0,
    "Plate"                                 : -1.0,
    "Target"                                : -1.0,
    "Field"                                 : -1.0,
    "EjBound"                               : -1,
    "SedBound"                              : -1,
    "DetDiaMinTrans"                        : -1,
    "DetDiaMaxTrans"                        : -1,
    "IsoDiaMinTrans"                        : -1,
    "IsoDiaMaxTrans"                        : -1,
    "IsoIntMinTrans"                        : -1,
    "IsoIntMaxTrans"                        : -1,
    "IsoDiaMinFlu"                          : -1,
    "IsoDiaMaxFlu"                          : -1,
    "IsoIntMinFlu"                          : -1,
    "IsoIntMaxFlu"                          : -1,
    # Date / time sentinels
    "Date"                                  : "01/01/1900",
    "Time"                                  : "00:00:00",
}

# ---------------------------------------------------------------------------
# Column renaming
# ---------------------------------------------------------------------------

