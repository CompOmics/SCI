LABS = [
    'BYU_timsTOF_1',
    'Broad_timsTOF_1',
    'KI_Astral_1',
    'LSMBO_Astral_1',
    'KI_Astral_2',
    'DTU_Astral_1',
    'IMP_Astral_1',
    'Wyss_Astral_1',
    'IMP_timsTOF_1'
]
LAB_NAME_MAPPING = {
    'DTU_Astral_1': 'Astral_1',
    'LSMBO_Astral_1': 'Astral_2',
    'KI_Astral_1': 'Astral_3',
    'Wyss_Astral_1': 'Astral_4',
    'IMP_Astral_1': 'Astral_5',
    'KI_Astral_2': 'Astral_6',
    'IMP_timsTOF_1': 'timsTOF_1',
    'Broad_timsTOF_1': 'timsTOF_2',
    'BYU_timsTOF_1': 'timsTOF_3'
}

LAB_TO_PLATENR = {
    'BYU_timsTOF_1': '?', # Unsure
    'Broad_timsTOF_1': 8,
    'KI_Astral_1': 10,
    'LSMBO_Astral_1': 13,
    'KI_Astral_2': 11,
    'DTU_Astral_1': 14,
    'IMP_Astral_1': 5,
    'Wyss_Astral_1': 4,
    'IMP_timsTOF_1': 6
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
N_ROWS = 24   # to