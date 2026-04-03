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
LAB_NAMES = list(LAB_NAME_MAPPING.keys())

CONDITIONS = [
    'K562', 
    'HEK293T', 
    'QC - Hela',
    'Site specific QC', 
    'Sorting control',
    'Processing control',
    'Cellenone blank',
    'Site specific blanks',
]

LABS = [
    'Astral_1',
    'Astral_2',
    'Astral_3',
    'Astral_4',
    'Astral_5',
    'Astral_6',
    'timsTOF_1',
    'timsTOF_2',
    'timsTOF_3'
]