from .constants import _X_TO_COL, N_ROWS, LAB_NAME_MAPPING
import pandas as pd
import numpy as np

########################
### FILENAME PARSING ###
########################
def filename_DTU(filename):
    return '_'.join(filename.split('_')[6:-1])

def filename_LSMBO(filename):
    return '_'.join(filename.split('_')[6:-1])

def filename_KI(filename):
    filename = filename.replace('7_8min', '140SPD')
    return '_'.join(filename.split('_')[6:-1])

def filename_Wyss(filename):
    if filename.endswith('.RAW'):
        return '_'.join(filename.split('_')[6:]).split('.')[0]

    filename = '_'.join(filename.split('_')[6:-1]).split('.')[0]
    if filename.endswith('_'):
        return filename[:-1]
    return filename

def filename_Wyss_sdrf(filename):
    return '_'.join(filename.split('_')[6:-1]).split('.')[0]

def filename_IMP_Astral(filename):
    return '_'.join(filename.split('_')[6:-2])

def filename_IMP_timsTOF(filename):
    return '_'.join(filename.split('_')[9:-3])

def filename_Broad(filename):
    if filename.endswith('.RAW'):
        if 'redo' in filename:
            return '_'.join(filename.split('_')[6:-5])
        return '_'.join(filename.split('_')[6:-4])

    if 'redo' in filename:
        return '_'.join(filename.split('_')[6:-2])
    return '_'.join(filename.split('_')[6:-1])

def filename_BYU(filename):
    return '_'.join(filename.split('_')[6:-4])

def get_institute(filename):
    if 'Wyss' in filename:
        return 'Wyss_Astral_1'
    if 'Broad' in filename:
        return 'Broad_timsTOF_1'
    if 'Schoof' in filename:
        return 'DTU_Astral_1'
    if 'LSMBO' in filename:
        return 'LSMBO_Astral_1'
    if 'RM_Astral' in filename:
        return 'IMP_Astral_1'
    if 'RM_timsTOF' in filename:
        return 'IMP_timsTOF_1'
    if 'BYU' in filename:
        return 'BYU_timsTOF_1'
    if '20251007_KI_AV' in filename:
        return 'KI_Astral_1'
    if '20260113_KI_AV' in filename:
        return 'KI_Astral_2'
    
FILENAMECONVERTERS = {
    'DTU_Astral_1': filename_DTU,
    'LSMBO_Astral_1': filename_LSMBO,
    'KI_Astral_1': filename_KI,
    'KI_Astral_2': filename_KI,
    'Wyss_Astral_1': filename_Wyss,
    'IMP_Astral_1': filename_IMP_Astral,
    'IMP_timsTOF_1': filename_IMP_timsTOF,
    'Broad_timsTOF_1': filename_Broad,
    'BYU_timsTOF_1': filename_BYU
}

##############################
### Plate metadata parsing ###
##############################

def xy_to_rc(x: str, y: int) -> str:
    """
    Convert physical (X, Y) → table (row_label, col_label).

    Parameters
    ----------
    x : str   Column letter, e.g. 'A' or 'a'
    y : int   Row number 1–24 (1 = bottom of plate)

    Returns
    -------
    (row_label, col_label) e.g. ('R24', 'C1')
    """
    x = x.upper()
    if x not in _X_TO_COL:
        raise ValueError(f"X must be a letter A–P, got '{x}'")
    if not (1 <= y <= N_ROWS):
        raise ValueError(f"Y must be 1–{N_ROWS}, got {y}")

    col_idx = _X_TO_COL[x]
    row_idx = N_ROWS + 1 - y          # Y=1 → R24,  Y=24 → R1

    return f"R{row_idx}-C{col_idx}"

def compute_injection_order(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Compute normalised injection order [0, 1] within each plate,
    based on comment[Gradient time]_datetime. Mirrors the R implementation.
    """
    df = metadata.copy()
    df['_datetime'] = pd.to_datetime(df['comment[Gradient time]_datetime'],
                                     errors='coerce')

    def _plate_order(group):
        dt           = group['_datetime']
        n_observed   = dt.notna().sum()
        ranks        = dt.rank(method='average', na_option='keep')
        if n_observed > 1:
            order = (ranks - 1) / (n_observed - 1)
        else:
            order = ranks * np.nan  # single cell → will be imputed to 0.5
        # NaN from rank → 0.5
        order = order.where(order.notna(), 0.5)
        # Missing datetime → plate median, fallback 0.5
        plate_median = order[dt.notna()].median()
        plate_median = plate_median if pd.notna(plate_median) else 0.5
        order        = order.where(dt.notna(), plate_median)
        return order

    df['injection_order'] = (df.groupby('plate_id', group_keys=False)
                               .apply(_plate_order))
    df['injection_order'] = df['injection_order'].fillna(0.5)
    df = df.drop(columns=['_datetime'])
    return df

def parse_metadata_columns(sdrf):
    # Sample name to lower case to prevent sample name mismatches due to differing capitalizations
    sdrf['sample_name'] = sdrf['sample_name'].apply(lambda x: x.lower())
    
    # Time related metadata parsing
    parsed = pd.to_datetime(sdrf['comment[plate arrival]'], format='%Y%m%d-%H%M%S', exact=False)
    sdrf['comment[plate arrival]_date'] = parsed.dt.floor('D')  # date only, as datetime
    sdrf['comment[plate arrival]_time'] = parsed.dt.hour + parsed.dt.minute / 60 + parsed.dt.second / 3600  # decimal hours
    sdrf['comment[plate arrival]_datetime'] = parsed

    parsed = pd.to_datetime(sdrf['comment[Gradient time]'], format='%Y%m%d-%H%M%S', exact=False)
    sdrf['comment[Gradient time]_date'] = parsed.dt.floor('D')  # date only, as datetime
    sdrf['comment[Gradient time]_time'] = parsed.dt.hour + parsed.dt.minute / 60 + parsed.dt.second / 3600  # decimal hours
    sdrf['comment[Gradient time]_datetime'] = parsed
    
    sdrf['instrument_type'] = sdrf['plate_id'].apply(lambda x: x.split('_')[1])
    sdrf['sample_id'] = sdrf.apply(lambda x: x['sample_name'] + '|' + x['plate_id'], axis=1)
    sdrf['source name'] = sdrf['source name'].apply(lambda x: x.replace('Cell A', 'K562').replace('Cell B', 'HEK293T'))
    sdrf['plate_idx'] = sdrf['plate_id'].apply(lambda x: LAB_NAME_MAPPING[x])
    
    # Parse injection order
    sdrf = compute_injection_order(sdrf)
    
    _ = sdrf.pop('comment[plate arrival].1')
    return sdrf

