import os
from glob import glob
import pandas as pd

from .parse import xy_to_rc, FILENAMECONVERTERS
from .constants import LAB_TO_PLATENR

def read_SDRF(root):
    sdrf_list = []

    for sdrf_path in glob(os.path.join(root, '*', 'SDRF.xlsx')):
        sdrf_df = pd.read_excel(sdrf_path)
        plate_id = os.path.basename(os.path.dirname(sdrf_path))
        sdrf_df['plate_id'] = plate_id
        sdrf_df['institute_name'] = plate_id.split('_')[0]
        sdrf_list.append(sdrf_df)
    sdrf = pd.concat(sdrf_list, axis=0)

    sdrf['sample_name'] = sdrf.apply(
        lambda x: FILENAMECONVERTERS[
            x['plate_id']
        ](
            x['comment[data file]']
        ).replace(
            'control',
            'Control'
        ),
        axis=1
    )
    sdrf['comment[Plate Number]'] = sdrf['plate_id'].apply(
        lambda x: LAB_TO_PLATENR[x]
    )
    sdrf['WellName'] = sdrf.apply(lambda x: xy_to_rc(x['comment[X position]'], x['comment[Y position]']), axis=1)
    XPos_mapper = {l: i+1 for i, l in enumerate('ABCDEFGHIJKLMNOP')}
    sdrf['comment[X position]'] = sdrf['comment[X position]'].apply(
        lambda x: XPos_mapper[x]
    )
    
    sdrf['Reinjection'] = False

    # Annotate Wyss reinjections
    reinjection_1 = pd.read_csv(os.path.join(root, 'Wyss_Astral_1', 'CSI_Plate4run1.csv'), skiprows=1)
    reinjection_2 = pd.read_csv(os.path.join(root, 'Wyss_Astral_1', 'CSI_Plate4run3.csv'), skiprows=1)
    reinjection_2['sample_name'] = reinjection_2['File Name'].apply(lambda x: "_".join(x.split('_')[6:-1]))
    wyss_reinjection = reinjection_2['sample_name'].tolist()
    sdrf['Reinjection'] = sdrf.apply(
        lambda x: True if x['plate_id']=='Wyss_Astral_1' and x['sample_name'] in wyss_reinjection else x['Reinjection'], 
        axis=1
    )

    # Filter out Broad institute 'redo'
    broad_reinjection = sdrf.loc[
        (sdrf['plate_id']=='Broad_timsTOF_1') &
        (sdrf['comment[data file]'].apply(lambda x: 'redo' in x)),
        'sample_name'
    ].tolist()
    sdrf['Reinjection'] = sdrf.apply(
        lambda x: True if x['plate_id']=='Broad_timsTOF_1' and x['sample_name'] in broad_reinjection else x['Reinjection'],
        axis=1
    )

    # Remove duplicates
    sdrf = (
        sdrf.sort_values('Reinjection', ascending=False)
        .drop_duplicates(subset=['sample_name', 'plate_id'], keep='first')
    )
    
    # BYU pipetted the QC in the blank wells and reversely so.
    mask = sdrf['plate_id'] == 'BYU_timsTOF_1'
    sdrf.loc[mask, 'source name'] = sdrf.loc[mask, 'source name'].apply(
        lambda x: 'Site specific QC' if x == 'Site specific blanks'
        else 'Site specific blanks' if x == 'Site specific QC'
        else x
    )
    
    # Parse more columns
    return sdrf

def read_cellenONE(folder):
    isolated_df_list = []
    for plate_path in glob(os.path.join(folder, '*')):
        plate_nr = int(os.path.basename(plate_path).split('_')[1])
        for sample_type_folder in glob(os.path.join(plate_path, '*')):
            sample_type = os.path.basename(sample_type_folder).split('_')[0]

            try:
                found_path = True
                isolated_path = list(glob(os.path.join(sample_type_folder, 'Reordered_*_isolated.xls')))[0]
            except:
                found_path = False
                isolated_path = '' # For linting checks
                print('Not Found\nPlate {}\tSample {}\n'.format(plate_nr, sample_type))
            
            if found_path:
                isolated_df = pd.read_csv(isolated_path, sep='\t')
                isolated_df['comment[Plate Number]'] = plate_nr
                isolated_df['comment[Sample Type]'] = sample_type
                isolated_df['image_path'] = isolated_df['ImageFile'].apply(lambda x: os.path.join(isolated_path, x.split('"')[1]))
    
                isolated_df_list.append(isolated_df)

    df = pd.concat(isolated_df_list, ignore_index=True)
    df = df.loc[:,['Unnamed' not in col for col in isolated_df.columns]]
    df = df.rename(columns={'XPos': 'comment[X position]', 'YPos': 'comment[Y position]'})
    df = df.loc[:,['Unnamed' not in col for col in df.columns]]
    return df

def add_docx_annotations(sdrf):
    # Plate 1
    plate1_notes    = "HEK dispensed instead of K562 in wells C13-N13"

    # Plate 2
    plate2_dying    = ['R19-C10', 'R9-C8']
    plate2_dividing = ['R9-C4', 'R10-C13']

    # Plate 3
    plate3_dying    = ['R16-C3', 'R19-C5', 'R6-C8', 'R10-C10', 'R11-C10', 'R12-C11', 'R9-C14']
    plate3_dividing = ['R16-C5', 'R14-C10', 'R6-C11']

    # Plate 4
    plate4_dying    = ['R5-C3', 'R4-C14', 'R10-C7', 'R10-C6', 'R4-C7', 'R9-C9']
    plate4_dividing = ['R14-C8', 'R18-C10', 'R6-C14']

    # Plate 5
    plate5_dying    = ['R9-C5', 'R7-C7', 'R5-C8']
    plate5_onecell  = ['R2-C11']  # one cell only

    # Plate 6
    plate6_dying    = ['R19-C5', 'R8-C12']
    plate6_dividing = ['R13-C6', 'R5-C13', 'R12-C14']

    # Plate 7
    plate7_dying    = ['R13-C3', 'R18-C6']
    plate7_dividing = ['R20-C10', 'R10-C11', 'R12-C12']

    # Plate 8
    plate8_dying    = ['R17-C4', 'R9-C3']
    plate8_dividing = ['R15-C12', 'R9-C9']
    plate8_doublet  = ['R20-C5']

    # Plate 9
    plate9_dying    = ['R7-C4', 'R6-C8', 'R11-C9']
    plate9_dividing = ['R13-C3', 'R8-C12', 'R12-C13']

    # Plate 10
    plate10_dying    = ['R11-C7']
    plate10_dividing = ['R20-C11', 'R6-C6', 'R4-C9']

    # Plate 11
    plate11_notes = "K562 less viable on that day (81.4%)"

    # Plate 12
    plate12_notes = "K562 less viable on that day (81.4%)"

    # Plate 13
    plate13_dying    = ['R17-C3', 'R14-C6', 'R21-C8', 'R13-C9']
    plate13_dividing = ['R16-C12']

    # Plate 14
    plate14_dying    = ['R7-C6', 'R6-C12']
    plate14_dividing = ['R4-C8']

    sdrf['Doublet'] = False
    sdrf['Dividing'] = False
    sdrf['Dying'] = False
    sdrf['CellInBlank'] = False

    def annotate(col, annotations, plate_id):
        if plate_id not in sdrf['comment[Plate Number]'].unique():
            return
        mask = (sdrf['comment[Plate Number]'] == plate_id) & (sdrf['WellName'].isin(annotations))
        sdrf.loc[mask, col] = True

    # Plate 2
    annotate('Dying',    plate2_dying,    2)
    annotate('Dividing', plate2_dividing, 2)

    # Plate 3
    annotate('Dying',    plate3_dying,    3)
    annotate('Dividing', plate3_dividing, 3)

    # Plate 4
    annotate('Dying',    plate4_dying,    4)
    annotate('Dividing', plate4_dividing, 4)

    # Plate 5
    annotate('Dying',       plate5_dying,   5)
    annotate('CellInBlank', plate5_onecell, 5)

    # Plate 6
    annotate('Dying',    plate6_dying,    6)
    annotate('Dividing', plate6_dividing, 6)

    # Plate 7
    annotate('Dying',    plate7_dying,    7)
    annotate('Dividing', plate7_dividing, 7)

    # Plate 8
    annotate('Dying',    plate8_dying,    8)
    annotate('Dividing', plate8_dividing, 8)
    annotate('Doublet',  plate8_doublet,  8)

    # Plate 9
    annotate('Dying',    plate9_dying,    9)
    annotate('Dividing', plate9_dividing, 9)

    # Plate 10
    annotate('Dying',    plate10_dying,    10)
    annotate('Dividing', plate10_dividing, 10)

    # Plate 13
    annotate('Dying',    plate13_dying,    13)
    annotate('Dividing', plate13_dividing, 13)

    # Plate 14
    annotate('Dying',    plate14_dying,    14)
    annotate('Dividing', plate14_dividing, 14)
    return sdrf

