import pandas as pd
from tqdm import tqdm
from glob import glob
import os
import numpy as np

from .proteobench_utils import _LOAD_FUNCTIONS
from ..constants.metadata import QC_RUNS


def read_data_diann(
    filename,
    converter,
    metadata
):
    report = pd.read_parquet(filename)
    report['path'] = report['Run']
    report['Run'] = report['Run'].apply(lambda x: converter(x+'.RAW'))

    metadata['Run'] = metadata['sample_name'].apply(lambda x: x.upper())
    metadata = metadata.rename(columns={'Intensity': 'CellenONE.Intensity'})
    
    report = report.merge(metadata, on='Run', how='left')
    
    return report

def read_data_msaid(
    foldername,
    converter,
    metadata,
    q_value_column='GLOBAL_Q_VALUE',
    q_value_value=.01,
    verbose=False
):
    report = _LOAD_FUNCTIONS['MSAID'](
        input_csv=os.path.join(foldername, 'precursors.tsv'),
        input_csv_secondary=os.path.join(foldername, 'fasta_mapping.tsv')
    )
    count_all = len(report)
    report['RAW_FILE_NAME'] = report['RAW_FILE_NAME'].apply(lambda x: converter(x+'.RAW'))

    if verbose:
        print('Filtering DECOYS: {} {:.2f}%'.format(
            report["DECOY"].sum(),
            report["DECOY"].sum() / count_all *100
        ))
    report = report[~report['DECOY']]

    if verbose:
        print('Filtering on qvalue {}: {} {:.2f}%'.format(
            q_value_column,
            (report[q_value_column]>q_value_value).sum(),
            (report[q_value_column]>q_value_value).sum() / count_all *100,
        ))

    report = report[report[q_value_column] < q_value_value]
    
    metadata['RAW_FILE_NAME'] = metadata['sample_name'].apply(lambda x: x.upper())
    metadata = metadata.rename(columns={'Intensity': 'CellenONE.Intensity'})
    
    report = report.merge(metadata, on='RAW_FILE_NAME', how='left')
    return report

def read_data_alphadia(
    filename,
    converter,
    metadata=None
):
    report = _LOAD_FUNCTIONS['AlphaDIA'](
        input_csv=filename
    )

    id_columns = ['sequence', 'mods', 'mod_sites', 'charge', 'genes', 'proforma']
    raw_file_col_mapper = {
        col:  converter(col+'.RAW') 
        for col in report.columns 
        if col not in id_columns
    }
    report = report.rename(columns=raw_file_col_mapper)
    return report

def read_data_spectronaut(
    foldername,
    converter,
    metadata
):
    result_path = list(glob(
        os.path.join(foldername, '*BGS Factory Report (Normal).tsv')
    ))
    if len(result_path) != 1:
        raise Exception(f'File path ambiguous. Check files at: {foldername}')
    result_path = result_path[0]

    report = _LOAD_FUNCTIONS['Spectronaut'](
        input_csv=result_path
    )
    
    report['R.FileName'] = report['R.FileName'].apply(lambda x: converter(x+'.RAW'))
    
    metadata['R.FileName'] = metadata['sample_name'].apply(lambda x: x.upper())
    metadata = metadata.rename(columns={'Intensity': 'CellenONE.Intensity'})

    report = report.merge(metadata, on='R.FileName', how='left')
    return report


# PEAKS
def map_filename(col, mapper):
    new_col = col
    for k in reversed(mapper): # Reverse to no map twice (first 11, then 1)
        if k in new_col:
            new_col = new_col.replace(k, mapper[k])
    if new_col == col:
        return {}
    return {col: new_col}

def create_filename_mapper(prefix, suffix, path_to_raw, converter):
    filename_mapper = {
        f'{prefix}{i+1}{suffix}': converter(
            os.path.basename(filename) + '.RAW'
            )+' ' for i, filename in enumerate(
            sorted(list(glob(os.path.join(path_to_raw, '*.raw'))))
        )
    }
    if filename_mapper == {}:
            filename_mapper = {
                f'{prefix}{i+1}{suffix}': converter(
                    os.path.basename(filename) + '.RAW'
                    )+' ' for i, filename in enumerate(
                    sorted(list(glob(os.path.join(path_to_raw, '*.d'))))
                )
            }
    return filename_mapper

def fix_column_names(df, mapper):
    column_mapper = {}
    for col in df.columns:
        column_mapper.update(
            map_filename(col, mapper)
        )
    df = df.rename(columns=column_mapper)
    return df

def read_data_peaks(
    filename,
    converter,
    path_to_raw,
    metadata=None,
    prefix='Sample ',
    suffix=' '

):
    mapper = create_filename_mapper(
        prefix=prefix,
        suffix=suffix,
        path_to_raw=path_to_raw,
        converter=converter
    )
    if mapper == {}:
        raise Exception('Mapper empty, check suffix and prefixes.')
    df = fix_column_names(
        df=_LOAD_FUNCTIONS['PEAKS'](filename),
        mapper=mapper
    )

    quant_cols = [i for i in df.columns if 'Normalized Area' in i]
    for run in QC_RUNS:
        if not any([run+' ' in col for col in quant_cols]):
            df[f'{run} Normalized Area'] = np.nan
    return df

def read_data_fragpipe_diannquant(
        foldername,
        converter,
        metadata
):
    def fill_empties(report):
        missing = [run for run in QC_RUNS if run not in report['Run'].unique()]
        if len(missing) == 0:
            return report
        new_dtps = []
        row = report.loc[0,:]
        for m in missing:
            new_row = row.to_dict()
            new_row['Run'] = m
            new_row['Precursor.Quantity'] = np.nan
            new_row['Precursor.Normalised'] = np.nan
            new_dtps.append(new_row)

        new_rep = pd.concat([report, pd.DataFrame(new_dtps)], axis=0, ignore_index=True)
        return new_rep

    path_to_results = os.path.join(
        foldername,
        'dia-quant-output',
        'report.tsv'
    )
    report = _LOAD_FUNCTIONS['FragPipe (DIA-NN quant)'](
        path_to_results
    )
    report['Run'] = report['Run'].apply(
        lambda x: converter(x.split('_uncalibrated')[0]+'.RAW')
    )
    report = fill_empties(report)

    metadata['Run'] = metadata['sample_name'].apply(lambda x: x.upper())
    metadata = metadata.rename(columns={'Intensity': 'CellenONE.Intensity'})
    
    report = report.merge(metadata, on='Run', how='left')
    return report