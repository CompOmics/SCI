"""Read search-engine result folders for the QC/Quant benchmark.

Two layers live here:

- **Per-engine report readers** (``read_data_diann``, ``read_data_spectronaut``,
  ``read_data_alphadia``, ``read_data_msaid``, ``read_data_peaks``,
  ``read_data_fragpipe_diannquant``) load one engine's output table and merge
  in the sample metadata.
- **Result-folder readers** (:func:`read_result_folder`,
  :func:`read_diann_firstpass`) pick the right per-engine reader and settings
  parser for one result folder and return a *table* dict with keys:

  - ``report``       : the parsed, filtered quantification table (DataFrame)
  - ``user_input``   : parsed search parameters (engine-specific object)
  - ``comment``      : human-readable settings label
  - ``input_format`` : ProteoBench input-format name (e.g. ``'DIA-NN'``)
  - ``lab_name``     : lab / plate identifier

  This table dict is the unit consumed by
  :func:`SCI_py.benchmark.proteobench.parse_tables_to_proteobench`.
"""

import os
from glob import glob

import numpy as np
import pandas as pd

from ..constants.samples import QC_RUNS
from ..metadata.filenames import FILENAME_CONVERTERS
from .constants import SEARCH_SETTINGS_NAMES
from .parameters.alphadia import extract_params as extract_params_alphadia
from .parameters.diann import extract_params as extract_params_diann
from .parameters.fragger import extract_params as extract_params_fragpipe
from .parameters.msaid import extract_params as extract_params_msaid
from .parameters.peaks import extract_params as extract_params_peaks
from .parameters.spectronaut import read_spectronaut_settings as extract_params_spectronaut
from .proteobench.input_formats import _LOAD_FUNCTIONS
from .variants import firstpass_comment


# ---------------------------------------------------------------------------
# Per-engine report readers
# ---------------------------------------------------------------------------

def read_data_diann(
    filename,
    converter,
    metadata
):
    """Read a DIA-NN ``report.parquet`` and merge in sample metadata.

    The original ``Run`` is kept in a ``path`` column, ``Run`` is converted to
    the sample name, and rows are left-merged with ``metadata`` on the
    upper-cased ``sample_name`` (metadata ``Intensity`` becomes
    ``CellenONE.Intensity``). No q-value or decoy filtering is done here.

    Parameters
    ----------
    filename : str
        Path to the DIA-NN parquet report.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    metadata : pandas.DataFrame
        Sample metadata with a ``sample_name`` column. A key column is added
        to it in place, so pass a copy.

    Returns
    -------
    pandas.DataFrame
        The DIA-NN report with metadata columns added.
    """
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
    """Read a CHIMERYS (MSAID) result folder, filter it and merge in sample metadata.

    Loads ``precursors.tsv`` with FASTA headers from ``fasta_mapping.tsv``,
    converts ``RAW_FILE_NAME`` to the sample name, removes decoys and keeps
    rows with ``q_value_column < q_value_value``, then left-merges with
    ``metadata`` on the upper-cased ``sample_name`` (metadata ``Intensity``
    becomes ``CellenONE.Intensity``).

    Parameters
    ----------
    foldername : str
        Folder containing ``precursors.tsv`` and ``fasta_mapping.tsv``.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    metadata : pandas.DataFrame
        Sample metadata with a ``sample_name`` column. A key column is added
        to it in place, so pass a copy.
    q_value_column : str, optional
        Q-value column to filter on, by default ``'GLOBAL_Q_VALUE'``.
    q_value_value : float, optional
        Q-value threshold (strictly below is kept), by default 0.01.
    verbose : bool, optional
        If True, print how many rows the decoy and q-value filters remove.

    Returns
    -------
    pandas.DataFrame
        Filtered precursor table with metadata columns added.
    """
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
    """Read an AlphaDIA precursor table and rename run columns to sample names.

    Every column other than the identifier columns (``sequence``, ``mods``,
    ``mod_sites``, ``charge``, ``genes``, ``proforma``) is treated as a run
    and renamed with ``converter``.

    Parameters
    ----------
    filename : str
        Path to the AlphaDIA output (e.g. ``precursors.parquet``).
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    metadata : pandas.DataFrame, optional
        Not used; accepted for a uniform reader signature.

    Returns
    -------
    pandas.DataFrame
        The AlphaDIA table with run columns named by sample.
    """
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
    """Read a Spectronaut result folder and merge in sample metadata.

    Loads the single ``*BGS Factory Report (Normal).tsv`` file in the folder,
    converts ``R.FileName`` to the sample name and left-merges with
    ``metadata`` on the upper-cased ``sample_name`` (metadata ``Intensity``
    becomes ``CellenONE.Intensity``).

    Parameters
    ----------
    foldername : str
        Spectronaut result folder.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    metadata : pandas.DataFrame
        Sample metadata with a ``sample_name`` column. A key column is added
        to it in place, so pass a copy.

    Returns
    -------
    pandas.DataFrame
        The Spectronaut report with metadata columns added.

    Raises
    ------
    Exception
        If the folder does not contain exactly one matching report file.
    """
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
    """Replace PEAKS sample labels in one column name using ``mapper``.

    Keys of ``mapper`` are checked in reverse order (so ``'Sample 11 '`` is
    replaced before ``'Sample 1 '``) and every key found in ``col`` is
    replaced by its value.

    Parameters
    ----------
    col : str
        Column name from a PEAKS table.
    mapper : dict
        PEAKS sample label -> sample name, as built by
        :func:`create_filename_mapper`.

    Returns
    -------
    dict
        ``{col: new_name}`` if anything was replaced, otherwise ``{}``.
    """
    new_col = col
    for k in reversed(mapper): # Reverse to no map twice (first 11, then 1)
        if k in new_col:
            new_col = new_col.replace(k, mapper[k])
    if new_col == col:
        return {}
    return {col: new_col}

def create_filename_mapper(prefix, suffix, path_to_raw, converter):
    """Build the PEAKS sample-label -> sample-name mapping from the raw files.

    PEAKS labels samples ``Sample 1``, ``Sample 2``, ... in the order of the
    raw files. The ``*.raw`` files in ``path_to_raw`` are sorted and numbered
    from 1; if there are none, ``*.d`` folders are used instead.

    Parameters
    ----------
    prefix, suffix : str
        Text around the sample number in PEAKS column names (e.g.
        ``'Sample '`` and ``' '``).
    path_to_raw : str
        Directory with the lab's ``.raw`` files or ``.d`` folders.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.

    Returns
    -------
    dict
        ``'{prefix}{n}{suffix}'`` -> sample name followed by a space. Empty
        if no raw files are found.
    """
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
    """Rename PEAKS sample-label columns to sample names.

    Parameters
    ----------
    df : pandas.DataFrame
        PEAKS table.
    mapper : dict
        PEAKS sample label -> sample name (see :func:`create_filename_mapper`).

    Returns
    -------
    pandas.DataFrame
        ``df`` with columns renamed via :func:`map_filename`.
    """
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
    """Read a PEAKS LFQ features table and name its sample columns by sample.

    Sample columns are renamed via :func:`create_filename_mapper`, and an
    all-NaN ``'<run> Normalized Area'`` column is added for every run in
    ``QC_RUNS`` that has no quantification column.

    Parameters
    ----------
    filename : str
        Path to the PEAKS ``lfq.features.csv``.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    path_to_raw : str
        Directory with the lab's raw files, used to number the samples.
    metadata : pandas.DataFrame, optional
        Not used; accepted for a uniform reader signature.
    prefix, suffix : str, optional
        Text around the sample number in PEAKS column names.

    Returns
    -------
    pandas.DataFrame
        The PEAKS table with sample-named columns.

    Raises
    ------
    Exception
        If no raw files are found, so the mapping is empty.
    """
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
    """Read a FragPipe (DIA-NN quant) result folder and merge in sample metadata.

    Loads ``dia-quant-output/report.tsv``, converts ``Run`` (with the
    ``_uncalibrated`` suffix removed) to the sample name, adds a placeholder
    row with NaN quantities for every run in ``QC_RUNS`` missing from the
    report, and left-merges with ``metadata`` on the upper-cased
    ``sample_name`` (metadata ``Intensity`` becomes ``CellenONE.Intensity``).

    Parameters
    ----------
    foldername : str
        FragPipe result folder.
    converter : callable
        Maps a raw-file name (with ``.RAW`` appended) to a sample name, e.g.
        ``FILENAME_CONVERTERS[lab]``.
    metadata : pandas.DataFrame
        Sample metadata with a ``sample_name`` column. A key column is added
        to it in place, so pass a copy.

    Returns
    -------
    pandas.DataFrame
        The DIA-NN-format report with metadata columns added.
    """
    def fill_empties(report):
        """Append a copy of the first row, with NaN quantities, for each QC run absent from ``report``."""
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

# ---------------------------------------------------------------------------
# Result-folder readers
# ---------------------------------------------------------------------------

def read_result_folder(result_folder, lab, engine, settings, *,
                       metadata, fasta_mapping, labname_to_raw):
    """Read one result folder and return a table dict.

    Parameters
    ----------
    result_folder : str
        Path to the ``<lab>/<engine>/<settings>`` folder.
    lab, engine, settings : str
        Entry identifiers, matching the folder structure.
    metadata : pandas.DataFrame
        Sample-level metadata; filtered to this lab via ``plate_id``.
    fasta_mapping : dict
        Accession -> entry name, used to backfill DIA-NN second-pass names.
    labname_to_raw : dict
        Lab -> directory of raw ``.raw``/``.d`` files (required by PEAKS).

    Returns
    -------
    dict or None
        Table dict, or ``None`` if the engine is not handled.
    """
    metadata_lab = metadata[metadata['plate_id'] == lab]
    converter = FILENAME_CONVERTERS[lab]
    comment = SEARCH_SETTINGS_NAMES[engine][settings]
    engine_label = engine

    if engine == 'AlphaDIA':
        report = read_data_alphadia(
            filename=os.path.join(result_folder, 'precursors.parquet'),
            converter=converter,
            metadata=metadata_lab.copy(),
        )
        user_input = extract_params_alphadia(os.path.join(result_folder, 'log.txt'))

    elif engine == 'Chimerys':
        report = read_data_msaid(
            foldername=result_folder,
            converter=converter,
            metadata=metadata_lab.copy(),
            q_value_column='GLOBAL_Q_VALUE',
            q_value_value=0.01,
            verbose=False,
        )
        user_input = extract_params_msaid(os.path.join(result_folder, 'experiment_details.csv'))
        engine_label = 'MSAID'

    elif engine == 'DIANN':
        report = read_data_diann(
            filename=os.path.join(result_folder, 'report.parquet'),
            converter=converter,
            metadata=metadata_lab.copy(),
        )
        report = report[(report['Global.Q.Value'] < 0.01) & (report['Decoy'] == 0)]
        report['Precursor.Quantity'] = report['Precursor.Normalised']
        if settings in ('4', '5'):  # 2nd-pass run: fill Protein.Names from FASTA
            report['Protein.Names'] = report['Protein.Group'].map(fasta_mapping)
        user_input = extract_params_diann(os.path.join(result_folder, 'report.log.txt'))
        engine_label = 'DIA-NN'

    elif engine == 'FragPipe':
        report = read_data_fragpipe_diannquant(
            foldername=result_folder,
            converter=converter,
            metadata=metadata.copy(),
        )
        report = report[(report['Global.Q.Value'] < 0.01) & (report['Decoy'] == 0)]
        report['Precursor.Quantity'] = report['Precursor.Normalised']
        with open(os.path.join(result_folder, 'fragpipe.workflow'), 'rb') as f:
            user_input = extract_params_fragpipe(f)
        engine_label = 'FragPipe (DIA-NN quant)'

    elif engine == 'PEAKS':
        report = read_data_peaks(
            filename=os.path.join(result_folder, 'lfq.features.csv'),
            converter=converter,
            path_to_raw=labname_to_raw[lab],
            metadata=metadata_lab.copy(),
        )
        user_input = extract_params_peaks(os.path.join(result_folder, 'parameters.txt'))

    elif engine == 'Spectronaut':
        report = read_data_spectronaut(
            foldername=result_folder,
            converter=converter,
            metadata=metadata_lab.copy(),
        )
        params_files = glob(os.path.join(result_folder, '*BGS Factory Settings.txt'))
        if len(params_files) != 1:
            raise ValueError(
                f'Expected 1 parameter file, found {len(params_files)}: {params_files}'
            )
        user_input = extract_params_spectronaut(params_files[0])

    else:
        return None

    return {
        'report': report,
        'user_input': user_input,
        'comment': comment,
        'input_format': engine_label,
        'lab_name': lab,
    }


def read_diann_firstpass(result_folder, lab, settings, *, metadata):
    """Read the DIA-NN first-pass report (``report-first-pass.parquet``).

    Parameters
    ----------
    result_folder : str
        Path to the DIA-NN ``<lab>/DIANN/<settings>`` folder.
    lab, settings : str
        Entry identifiers.
    metadata : pandas.DataFrame
        Sample-level metadata; filtered to this lab via ``plate_id``.

    Returns
    -------
    dict
        Table dict with a ``(1st pass)`` comment and ``DIA-NN`` input format.
    """
    metadata_lab = metadata[metadata['plate_id'] == lab]
    converter = FILENAME_CONVERTERS[lab]

    report = read_data_diann(
        filename=os.path.join(result_folder, 'report-first-pass.parquet'),
        converter=converter,
        metadata=metadata_lab.copy(),
    )
    report = report[(report['Global.Q.Value'] < 0.01) & (report['Decoy'] == 0)]
    report['Precursor.Quantity'] = report['Precursor.Normalised']
    user_input = extract_params_diann(os.path.join(result_folder, 'report.log.txt'))

    return {
        'report': report,
        'user_input': user_input,
        'comment': firstpass_comment(settings),
        'input_format': 'DIA-NN',
        'lab_name': lab,
    }
