"""Per-engine result-folder readers for the QC/Quant parsing pipeline.

Each reader turns one result folder into a *table* dict with keys:

- ``report``       : the parsed, filtered quantification table (DataFrame)
- ``user_input``   : parsed search parameters (engine-specific object)
- ``comment``      : human-readable settings label
- ``input_format`` : ProteoBench input-format name (e.g. ``'DIA-NN'``)
- ``lab_name``     : lab / plate identifier

This table dict is the unit consumed by ``parse_tables_to_proteobench``.
"""

import os
from glob import glob

from ..constants import search_settings_names
from ..read_quantqc import (
    read_data_msaid,
    read_data_diann,
    read_data_alphadia,
    read_data_spectronaut,
    read_data_peaks,
    read_data_fragpipe_diannquant,
)
from ..parameters.diann import extract_params as extract_params_diann
from ..parameters.spectronaut import read_spectronaut_settings as extract_params_spectronaut
from ..parameters.msaid import extract_params as extract_params_msaid
from ..parameters.alphadia import extract_params as extract_params_alphadia
from ..parameters.peaks import extract_params as extract_params_peaks
from ..parameters.fragger import extract_params as extract_params_fragpipe
from ...metadata.parse import FILENAME_CONVERTERS

from .variants import firstpass_comment


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
    comment = search_settings_names[engine][settings]
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
