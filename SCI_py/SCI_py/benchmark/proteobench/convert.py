"""Turn benchmark table dicts into ProteoBench datapoints.

:func:`get_datapoint` converts one search-engine report into a ProteoBench
datapoint using the SCI parse settings in ``io_parse_settings/SCI``;
:func:`parse_tables_to_proteobench` does this for a list of table dicts (as
returned by :mod:`SCI_py.benchmark.readers`) and adds lab/plate labels.
"""


from tqdm import tqdm

import pandas as pd
import pathlib

from .settings import ParseSettingsBuilder
from .quantscores import QuantScores
from .datapoint import QuantDatapoint

PARSE_SETTINGS_DIR = (pathlib.Path(__file__).parent / "io_parse_settings" / "SCI").resolve().as_posix()


def get_datapoint(
    input_df: pd.DataFrame,
    user_input,
    module_id='SCI',
    input_format='DIA-NN',
    precursor_column_name='precursor ion',
    default_cutoff_min_prec=3,
    per_condition=True,
    all_datapoints=None
):
    """Compute a ProteoBench datapoint for one search-engine report.

    Converts the report to the ProteoBench standard format with the parse
    settings for ``input_format``, computes quantification scores and builds
    the datapoint.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Search-engine report (copied, not modified).
    user_input
        Parsed search parameters; its ``__dict__`` is stored in the datapoint.
    module_id : str, optional
        Parse-settings module id, by default ``'SCI'``.
    input_format : str, optional
        ProteoBench input-format name, by default ``'DIA-NN'``.
    precursor_column_name : str, optional
        Precursor column used by :class:`QuantScores`, by default ``'precursor ion'``.
    default_cutoff_min_prec : int, optional
        Minimum-precursor cutoff passed to
        :meth:`QuantDatapoint.generate_datapoint`, by default 3.
    per_condition : bool, optional
        Passed to :meth:`QuantDatapoint.generate_datapoint`, by default True.
    all_datapoints : optional
        Not used.

    Returns
    -------
    tuple
        ``(datapoint, intermediate)``: the datapoint (pandas.Series) and the
        intermediate table from :meth:`QuantScores.generate_intermediate`
        (pandas.DataFrame).

    Raises
    ------
    KeyError
        If any step (parse settings, format conversion, scoring, datapoint) fails.
    """
    input_df = input_df.copy(deep=True)
    try:
        parse_settings = ParseSettingsBuilder(
            parse_settings_dir=PARSE_SETTINGS_DIR,
            module_id=module_id
        ).build_parser(input_format)
    
    except Exception as e:
        raise KeyError(e)

    try:
        standard_format, replicate_to_raw = parse_settings.convert_to_standard_format(input_df)
    except KeyError as e:
        raise KeyError(e)

    # Calculate quantification scores
    try:
        quant_score = QuantScores(
            precursor_column_name, parse_settings.species_expected_ratio(), parse_settings.species_dict()
        )
    except Exception as e:
        raise KeyError(f"Error generating quantification scores: {e}")

    # Generate intermediate data structure
    try:
        intermediate_metric_structure = quant_score.generate_intermediate(standard_format, replicate_to_raw)
    except Exception as e:
        raise KeyError(f"Error generating intermediate data structure: {e}")

    # Generate current data point
    try:
        current_datapoint = QuantDatapoint.generate_datapoint(
            intermediate_metric_structure, 
            input_format, 
            user_input.__dict__,
            default_cutoff_min_prec=default_cutoff_min_prec,
            per_condition=per_condition
        )
    except Exception as e:
        raise KeyError(f"Error generating datapoint: {e}")
    
    return current_datapoint, intermediate_metric_structure

def parse_tables_to_proteobench(
    tables,
    lab_name_mapping,
    default_cutoff_min_prec=5,
):
    """Convert a list of table dicts into ProteoBench datapoints.

    Parameters
    ----------
    tables : list of dict
        Table dicts with keys ``report``, ``user_input``, ``input_format``,
        ``lab_name`` and ``comment`` (see :mod:`SCI_py.benchmark.readers`).
    lab_name_mapping : dict
        Plate id (``lab_name``) -> lab name, stored in the ``lab`` field.
    default_cutoff_min_prec : int, optional
        Passed to :func:`get_datapoint`, by default 5.

    Returns
    -------
    tuple
        ``(datapoints, intermediates)``: a DataFrame with one row per table,
        and the list of intermediate tables. Each gets ``plate_id``, ``lab``
        and ``comments`` added (intermediates also ``search_engine``).
    """
    datapoints, intermediates = [], []
    for table in tqdm(tables):
        datapoint, intermediate = get_datapoint(
            input_df=table['report'],
            user_input=table['user_input'],
            input_format=table['input_format'],
            default_cutoff_min_prec=default_cutoff_min_prec,
            all_datapoints=None
        )
        datapoint['plate_id'] = table['lab_name']
        datapoint['lab'] = lab_name_mapping[table['lab_name']]
        datapoint['comments'] = table['comment']

        intermediate['plate_id'] = table['lab_name']
        intermediate['lab'] = lab_name_mapping[table['lab_name']]
        intermediate['comments'] = table['comment']
        intermediate['search_engine'] = table['input_format']

        datapoints.append(datapoint)
        intermediates.append(intermediate)
    
    datapoints = pd.concat(datapoints, axis=1).T
    intermediate = pd.concat(intermediates, axis=1, ignore_index=True)
    return datapoints, intermediates   