import os

from glob import glob
from tqdm import tqdm

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pathlib

from .parse_settings import ParseSettingsBuilder
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