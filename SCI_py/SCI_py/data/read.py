import pandas as pd
import os

from ..metadata.parse import FILENAME_CONVERTERS, get_institute

def read_spectronaut(input_csv: str) -> pd.DataFrame:
    """
    Load a Spectronaut output file.

    Parameters
    ----------
    input_csv : str
        The path to the Spectronaut output file.

    Returns
    -------
    pd.DataFrame
        The loaded dataframe.
    """
    input_data_frame = pd.read_csv(input_csv, low_memory=False, sep="\t")

    if input_data_frame["FG.Quantity"].dtype == object:
        try:
            input_csv.seek(0)
        except AttributeError:
            # if input_csv is a PathPosix object, it does not have a seek method
            # This can occur when the io util functions are used.
            # Should probably be fixed some way in the future
            pass
        input_data_frame = pd.read_csv(input_csv, low_memory=False, sep="\t", decimal=",")
    
    input_data_frame["FG.LabeledSequence"] = input_data_frame["FG.LabeledSequence"].str.strip("_")
    mapper_path = os.path.join(os.path.dirname(__file__), "mapper.csv")
    mapper_df = pd.read_csv(mapper_path).set_index("gene_name")
    mapper = mapper_df["description"].to_dict()
    input_data_frame["PG.ProteinGroups"] = input_data_frame["PG.ProteinGroups"].str.split(";")
    input_data_frame["PG.ProteinGroups"] = input_data_frame["PG.ProteinGroups"].map(
        lambda x: [mapper[protein] if protein in mapper.keys() else protein for protein in x]
    )
    input_data_frame["PG.ProteinGroups"] = input_data_frame["PG.ProteinGroups"].str.join(";")
    return input_data_frame

def load_spectronaut(path):
    report = read_spectronaut(
        input_csv=path
    )

    # Parse sample name and sample_ids
    report['plate_id'] = report['R.FileName'].progress_apply(
        get_institute
    )
    report['sample_name'] = report.progress_apply(
        lambda x: FILENAME_CONVERTERS[x['plate_id']](x['R.FileName']+'.RAW'),
        axis=1
    )
    report['sample_name'] = report['sample_name'].progress_apply(lambda x: x.lower())

    # Parse the modification annotation
    report['FG.PrecursorIon'] = report.progress_apply(
        lambda x: (
            str(x['FG.LabeledSequence'])
                .strip('_')
                .replace('[Oxidation (M)]', '[UNIMOD:35]')
                .replace('[Acetyl (Protein N-term)]', '[UNIMOD:1]')
                .replace('[Carbamidomethyl (C)]', '[UNIMOD:4]') + '/' + str(x['FG.Charge'])
        ),
        axis=1
    )
    return report

