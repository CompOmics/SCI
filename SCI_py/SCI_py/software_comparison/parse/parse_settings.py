"""All formats available for the module."""

from __future__ import annotations

import os
from collections import defaultdict
from typing import Any, Dict, List

import pandas as pd
import toml
import re

def count_chars(input_string: str, isalpha: bool = True, isupper: bool = True) -> int:
    """
    Count the number of characters in the string that match the given criteria.

    Parameters
    ----------
    input_string : str
        The input string.
    isalpha : bool, optional
        Whether to count alphabetic characters. Defaults to True.
    isupper : bool, optional
        Whether to count uppercase characters. Defaults to True.

    Returns
    -------
    int
        The count of characters that match the criteria.
    """
    if isalpha and isupper:
        return sum(1 for char in input_string if char.isalpha() and char.isupper())
    if isalpha:
        return sum(1 for char in input_string if char.isalpha())
    if isupper:
        return sum(1 for char in input_string if char.isupper())
    
def get_stripped_seq(input_string: str, isalpha: bool = True, isupper: bool = True) -> str:
    """
    Get a stripped version of the sequence containing only characters that match the given criteria.

    Parameters
    ----------
    input_string : str
        The input string.
    isalpha : bool, optional
        Whether to include alphabetic characters. Defaults to True.
    isupper : bool, optional
        Whether to include uppercase characters. Defaults to True.

    Returns
    -------
    str
        The stripped sequence.
    """
    if isalpha and isupper:
        return "".join(char for char in input_string if char.isalpha() and char.isupper())
    if isalpha:
        return "".join(char for char in input_string if char.isalpha())
    if isupper:
        return "".join(char for char in input_string if char.isupper())


def match_brackets(
    input_string: str,
    pattern: str = r"\[([^]]+)\]",
    isalpha: bool = True,
    isupper: bool = True,
) -> tuple:
    """
    Match and extract bracketed modifications from the string.

    Parameters
    ----------
    input_string : str
        The input string.
    pattern : str, optional
        The regular expression pattern for matching modifications. Defaults to "\\[([^]]+)\\]".
    isalpha : bool, optional
        Whether to match alphabetic characters. Defaults to True.
    isupper : bool, optional
        Whether to match uppercase characters. Defaults to True.

    Returns
    -------
    tuple
        A tuple containing the matched modifications and their positions.
    """
    matches = [(match.group(), match.start(), match.end()) for match in re.finditer(pattern, input_string)]
    positions = (count_chars(input_string[0 : m[1]], isalpha=isalpha, isupper=isupper) for m in matches)
    mods = (m[0] for m in matches)
    return mods, positions


def to_lowercase(match) -> str:
    """
    Convert a match to lowercase.

    Parameters
    ----------
    match : re.Match
        The match object from a regular expression.

    Returns
    -------
    str
        The lowercase version of the matched string.
    """
    return match.group(0).lower()


def get_proforma_bracketed(
    input_string: str,
    before_aa: bool = True,
    isalpha: bool = True,
    isupper: bool = True,
    pattern: str = r"\[([^]]+)\]",
    modification_dict: dict = {
        "+57.0215": "Carbamidomethyl",
        "+15.9949": "Oxidation",
        "-17.026548": "Gln->pyro-Glu",
        "-18.010565": "Glu->pyro-Glu",
        "+42": "Acetyl",
    },
) -> str:
    """
    Generate a proforma string with bracketed modifications.

    Parameters
    ----------
    input_string : str
        The input sequence string.
    before_aa : bool, optional
        Whether to add the modification before the amino acid. Defaults to True.
    isalpha : bool, optional
        Whether to include alphabetic characters. Defaults to True.
    isupper : bool, optional
        Whether to include uppercase characters. Defaults to True.
    pattern : str, optional
        The regular expression pattern for matching modifications. Defaults to "\\[([^]]+)\\]".
    modification_dict : dict, optional
        A dictionary of modifications and their names.

    Returns
    -------
    str
        The proforma sequence with bracketed modifications.
    """
    input_string = re.sub(pattern, to_lowercase, input_string)
    modifications, positions = match_brackets(input_string, pattern=pattern, isalpha=isalpha, isupper=isupper)
    new_modifications = []

    for m in modifications:
        if m in modification_dict:
            new_modifications.append(modification_dict[m])
        else:
            new_modifications.append(m)

    modifications = new_modifications
    pos_mod_dict = dict(zip(positions, modifications))

    stripped_seq = get_stripped_seq(input_string, isalpha=isalpha, isupper=isupper)

    new_seq = ""
    for idx, aa in enumerate(stripped_seq):
        if before_aa:
            new_seq += aa
        if idx in pos_mod_dict:
            if idx == 0:
                new_seq += f"[{pos_mod_dict[idx]}]-"
            elif idx == len(stripped_seq):
                new_seq += f"-[{pos_mod_dict[idx]}]"
            else:
                new_seq += f"[{pos_mod_dict[idx]}]"
        if not before_aa:
            new_seq += aa

    return new_seq



class ParseSettingsBuilder:
    """
    Class to build the parser settings for a given input format.

    Parameters
    ----------
    parse_settings_dir : str
        The directory containing the parse settings files, by default None.
    module_id : str
        The ID of the module used to fetch the specific parse settings.
    """

    def __init__(self, parse_settings_dir: str, module_id: str):
        """
        Initialize the ParseSettingsBuilder object.

        Parameters
        ----------
        parse_settings_dir : str
            The directory containing the parse settings files.
        module_id : str
            The ID of the module used to fetch the specific parse settings.
        """

        self.PARSE_SETTINGS_TOMLS = toml.load(
            os.path.join(os.path.dirname(__file__), "io_parse_settings", "parse_settings_files.toml")
        )

        try:
            self.PARSE_SETTINGS_FILES = {
                key: os.path.join(parse_settings_dir, value)
                for key, value in self.PARSE_SETTINGS_TOMLS[module_id].items()
            }
        except KeyError:
            raise KeyError(
                f"Invalid module ID: {module_id}. Valid modules with configured parse settings are: {self.PARSE_SETTINGS_TOMLS.keys()}"
            )

        self.PARSE_SETTINGS_FILES_MODULE = os.path.join(parse_settings_dir, "module_settings.toml")
        self.INPUT_FORMATS = list(self.PARSE_SETTINGS_FILES.keys())
        self.MODULE_ID = module_id

        # Check if all files are present
        missing_files = [file for file in self.PARSE_SETTINGS_FILES.values() if not os.path.isfile(file)]
        if not os.path.isfile(self.PARSE_SETTINGS_FILES_MODULE):
            missing_files.append(self.PARSE_SETTINGS_FILES_MODULE)

        if missing_files:
            raise FileNotFoundError(f"The following parse settings files are missing: {missing_files}")

    def build_parser(self, input_format: str) -> object:
        """
        Build the parser for a given input format using the corresponding TOML files.

        Parameters
        ----------
        input_format : str
            The input format to build the parser for (e.g., "MaxQuant", "Sage").

        Returns
        -------
        ParseSettings
            The parser for the specified input format.
        """
        toml_file = self.PARSE_SETTINGS_FILES[input_format]
        parse_settings = toml.load(toml_file)
        parse_settings_module = toml.load(self.PARSE_SETTINGS_FILES_MODULE)

        parser = MODULE_TO_CLASS[self.MODULE_ID](parse_settings, parse_settings_module)
        if "modifications_parser" in parse_settings.keys():
            parser.add_modification_parser(ParseModificationSettings(parse_settings))

        return parser


class ParseSettingsQuant:
    """
    Structure that contains all the parameters used to parse
    the given benchmark run output depending on the software tool used.

    Parameters
    ----------
    parse_settings : Dict[str, Any]
        The settings for parsing, typically loaded from a TOML file.
    parse_settings_module : Dict[str, Any]
        Module-specific settings, typically loaded from a TOML file.
    """

    def __init__(self, parse_settings: Dict[str, Any], parse_settings_module: Dict[str, Any]):
        """
        Initialize the ParseSettings object with the parameters from the TOML files.

        Parameters
        ----------
        parse_settings : Dict[str, Any]
            The settings for parsing, typically loaded from a TOML file.
        parse_settings_module : Dict[str, Any]
            Module-specific settings, typically loaded from a TOML file.
        """
        self.mapper = parse_settings["mapper"]
        self.condition_mapper = parse_settings["condition_mapper"]
        self.run_mapper = parse_settings["run_mapper"]
        self.decoy_flag = parse_settings["general"]["decoy_flag"]
        self._species_dict = parse_settings["species_mapper"]
        self.contaminant_flag = parse_settings["general"]["contaminant_flag"]
        self.min_count_multispec = parse_settings_module["general"]["min_count_multispec"]
        self.analysis_level = parse_settings_module["general"]["level"]
        self._species_expected_ratio = parse_settings_module["species_expected_ratio"]
        self.modification_parser = None

    def species_dict(self) -> Dict[str, str]:
        """
        Get the species dictionary.

        Returns
        -------
        Dict[str, str]
            A dictionary of species mappings.
        """
        return self._species_dict

    def species_expected_ratio(self) -> float:
        """
        Get the expected species ratio.

        Returns
        -------
        float
            The expected ratio of species.
        """
        return self._species_expected_ratio

    def add_modification_parser(self, parser: ParseModificationSettings):
        """
        Add a modification parser to the settings.

        Parameters
        ----------
        parser : object
            The modification parser to add.
        """
        self.modification_parser = parser

    def _filter_decoys(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter out decoy proteins from the DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing protein data.

        Returns
        -------
        pd.DataFrame
            DataFrame with decoy proteins removed.
        """
        if "Reverse" in self.mapper:
            df_filtered = df[df["Reverse"] != self.decoy_flag].copy()
        else:
            df_filtered = df.copy()
        return df_filtered

    def _fix_colnames(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Fix column names in the DataFrame according to the mapper.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame with original column names.

        Returns
        -------
        pd.DataFrame
            DataFrame with standardized column names.
        """
        df.columns = [c.replace(".mzML.gz", ".mzML") for c in df.columns]
        return df

    def _mark_contaminants(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Mark contaminant proteins in the DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing protein data.

        Returns
        -------
        pd.DataFrame
            DataFrame with marked contaminants.
        """
        df["contaminant"] = df["Proteins"].str.contains(self.contaminant_flag)
        return df

    def _process_species_information(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate species information from the DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing protein data.

        Returns
        -------
        pd.DataFrame
            DataFrame with processed species information.
        """
        # Process species flags
        for flag, species in self._species_dict.items():
            df[species] = df["Proteins"].str.contains(flag)

        # Filter multi-species
        df["MULTI_SPEC"] = df[list(self._species_dict.values())].sum(axis=1) > self.min_count_multispec
        return df[df["MULTI_SPEC"] == False]

    def _validate_and_rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and rename columns according to the mapper.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame with original column names.

        Returns
        -------
        pd.DataFrame
            DataFrame with validated and renamed columns.
        """
        if not all(k in df.columns for k in self.mapper.keys()):
            raise ValueError(
                f"Columns {set(self.mapper.keys()).difference(set(df.columns))} not found in input dataframe."
                " Please check input file and selected software tool."
            )
        df.rename(columns=self.mapper, inplace=True)
        return df

    def _create_replicate_mapping(self) -> Dict[int, List[str]]:
        """
        Create replicate mapping for the DataFrame.

        Returns
        -------
        Dict[int, List[str]]
            A dictionary mapping replicates to raw data.
        """
        replicate_to_raw = defaultdict(list)
        for k, v in self.condition_mapper.items():
            replicate_to_raw[v].append(k)
        return replicate_to_raw

    def _handle_data_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Handle different data formats and convert to standard format.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame in either long or short format.

        Returns
        -------
        pd.DataFrame
            DataFrame converted to standard format.
        """
        # If "Raw file" is in mapper values, data is already in long format - skip melting
        if "Raw file" not in self.mapper.values():
            melt_vars = self.condition_mapper.keys()
            df_melted = df.melt(
                id_vars=list(set(df.columns).difference(set(melt_vars))),
                value_vars=melt_vars,
                var_name="Raw file",
                value_name="Intensity",
            )
        else:
            df_melted = df.copy()

        df_melted["replicate"] = df_melted["Raw file"].map(self.condition_mapper)
        return pd.concat([df_melted, pd.get_dummies(df_melted["Raw file"])], axis=1)

    def _process_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate modified sequences using the modification parser.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing modification data.

        Returns
        -------
        pd.DataFrame
            DataFrame with processed modifications.
        """
        if self.modification_parser is not None:
            return self.modification_parser.convert_to_proforma(df, self.analysis_level)
        return df

    def _filter_zero_intensities(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter out rows with zero or negative intensities.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing intensity data.

        Returns
        -------
        pd.DataFrame
            DataFrame with zero/negative intensities removed.
        """
        zero_intensity_count = len(df[df["Intensity"] <= 0])
        if zero_intensity_count > 0:
            pass
            # print(f"WARNING: {zero_intensity_count} rows with 0 intensity were removed.")
        return df[df["Intensity"] > 0]

    def _format_by_analysis_level(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Format DataFrame according to the analysis level.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame to be formatted.

        Returns
        -------
        pd.DataFrame
            Formatted DataFrame according to analysis level.
        """
        if self.analysis_level == "ion":
            if "proforma" in df.columns and "Charge" in df.columns:
                df["precursor ion"] = df["proforma"] + "/" + df["Charge"].astype(str)
            return df
        elif self.analysis_level == "peptidoform":
            if "proforma" in df.columns:
                df["peptidoform"] = df["proforma"]
            return df
        else:
            raise ValueError(f"Analysis level '{self.analysis_level}' not supported.")

    def convert_to_standard_format(self, df: pd.DataFrame) -> tuple[pd.DataFrame, Dict[int, List[str]]]:
        """
        Convert a software tool output into a generic format supported by the module.

        Steps:
        1. Validate and rename columns
        2. Create replicate mapping
        3. Filter decoys
        4. Fix column names
        5. Mark contaminants
        6. Process species information
        7. Handle data format (long vs short)
        8. Process modifications if needed
        9. Filter zero intensities
        10. Format based on analysis level

        Parameters
        ----------
        df : pd.DataFrame
            The input DataFrame to convert.

        Returns
        -------
        tuple[pd.DataFrame, Dict[int, List[str]]]
            The converted DataFrame and a dictionary mapping replicates to raw data.
        """
        df = self._validate_and_rename_columns(df)
        replicate_to_raw = self._create_replicate_mapping()
        df = self._filter_decoys(df)
        df = self._fix_colnames(df)
        df = self._mark_contaminants(df)
        df = self._process_species_information(df)
        df = self._process_modifications(df)
        df_melted = self._handle_data_format(df)
        df_melted = self._filter_zero_intensities(df_melted)
        return self._format_by_analysis_level(df_melted), replicate_to_raw


class ParseModificationSettings:
    """
    Settings for parsing modifications in protein data.

    Parameters
    ----------
    parse_settings : Dict[str, Any]
        Dictionary containing modification-specific parsing settings.
    """

    def __init__(self, parse_settings: Dict[str, Any]):
        """
        Initialize ParseModificationSettings.

        Parameters
        ----------
        parse_settings : Dict[str, Any]
            Dictionary containing modification-specific parsing settings.
        """
        self.modifications_mapper = parse_settings["modifications_parser"]["modification_dict"]
        self.modifications_isalpha = parse_settings["modifications_parser"]["isalpha"]
        self.modifications_isupper = parse_settings["modifications_parser"]["isupper"]
        self.modifications_before_aa = parse_settings["modifications_parser"]["before_aa"]
        self.modifications_pattern = parse_settings["modifications_parser"]["pattern"]
        self.modifications_pattern = rf"{self.modifications_pattern}"
        self.modifications_parse_column = parse_settings["modifications_parser"]["parse_column"]

    def convert_to_proforma(self, df: pd.DataFrame, analysis_level: str) -> tuple[pd.DataFrame, Dict[int, List[str]]]:
        """
        Convert modifications to ProForma format.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing modification data.
        analysis_level : str
            The level of analysis to perform.

        Returns
        -------
        tuple[pd.DataFrame, Dict[int, List[str]]]
            The converted DataFrame and a dictionary mapping replicates to raw data.
        """
        df["proforma"] = df[self.modifications_parse_column].apply(
            get_proforma_bracketed,
            before_aa=self.modifications_before_aa,
            isalpha=self.modifications_isalpha,
            isupper=self.modifications_isupper,
            pattern=self.modifications_pattern,
            modification_dict=self.modifications_mapper,
        )

        if analysis_level == "ion":
            try:
                df["precursor ion"] = df["proforma"] + "/" + df["Charge"].astype(str)
            except KeyError as e:
                raise KeyError(
                    "Not all columns required for making the precursor ion are available."
                    " Is the charge available in the input file?"
                ) from e

            return df

        elif analysis_level == "peptidoform":
            if "proforma" in df.columns:
                df["peptidoform"] = df["proforma"]
            return df


MODULE_TO_CLASS = {
    'SCI': ParseSettingsQuant,
}
