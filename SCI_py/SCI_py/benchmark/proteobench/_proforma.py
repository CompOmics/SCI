"""ProForma sequence helpers shared by the ProteoBench input-format and settings parsers.

Copied from ProteoBench (https://github.com/Proteobench/ProteoBench).
"""

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
