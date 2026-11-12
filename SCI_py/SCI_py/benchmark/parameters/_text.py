"""Text helpers shared by the PEAKS and Spectronaut parameter parsers."""

import re
from typing import List, Optional


def clean_text(text: str) -> str:
    """
    Clean the input text by removing leading and trailing spaces, colons, commas, or tabs.

    Parameters
    ----------
    text : str
        The text to be cleaned.

    Returns
    -------
    str
        The cleaned text.
    """
    text = re.sub(r"^[\s:,\t]+|[\s:,\t]+$", "", text)
    return text


def extract_value_regex(lines: List[str], search_term: str) -> Optional[str]:
    """
    Extract the value associated with a search term using regular expressions.

    Parameters
    ----------
    lines : List[str]
        The list of lines to search through.
    search_term : str
        The regular expression to search for in the lines.

    Returns
    -------
    Optional[str]
        The extracted value, or None if the search term is not found.
    """
    return next((clean_text(re.split(search_term, line)[1]) for line in lines if re.search(search_term, line)), None)
