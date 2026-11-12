"""FASTA header parsing for the QC/Quant parsing pipeline."""

import re

# Matches UniProt-style headers such as ">sp|P00698|LYSC_CHICK ...",
# capturing the accession (group 1) and the entry name (group 2).
_HEADER_PATTERN = re.compile(r'^>\w+\|([^|]+)\|(\S+)')


def parse_fasta_headers(fasta_path):
    """Map accession to entry name from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Accession -> entry name, e.g. ``{'P00698': 'LYSC_CHICK'}``.

    Notes
    -----
    Used to backfill ``Protein.Names`` for DIA-NN second-pass reports,
    which do not always carry the entry name column.
    """
    mapping = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                m = _HEADER_PATTERN.match(line)
                if m:
                    mapping[m.group(1)] = m.group(2)
    return mapping
