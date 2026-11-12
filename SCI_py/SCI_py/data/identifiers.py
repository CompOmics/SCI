from __future__ import annotations


def simplify_protein_group_ids(protein_ids) -> list[str]:
    """
    Simplify UniProt-style protein(-group) identifiers of the form
    ``db|accession|entry_name`` (``;``-joined for protein groups) down to just the
    accession(s), e.g.::

        'sp|A0A075B6H7|KV37_HUMAN;sp|P01624|KV315_HUMAN' -> 'A0A075B6H7;P01624'

    A member that doesn't match the ``db|accession|entry_name`` pattern (exactly two
    ``|``) is left unchanged.

    Parameters
    ----------
    protein_ids:
        Iterable of protein(-group) identifier strings.

    Returns
    -------
    List of simplified identifier strings, same length and order as `protein_ids`.
    """
    simplified = []
    for group in protein_ids:
        accessions = []
        for member in group.split(';'):
            parts = member.split('|')
            accessions.append(parts[1] if len(parts) == 3 else member)
        simplified.append(';'.join(accessions))
    return simplified
