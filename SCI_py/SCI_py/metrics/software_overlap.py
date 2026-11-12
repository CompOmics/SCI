"""Computations for the overlap between search engines (QC/Quant benchmark).

Turns ProteoBench intermediate tables into the inputs for the UpSet and box
plots in :mod:`SCI_py.plot.software_overlap`.

- :func:`intermediate_to_upsetplotdf` -- precursor-by-tool presence table for UpSet plots
- :func:`get_boxplot_data` -- values for the precursors in one UpSet intersection
"""
import pandas as pd

from ..constants.samples import SAMPLE_TYPE_A, SAMPLE_TYPE_B

_ALL_SAMPLES = SAMPLE_TYPE_A + SAMPLE_TYPE_B
_KEEP_COLS = [
    'log2_A_vs_B', 'nr_observed', 'species',
    'precursor ion', 'lab', 'search_engine', 'AbsQuantError',
]


def intermediate_to_upsetplotdf(intermediate, min_prec=5, per_condition=False):
    """Filter to precursors observed in >= min_prec samples; return UpSet df and filtered df.

    Parameters
    ----------
    intermediate:
        Long-format DataFrame for a single lab with precursor-level quantification.
    min_prec:
        Minimum number of non-missing observations required.
    per_condition:
        If True, require min_prec in both condition A and condition B separately.
    """
    sample_cols = [c for c in intermediate.columns if c in _ALL_SAMPLES]
    if sample_cols:
        if per_condition:
            a_cols = [c for c in intermediate.columns if c in SAMPLE_TYPE_A]
            b_cols = [c for c in intermediate.columns if c in SAMPLE_TYPE_B]
            mask = (
                intermediate[a_cols].notna().sum(axis=1) >= min_prec
            ) & (
                intermediate[b_cols].notna().sum(axis=1) >= min_prec
            )
        else:
            mask = intermediate[sample_cols].notna().sum(axis=1) >= min_prec
    else:
        mask = intermediate['nr_observed'] >= min_prec

    keep     = [c for c in _KEEP_COLS if c in intermediate.columns]
    filtered = (
        intermediate.loc[mask, keep]
        .reset_index(drop=True)
        .dropna(subset=['precursor ion', 'search_engine', 'AbsQuantError'])
    )
    import upsetplot  # imported here: it pulls in matplotlib.pyplot

    ups = upsetplot.from_contents(
        filtered.groupby('search_engine')['precursor ion'].apply(set).to_dict()
    )
    return ups, filtered


def _build_selection_matrix(n_tools):
    """Tuples for unique-per-tool (n rows) + all-tools (1 row)."""
    matrix = [[False] * n_tools for _ in range(n_tools)]
    for i in range(n_tools):
        matrix[i][i] = True
    matrix.append([True] * n_tools)
    return [tuple(row) for row in matrix]


def get_boxplot_data(ups, intermediate, column='AbsQuantError'):
    """Per-tool quantification error for unique-to-tool vs all-tools precursor subsets.

    Returns a long-format DataFrame with columns ``column``, ``tool``,
    ``Quantified per software``.
    """
    tools    = list(ups.index.names)
    prec_map = ups.reset_index().groupby(tools)['id'].apply(set)
    sel      = _build_selection_matrix(len(tools))
    rows = []
    for i, tool in enumerate(tools):
        for label, idx in [('unique to tool', sel[i]), ('all tools', sel[-1])]:
            try:
                precursors = prec_map.loc[idx]
            except KeyError:
                precursors = set()
            vals = intermediate.loc[
                (intermediate['search_engine'] == tool) &
                (intermediate['precursor ion'].isin(precursors)),
                column
            ].tolist()
            rows.append({'values': vals, 'tool': tool, 'subset': label})
    return (
        pd.concat([pd.DataFrame(r) for r in rows], ignore_index=True)
        .rename(columns={'values': column, 'subset': 'Quantified per software'})
    )
