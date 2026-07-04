"""Variant logic for the QC/Quant parsing pipeline.

Each result folder yields one ``main`` variant and, depending on engine and
settings, optionally a ``ticnorm`` (TIC-normalised) and a ``firstpass``
(DIA-NN first-pass) variant.
"""

from ..constants import search_settings_names


def firstpass_comment(settings):
    """Derive the first-pass comment label from the DIA-NN settings base name."""
    base = search_settings_names['DIANN'][settings]
    if '(2nd pass)' in base:
        return base.replace('(2nd pass)', '(1st pass)')
    return base + ' (1st pass)'


def needs_ticnorm(engine, settings):
    """Return True if this (engine, settings) should generate a TIC-norm variant.

    - PEAKS: all settings except those already TIC-normalised (settings 4).
    - Chimerys: all settings.
    - Spectronaut: only the ``method-evaluation`` setting (settings 5).
    - All other engines: no TIC-norm variant.
    """
    comment = search_settings_names[engine][settings]
    if engine == 'PEAKS':
        return 'TIC-norm' not in comment  # settings 4 already has TIC-norm in the name
    if engine == 'Chimerys':
        return True
    if engine == 'Spectronaut':
        return comment == 'The real method-evaluation'  # only settings 5
    return False


def has_firstpass(engine, settings):
    """Return True if this entry has a DIA-NN first-pass report."""
    return engine == 'DIANN' and settings in ['4', '5']


def apply_tic_norm(table):
    """Return a copy of a parsed table with TIC-normalised quantities.

    The quantity column normalised depends on ``input_format``:
    PEAKS ``* Normalized Area``, MSAID ``QUANTIFICATION`` (per ``SAMPLE_NAME``),
    Spectronaut ``FG.Quantity`` (per ``R.FileName``). The comment gains a
    ``TIC-norm`` suffix. Other formats are returned with only the comment change.
    """
    report = table['report'].copy()
    fmt = table['input_format']

    if fmt == 'PEAKS':
        norm_cols = [c for c in report.columns if c.endswith(' Normalized Area')]
        report[norm_cols] = report[norm_cols] / report[norm_cols].sum()

    elif fmt == 'MSAID':
        tic = report.groupby('SAMPLE_NAME')['QUANTIFICATION'].sum()
        report['QUANTIFICATION'] = report['QUANTIFICATION'] / report['SAMPLE_NAME'].map(tic)

    elif fmt == 'Spectronaut':
        tic = report.groupby('R.FileName')['FG.Quantity'].sum()
        report['FG.Quantity'] = report['FG.Quantity'] / report['R.FileName'].map(tic)

    return {**table, 'report': report, 'comment': table['comment'] + ' TIC-norm'}
