"""Entry enumeration for the QC/Quant parsing pipeline."""

import os
from glob import glob

from .constants import SEARCH_SETTINGS_NAMES
from .variants import needs_ticnorm, has_firstpass


def enumerate_entries(root_qcquant):
    """Walk the QC/Quant directory and list all entries to process.

    The directory is expected to follow ``<root>/<lab>/<engine>/<settings>``.
    Folders whose engine or settings are not registered in
    ``SEARCH_SETTINGS_NAMES`` are skipped.

    Parameters
    ----------
    root_qcquant : str
        Root directory containing the per-lab result folders.

    Returns
    -------
    list of tuple
        ``(result_folder, lab, engine, settings, variant)`` tuples, where
        ``variant`` is ``'main'`` plus, where applicable, ``'ticnorm'`` and
        ``'firstpass'``.
    """
    entries = []

    for result_folder in sorted(glob(os.path.join(root_qcquant, '*', '*', '*'))):
        parts = result_folder.split('/')
        lab, engine, settings = parts[-3], parts[-2], parts[-1]

        if engine not in SEARCH_SETTINGS_NAMES:
            continue
        if settings not in SEARCH_SETTINGS_NAMES[engine]:
            continue

        entries.append((result_folder, lab, engine, settings, 'main'))

        if needs_ticnorm(engine, settings):
            entries.append((result_folder, lab, engine, settings, 'ticnorm'))

        if has_firstpass(engine, settings):
            entries.append((result_folder, lab, engine, settings, 'firstpass'))

    return entries
