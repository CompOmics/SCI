"""Entry enumeration for the QC/Quant parsing pipeline."""

import os
from glob import glob

from ..constants import search_settings_names
from .variants import needs_ticnorm, has_firstpass


def enumerate_entries(root_qcquant):
    """Walk the QC/Quant directory and list all entries to process.

    The directory is expected to follow ``<root>/<lab>/<engine>/<settings>``.
    Folders whose engine or settings are not registered in
    ``search_settings_names`` are skipped.

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

        if engine not in search_settings_names:
            continue
        if settings not in search_settings_names[engine]:
            continue

        entries.append((result_folder, lab, engine, settings, 'main'))

        if needs_ticnorm(engine, settings):
            entries.append((result_folder, lab, engine, settings, 'ticnorm'))

        if has_firstpass(engine, settings):
            entries.append((result_folder, lab, engine, settings, 'firstpass'))

    return entries
