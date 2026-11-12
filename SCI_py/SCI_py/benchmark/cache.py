"""Per-entry pickle cache for the QC/Quant parsing pipeline.

Each parsed entry is stored as
``{cache_dir}/{lab}__{engine}__{settings}__{variant}.pkl``. The cached value is
a dict with keys ``datapoint`` (DataFrame) and ``intermediates`` (list of
DataFrames). Caching lets the pipeline add new entries without re-parsing the
whole dataset; delete a ``.pkl`` file to force that entry to be re-parsed.
"""

import os
import pickle


class EntryCache:
    """Read/write access to the per-entry parse cache under ``cache_dir``."""

    def __init__(self, cache_dir):
        """Store ``cache_dir`` and create it if it does not exist."""
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)

    def path(self, lab, engine, settings, variant):
        """Return the cache file path for one entry."""
        key = f'{lab}__{engine}__{settings}__{variant}'
        return os.path.join(self.cache_dir, f'{key}.pkl')

    def is_cached(self, lab, engine, settings, variant):
        """Return True if this entry is already cached."""
        return os.path.exists(self.path(lab, engine, settings, variant))

    def save(self, lab, engine, settings, variant, datapoint, intermediates):
        """Write one parsed entry to the cache."""
        with open(self.path(lab, engine, settings, variant), 'wb') as f:
            pickle.dump({'datapoint': datapoint, 'intermediates': intermediates}, f)

    def load(self, lab, engine, settings, variant):
        """Load one cached entry, or return ``None`` if it is not cached."""
        path = self.path(lab, engine, settings, variant)
        if not os.path.exists(path):
            return None
        with open(path, 'rb') as f:
            return pickle.load(f)
