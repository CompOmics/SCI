"""Parse search-engine parameter/settings files into :class:`ProteoBenchParameters`.

One module per engine; each exposes an ``extract_params``-style function
(Spectronaut: ``read_spectronaut_settings``). Copied from ProteoBench.

- :mod:`.diann`, :mod:`.spectronaut`, :mod:`.alphadia`, :mod:`.peaks`,
  :mod:`.fragger` (FragPipe), :mod:`.msaid` (CHIMERYS)
- :mod:`.base` -- the :class:`ProteoBenchParameters` container
"""

from .base import ProteoBenchParameters

__all__ = ["ProteoBenchParameters"]
