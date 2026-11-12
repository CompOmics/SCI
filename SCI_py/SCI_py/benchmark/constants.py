"""Constants for the QC/Quant search-engine benchmark.

``SEARCH_SETTINGS_NAMES`` maps, per search engine, the settings-folder number
(the ``<settings>`` level in ``<lab>/<engine>/<settings>``) to a human-readable
settings label. The label is used as the ``comment`` of each benchmark entry.
"""

SEARCH_SETTINGS_NAMES = {
    'DIANN': {
        '1': 'Library-free',
        '2': 'MBR',
        '3': 'Library-free (recommended settings)',
        '4': 'Default (2nd pass)',
        '5' : 'DIANN 2.6'
    },
    'Chimerys': {
        '1': '5.0.0 MS2(Area)',
        '2': '5.0.0 MS1(Area)',
        '3': '5.0.0 MS1(Area)-MBR',
        '4': '4.4.1 MS2(Area)'
    },
    'AlphaDIA': {
        '1': 'Default',
    },
    'Spectronaut': {
        '1': 'DirectDIA+',
        '2': 'DirectDIA+ method-evaluation',
        '3': 'Library (Astral + timsTOF)',
        '4': 'Library',
        '5': 'The real method-evaluation'
    },
    'PEAKS':{
        '1': 'Recommended',
        '2': 'MBR',
        '3': 'Optimized isolation settings / MaxLFQ',
        '4': 'MaxLFQ / MBR / TIC-norm'
        },
    'FragPipe': {
        '1': 'standard',
        '2': 'DIAUmpire'
    }
}
