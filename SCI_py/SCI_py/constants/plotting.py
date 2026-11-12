"""Colours and display orders used by the plotting functions.

- ``CONDITION_COLORS`` / ``LAB_COLORS`` / ``SOFTWARE_COLORS`` -- hex colour per
  condition (``source name``), anonymised lab and search engine
- ``CONDITION_ORDER`` / ``LAB_ORDER`` -- order in which conditions and labs are shown
"""
from .labs import CONDITIONS, LABS

# Colors
CONDITION_COLORS = {
    'K562':                 '#9B72CF',   # muted purple
    'HEK293T':              '#F4A7B9',   # soft rose pink
    'QC - Hela':            '#F9D776',   # warm muted yellow
    'Site specific QC':     '#2E8B57',   # sage green
    'Cellenone blank':      '#E8735A',   # terracotta red
    'Site specific blanks': '#6DBF8A',   # sea green   # dusty rose / soft coral
    'Processing control':   '#5B8DB8',   # steel blue
    'Sorting control':      '#A8CCE8',   # powder blue
    'QC - HYE A':           "#D87A21",
    'QC - HYE B':           "#AF6D09"
}
LAB_COLORS = {
    'Astral_1':  '#6aaed6',
    'Astral_2':  '#f0b84a',
    'Astral_3':  '#5bbfaa',
    'Astral_4':  '#f0956e',
    'Astral_5':  '#dda0cc',
    'Astral_6':  '#d4b48c',
    'timsTOF_1': '#fcd4ee',
    'timsTOF_2': '#b8b8b8',
    'timsTOF_3': '#f0e87a'
}

# Orderings
# Same values as SCI_py.constants.labs, kept as separate lists.
CONDITION_ORDER = list(CONDITIONS)
LAB_ORDER = list(LABS)

SOFTWARE_COLORS = {
    "MaxQuant": "#8bc6fd",
    "AlphaPept": "#17212b",
    "ProlineStudio": "#8b26ff",
    "MSAngel": "#C0FA7D",
    "FragPipe": "#F89008",
    "i2MassChroQ": "#108E2E",
    "Sage": "#E43924",
    "WOMBAT": "#663200",
    "DIA-NN": "#d42f2f",
    "AlphaDIA": "#1D2732",
    "Custom": "#000000",
    "Spectronaut": "#007548",
    "FragPipe (DIA-NN quant)": "#F89008",
    "MSAID": "#bfef45",
    "MetaMorpheus": "#637C7A",
    "Proteome Discoverer": "#911eb4",
    "PEAKS": "#f032e6",
    "quantms": "#f5e830",
}