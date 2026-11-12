"""
Figure export helpers shared by all ``SCI_py.plot`` modules.

``save_figure`` writes every figure as PNG, SVG and PDF; the legend helpers
save stand-alone legend figures. Importing this module also sets matplotlib
font types so text stays editable in PDF/SVG output.
"""
import os
import matplotlib
from matplotlib import pyplot as plt
import warnings

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'  # keep SVG text as editable <text> elements


def save_figure(fig, path: str) -> None:
    """
    Save a figure as PNG, SVG and PDF.

    Each format goes into its own subfolder next to ``path``:
    ``<dir>/png/<name>.png``, ``<dir>/svg/<name>.svg`` and
    ``<dir>/pdf/<name>.pdf``. The extension of ``path`` is ignored.
    Figures are saved with a tight bounding box and a transparent background.

    Parameters
    ----------
    fig:
        Matplotlib figure to save.
    path:
        Target path; only its directory and base name are used.
    """
    base, _ = os.path.splitext(os.path.basename(path))
    root    = os.path.dirname(path)

    for fmt in ('png', 'svg', 'pdf'):
        out = os.path.join(root, fmt, f'{base}.{fmt}')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, bbox_inches='tight', transparent=True, format=fmt)


def save_legend(ax, save_dir: str, filename: str = 'legend.png') -> None:
    """
    Save the legend of an Axes as a stand-alone figure.

    Handles and labels are taken from ``ax``; the legend title is ``'Lab'``.

    Parameters
    ----------
    ax:
        Axes whose legend entries are copied.
    save_dir:
        Directory to save the legend figure in (via ``save_figure``).
    filename:
        Output filename.
    """
    handles, labels = ax.get_legend_handles_labels()
    fig, ax_leg = plt.subplots(figsize=(3, len(labels) * 0.35 + 0.3))
    ax_leg.legend(handles, labels, title='Lab', loc='center', frameon=False)
    ax_leg.axis('off')
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()
    
    save_figure(fig, os.path.join(save_dir, filename))
    plt.close(fig)


def save_legend_from_mapping(
    color_mapping: dict,
    save_dir:      str,
    title:         str           = 'Legend',
    filename:      str           = 'legend.png',
) -> None:
    """
    Create and save a standalone legend from a color mapping.

    Parameters
    ----------
    color_mapping:
        Dict mapping label → color (e.g. LAB_COLORS).
    save_dir:
        Directory to save the legend figure.
    title:
        Legend title.
    filename:
        Output filename.
    """
    from matplotlib.patches import Patch

    handles = [Patch(facecolor=color, label=label)
               for label, color in color_mapping.items()]

    fig, ax = plt.subplots(figsize=(3, len(handles) * 0.35 + 0.3))
    ax.legend(handles, color_mapping.keys(), title=title, loc='center', frameon=False)
    ax.axis('off')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()

    save_figure(fig, os.path.join(save_dir, filename))
    plt.close(fig)