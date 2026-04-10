import os
from matplotlib import pyplot as plt
import warnings


def _save(fig, path: str) -> None:
    base, _ = os.path.splitext(os.path.basename(path))
    root    = os.path.dirname(path)

    for fmt in ('png', 'svg', 'pdf'):
        out = os.path.join(root, fmt, f'{base}.{fmt}')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, bbox_inches='tight', transparent=True, format=fmt)


def _save_legend(ax, save_dir: str, filename: str = 'legend.png') -> None:

    fig, ax_leg = plt.subplots(figsize=(3, len(labels) * 0.35 + 0.3))
    ax_leg.legend(handles, labels, title='Lab', loc='center', frameon=False)
    ax_leg.axis('off')
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()
    
    _save(fig, os.path.join(save_dir, filename))
    plt.close(fig)


def _save_legend_from_mapping(
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

    _save(fig, os.path.join(save_dir, filename))
    plt.close(fig)