import os
from matplotlib import pyplot as plt

def _save(fig, path: str) -> None:
    base, _ = os.path.splitext(os.path.basename(path))
    root    = os.path.dirname(path)

    for fmt in ('png', 'svg', 'pdf'):
        out = os.path.join(root, fmt, f'{base}.{fmt}')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, bbox_inches='tight', transparent=True, format=fmt)


def _save_legend(ax, save_dir: str, filename: str = 'legend.png') -> None:
    handles, labels = ax.get_legend_handles_labels()
    fig, ax_leg = plt.subplots(figsize=(3, len(labels) * 0.35 + 0.3))
    ax_leg.legend(handles, labels, title='Lab', loc='center', frameon=False)
    ax_leg.axis('off')
    fig.tight_layout()
    _save(fig, os.path.join(save_dir, filename))
    plt.close(fig)