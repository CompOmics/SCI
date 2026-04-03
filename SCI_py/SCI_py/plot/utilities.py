import os
from matplotlib import pyplot as plt

def _save(fig, path: str) -> None:
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    fig.savefig(path, bbox_inches='tight', transparent=True)


def _save_legend(ax, save_dir: str, filename: str = 'legend.png') -> None:
    handles, labels = ax.get_legend_handles_labels()
    fig, ax_leg = plt.subplots(figsize=(3, len(labels) * 0.35 + 0.3))
    ax_leg.legend(handles, labels, title='Lab', loc='center', frameon=False)
    ax_leg.axis('off')
    fig.tight_layout()
    _save(fig, os.path.join(save_dir, filename))
    plt.close(fig)