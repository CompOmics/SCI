from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from ..constants.plotting import CONDITION_COLORS, LAB_COLORS, LAB_ORDER
from .utilities import _save

def plot_rt_stability(
    rt_props,
    color_by:        str             = 'source name',
    show_individual: bool            = True,
    title:           str             = '',
    save_dir:        str | None      = None,
    filename:        str             = 'rt_stability.png',
    show:            bool            = True,
):
    """
    Plot retention time stability as proportion curves per sample and group.

    Parameters
    ----------
    rt_props:
        Output of ``compute_rt_proportions``.
    color_by:
        ``'source name'`` — colour by condition.
        ``'plate_idx'``   — colour by lab.
    show_individual:
        If ``True``, plot individual sample curves behind group means.
    title:
        Figure title prefix.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    if color_by not in ('source name', 'plate_idx'):
        raise ValueError("color_by must be 'source name' or 'plate_idx'")

    palette  = CONDITION_COLORS if color_by == 'source name' else LAB_COLORS
    groups   = rt_props[color_by].unique()
    sample_color = (
        rt_props[['sample_id', color_by]]
        .drop_duplicates()
        .set_index('sample_id')[color_by]
        .to_dict()
    )

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)

    if show_individual:
        for sample in sorted(rt_props['sample_id'].unique()):
            sub   = rt_props[rt_props['sample_id'] == sample].sort_values('bin_mid')
            color = palette.get(sample_color.get(sample, ''), 'black')
            ax.plot(sub['bin_mid'], sub['proportion'],
                    color=color, linewidth=0.8, alpha=0.2)

    for group in groups:
        mean_curve = (
            rt_props[rt_props[color_by] == group]
            .groupby('bin_mid', observed=True)['proportion'].mean()
            .reset_index()
            .sort_values('bin_mid')
        )
        ax.plot(mean_curve['bin_mid'], mean_curve['proportion'],
                color=palette.get(group, 'black'), linewidth=2.5)

    ax.set(xlabel='Retention time (bin midpoints)', ylabel='Proportion',
           title=f'{title} RT stability'.strip())
    ax.spines[['top', 'right']].set_visible(False)

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')
        _save_rt_legend(groups, palette, color_by, save_dir, filename, show)

    plt.show() if show else plt.close(fig)


def _save_rt_legend(groups, palette, color_by, save_dir, filename, show):
    handles = [
        mlines.Line2D([], [], color=palette.get(g, 'black'), linewidth=2, label=g)
        for g in groups
    ]
    legend_filename = filename.replace('.png', '_legend.png')
    fig, ax = plt.subplots(figsize=(4, len(groups) * 0.4 + 0.5))
    ax.legend(handles=handles, title=color_by, loc='center', frameon=False, fontsize=9)
    ax.axis('off')
    _save(fig, f'{save_dir}/{legend_filename}')
    plt.show() if show else plt.close(fig)