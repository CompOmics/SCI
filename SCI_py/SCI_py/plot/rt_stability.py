from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

from ..constants.plotting import CONDITION_COLORS, LAB_COLORS
from .utilities import _save


BUILTIN_PALETTES = {
    'source name': CONDITION_COLORS,
    'plate_idx':   LAB_COLORS,
}


def plot_rt_stability(
    rt_props,
    color_by:        str          = 'source name',
    show_individual: bool         = True,
    show_mean:       bool         = True,
    title:           str          = '',
    cmap:            str          = 'viridis',
    palette:         dict | None  = None,
    save_dir:        str | None   = None,
    filename:        str          = 'rt_stability.png',
    show:            bool         = True,
):
    sample_meta = (
        rt_props[['sample_id', color_by]]
        .drop_duplicates()
        .set_index('sample_id')[color_by]
    )

    continuous = (
        color_by not in BUILTIN_PALETTES
        and pd.api.types.is_numeric_dtype(sample_meta)
    )

    if continuous:
        vmin, vmax   = sample_meta.min(), sample_meta.max()
        norm         = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap_obj     = plt.get_cmap(cmap)
        sample_color = {s: cmap_obj(norm(v)) for s, v in sample_meta.items()}
    else:
        _palette = (
            BUILTIN_PALETTES.get(color_by)
            or palette
            or dict(zip(
                sample_meta.unique(),
                plt.rcParams['axes.prop_cycle'].by_key()['color'],
            ))
        )
        sample_color = {s: _palette.get(v, 'black') for s, v in sample_meta.items()}
        groups       = sample_meta.unique()

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)

    if show_individual:
        for sample in sorted(rt_props['sample_id'].unique()):
            sub = rt_props[rt_props['sample_id'] == sample].sort_values('bin_mid')
            ax.plot(sub['bin_mid'], sub['proportion'],
                    color=sample_color.get(sample, 'black'), linewidth=0.8, alpha=0.2)

    if continuous:
        if show_mean:
            for val, grp in rt_props.groupby(color_by):
                mean_curve = (
                    grp.groupby('bin_mid', observed=True)['proportion']
                    .mean().reset_index().sort_values('bin_mid')
                )
                ax.plot(mean_curve['bin_mid'], mean_curve['proportion'],
                        color=cmap_obj(norm(val)), linewidth=2.5)

        sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.8, label=color_by)

    else:
        if show_mean:
            for group in groups:
                mean_curve = (
                    rt_props[rt_props[color_by] == group]
                    .groupby('bin_mid', observed=True)['proportion']
                    .mean().reset_index().sort_values('bin_mid')
                )
                ax.plot(mean_curve['bin_mid'], mean_curve['proportion'],
                        color=_palette.get(group, 'black'), linewidth=2.5)

    ax.set(xlabel='Retention time (bin midpoints)', ylabel='Proportion',
           title=f'{title} RT stability'.strip())
    ax.spines[['top', 'right']].set_visible(False)

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')
        if continuous:
            _save_rt_colorbar(cmap_obj, norm, color_by, save_dir, filename, show)
        else:
            _save_rt_legend(groups, _palette, color_by, save_dir, filename, show)

    plt.show() if show else plt.close(fig)


def _save_rt_colorbar(cmap_obj, norm, label, save_dir, filename, show):
    fig, ax = plt.subplots(figsize=(1.2, 4))
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm),
                 cax=ax, orientation='vertical', label=label)
    _save(fig, f'{save_dir}/{filename.replace(".png", "_legend.png")}')
    plt.show() if show else plt.close(fig)


def _save_rt_legend(groups, palette, color_by, save_dir, filename, show):
    handles = [mlines.Line2D([], [], color=palette.get(g, 'black'), linewidth=2, label=g)
               for g in groups]
    fig, ax = plt.subplots(figsize=(4, len(groups) * 0.4 + 0.5))
    ax.legend(handles=handles, title=color_by, loc='center', frameon=False, fontsize=9)
    ax.axis('off')
    _save(fig, f'{save_dir}/{filename.replace(".png", "_legend.png")}')
    plt.show() if show else plt.close(fig)