"""
Retention-time distribution (RT stability) curves per sample.

Used in Notebooks/2_DatasetStatistics/02_QC_Metrics and 03_Drift_over_time.
The per-sample RT proportions are computed in ``SCI_py.metrics.rt_stability``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

from ..constants.plotting import CONDITION_COLORS, LAB_COLORS
from .export import save_figure


_BUILTIN_PALETTES = {
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
    """
    Plot per-sample retention-time distributions, with group mean curves.

    Each sample is a thin line of the proportion of identifications per RT
    bin. If ``color_by`` is a numeric column (and not ``source name`` or
    ``plate_idx``), samples are coloured with a colormap and a colorbar is
    added; otherwise they are coloured by category.

    Parameters
    ----------
    rt_props:
        Long-format DataFrame with ``sample_id``, ``bin_mid``,
        ``proportion`` and the ``color_by`` column (as returned by
        ``compute_rt_proportions``).
    color_by:
        Column used to colour samples and group the mean curves.
    show_individual:
        Whether to draw the per-sample curves.
    show_mean:
        Whether to draw the mean curve per group (or per value for a
        numeric ``color_by``).
    title:
        Prefix for the plot title (``'<title> RT stability'``).
    cmap:
        Colormap name for a numeric ``color_by``.
    palette:
        Dict mapping category to colour. Only used when ``color_by`` is not
        ``source name`` or ``plate_idx`` (those have built-in palettes).
        Falls back to the matplotlib colour cycle.
    save_dir:
        Directory for saved figures. ``None`` disables saving. When saving,
        a separate legend (or colorbar) figure is saved as
        ``<filename>_legend``.
    filename:
        Filename for the main figure.
    show:
        Whether to call ``plt.show()``.
    """
    sample_meta = (
        rt_props[['sample_id', color_by]]
        .drop_duplicates()
        .set_index('sample_id')[color_by]
    )

    continuous = (
        color_by not in _BUILTIN_PALETTES
        and pd.api.types.is_numeric_dtype(sample_meta)
    )

    if continuous:
        vmin, vmax   = sample_meta.min(), sample_meta.max()
        norm         = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap_obj     = plt.get_cmap(cmap)
        sample_color = {s: cmap_obj(norm(v)) for s, v in sample_meta.items()}
    else:
        _palette = (
            _BUILTIN_PALETTES.get(color_by)
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
        save_figure(fig, f'{save_dir}/{filename}')
        if continuous:
            _save_rt_colorbar(cmap_obj, norm, color_by, save_dir, filename, show)
        else:
            _save_rt_legend(groups, _palette, color_by, save_dir, filename, show)

    plt.show() if show else plt.close(fig)


def _save_rt_colorbar(cmap_obj, norm, label, save_dir, filename, show):
    """Save a stand-alone vertical colorbar for a numeric ``color_by``."""
    fig, ax = plt.subplots(figsize=(1.2, 4))
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm),
                 cax=ax, orientation='vertical', label=label)
    save_figure(fig, f'{save_dir}/{filename.replace(".png", "_legend.png")}')
    plt.show() if show else plt.close(fig)


def _save_rt_legend(groups, palette, color_by, save_dir, filename, show):
    """Save a stand-alone line legend for a categorical ``color_by``."""
    handles = [mlines.Line2D([], [], color=palette.get(g, 'black'), linewidth=2, label=g)
               for g in groups]
    fig, ax = plt.subplots(figsize=(4, len(groups) * 0.4 + 0.5))
    ax.legend(handles=handles, title=color_by, loc='center', frameon=False, fontsize=9)
    ax.axis('off')
    save_figure(fig, f'{save_dir}/{filename.replace(".png", "_legend.png")}')
    plt.show() if show else plt.close(fig)