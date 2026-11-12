"""
Coefficient of variation (CV) plots per lab and condition.

Part of the counts and CVs analysis (Notebooks/2_DatasetStatistics/01_Counts_CVs).
CVs are computed in ``SCI_py.metrics.cvs``.
"""
from __future__ import annotations

import os
import matplotlib.pyplot as plt
import seaborn as sns

from ..constants.plotting import CONDITION_ORDER, LAB_ORDER, LAB_COLORS, CONDITION_COLORS
from ..metrics.cvs import filter_common_proteome
from .export import save_figure, save_legend_from_mapping

def plot_cvs_overview(
    cvs,
    conditions:       list[str]    = CONDITION_ORDER,
    min_n:            int          = 10,
    common_proteome:  bool         = False,
    color_by:         str          = 'source name',
    plot_type:        str          = 'violin',
    ylabel:           str          = 'CV (%)',
    split_conditions: bool         = False,
    ylim:             float | None = None,
    save_dir:         str | None   = None,
    filename:         str          = 'cvs_overview.png',
    show:             bool         = True,
):
    """
    Plot coefficient of variation distributions across labs and conditions.

    Parameters
    ----------
    cvs:
        DataFrame with columns ``feature_id``, ``lab``, ``condition``,
        ``cv``, and ``n`` (as returned by ``compute_cv_matrix``).
    conditions:
        Conditions to include. Defaults to ``CONDITION_ORDER``.
    min_n:
        Minimum number of observations required to include a CV value.
    common_proteome:
        If ``True``, restrict to features quantified in all labs per condition.
    color_by:
        ``'source name'`` — colour by condition, group by lab.
        ``'lab'``         — colour by lab, group by condition.
        Ignored when ``split_conditions=True``.
    plot_type:
        ``'violin'`` or ``'box'``.
    ylabel:
        Y-axis label.
    split_conditions:
        If ``True``, one subplot per condition with labs on the x-axis.
    ylim:
        Upper y-axis limit in percent.  Lower limit is always 0.
    save_dir:
        Directory for saved figures.  ``None`` disables saving.
    filename:
        Filename for the main figure.
    show:
        Whether to call ``plt.show()``.
    """
    if color_by not in ('lab', 'source name'):
        raise ValueError("color_by must be 'lab' or 'source name'")
    if plot_type not in ('violin', 'box'):
        raise ValueError("plot_type must be 'violin' or 'box'")

    # --- Filtering ---
    condition_order = [c for c in CONDITION_ORDER if c in conditions]
    cvs = cvs[cvs['condition'].isin(condition_order) & (cvs['n'] >= min_n)].copy()

    if common_proteome:
        cvs = filter_common_proteome(cvs, condition_order)

    # --- Plot dispatch ---
    _plot_fn = _violin if plot_type == 'violin' else _box

    if split_conditions:
        fig, axes = _plot_split(cvs, condition_order, ylabel, ylim,
                                common_proteome, _plot_fn)
    else:
        fig, axes = _plot_combined(cvs, condition_order, color_by,
                                   ylabel, ylim, _plot_fn)

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
        save_legend_from_mapping(LAB_COLORS, save_dir, title='Lab', filename='legend_cvs.png')

    plt.show() if show else plt.close(fig)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _violin(data, x, y, order, palette, hue=None, hue_order=None, ax=None):
    """Draw a seaborn violin plot (quartile lines), coloured by ``x`` when no ``hue`` is given."""
    kwargs = dict(width=0.8, inner='quartile')
    if hue:
        sns.violinplot(data=data, x=x, y=y, hue=hue, hue_order=hue_order,
                       order=order, palette=palette, ax=ax, **kwargs)
    else:
        sns.violinplot(data=data, x=x, y=y, hue=x, hue_order=order,
                       order=order, palette=palette, legend=False,
                       ax=ax, **kwargs)


def _box(data, x, y, order, palette, hue=None, hue_order=None, ax=None):
    """Draw a seaborn box plot (small fliers), coloured by ``x`` when no ``hue`` is given."""
    kwargs = dict(width=0.6, flierprops=dict(marker='o', markersize=2, alpha=0.4))
    if hue:
        sns.boxplot(data=data, x=x, y=y, hue=hue, hue_order=hue_order,
                    order=order, palette=palette, ax=ax, **kwargs)
    else:
        sns.boxplot(data=data, x=x, y=y, hue=x, hue_order=order,
                    order=order, palette=palette, legend=False,
                    ax=ax, **kwargs)


def _plot_split(cvs, condition_order, ylabel, ylim, common_proteome, plot_fn):
    """Draw one CV panel per condition with labs on the x-axis; return ``(fig, axes)``."""
    n     = len(condition_order)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 5),
                              sharey=True, constrained_layout=True)
    axes = [axes] if n == 1 else list(axes)

    for ax, condition in zip(axes, condition_order):
        sub        = cvs[cvs['condition'] == condition]
        n_features = sub['feature_id'].nunique()

        plot_fn(data=sub, x='lab', y='cv', order=LAB_ORDER,
                palette=LAB_COLORS, ax=ax)

        ax.set(title  = condition,
               xlabel = '',
               ylabel = ylabel if ax == axes[0] else '')
        ax.set_xticklabels([])
        ax.set_xlabel(f'n={n_features}', labelpad=6, fontsize=9, color='grey')
        if ylim is not None:
            ax.set_ylim(0, ylim)

    proteome_label = 'common proteome' if common_proteome else 'all features'
    fig.suptitle(f'{ylabel} — {proteome_label}', fontweight='bold')
    return fig, axes


def _plot_combined(cvs, condition_order, color_by, ylabel, ylim, plot_fn):
    """Draw all CVs on one axis, grouped and coloured per ``color_by``; return ``(fig, [ax])``."""
    if color_by == 'source name':
        x, hue, order, hue_order, palette = (
            'lab', 'condition', LAB_ORDER, condition_order, CONDITION_COLORS
        )
    else:
        x, hue, order, hue_order, palette = (
            'condition', 'lab', condition_order, LAB_ORDER, LAB_COLORS
        )

    fig, ax = plt.subplots(1, 1, figsize=(14, 5), constrained_layout=True)
    plot_fn(data=cvs, x=x, y='cv', order=order, palette=palette,
            hue=hue, hue_order=hue_order, ax=ax)

    ax.set(xlabel='', ylabel=ylabel, title=ylabel)
    ax.tick_params(axis='x', rotation=45)
    ax.get_legend().remove()
    if ylim is not None:
        ax.set_ylim(0, ylim)

    return fig, [ax]