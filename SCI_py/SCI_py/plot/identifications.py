from __future__ import annotations
import os
import matplotlib.pyplot as plt
import seaborn as sns
from .constants import CONDITION_ORDER, LAB_ORDER, LAB_COLORS
from .utilities import _save, _save_legend

def _metric_label(count_col: str) -> str:
    _METRIC_LABELS = {
        'n_proteins'  : 'Number of proteins',
        'n_precursors': 'Number of precursors',
    }
    return _METRIC_LABELS.get(count_col, count_col.replace('_', ' ').title())

def plot_counts_overview(
    adata,
    count_col:        str,
    conditions:       list[str]     = CONDITION_ORDER,
    split_conditions: bool          = False,
    show_xlabels:     bool          = False,
    save_dir:         str | None    = None,
    filename:         str           = 'counts_overview.png',
    title:            str | None    = None,
    show:             bool          = True,
):
    """
    Plot per-cell identification counts grouped by lab.

    Parameters
    ----------
    adata:
        AnnData object with obs columns ``source name``, ``plate_idx``,
        and ``count_col``.
    count_col:
        Column in ``adata.obs`` to plot on the y-axis
        (e.g. ``'n_proteins'``, ``'n_precursors'``).
    conditions:
        Source names to include.  Defaults to ``CONDITION_ORDER``.
    split_conditions:
        If ``True``, produce one subplot per condition (source name).
        If ``False`` (default), all conditions share a single axis with
        ``source name`` on the x-axis.
    show_xlabels:
        If ``True``, shows the xlabels when split_conditions is also set to True.
        If ``False`` (default), xlables are stripped. Colors are still indicative.
    save_dir:
        Directory to save figures.  If ``None``, figures are not saved.
    filename:
        Filename for the main figure (used only when ``save_dir`` is set).
    title:
        Custom figure title.  Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    ylabel          = _metric_label(count_col)
    condition_order = [c for c in CONDITION_ORDER if c in conditions]
    data            = adata.obs[adata.obs['source name'].isin(conditions)]

    _boxplot_kwargs = dict(
        hue       = 'plate_idx',
        hue_order = LAB_ORDER,
        palette   = LAB_COLORS,
        width     = 0.6,
        flierprops= dict(marker='o', markersize=2, alpha=0.4),
    )

    if split_conditions:
        n     = len(condition_order)
        fig, axes = plt.subplots(1, n, figsize=(4 * n, 5), sharey=False)
        axes  = [axes] if n == 1 else list(axes)

        for ax, condition in zip(axes, condition_order):
            sns.boxplot(
                data  = data[data['source name'] == condition],
                x     = 'plate_idx',
                y     = count_col,
                order = LAB_ORDER,
                ax    = ax,
                **_boxplot_kwargs,
            )
            
            if show_xlabels:
                ax.set(title=condition, xlabel='Lab', ylabel=ylabel if ax == axes[0] else '')
                ax.tick_params(axis='x', rotation=90)
            else:
                ax.set(title=condition, xlabel='', ylabel=ylabel if ax == axes[0] else '')
                ax.set_xticklabels([])

        fig.suptitle(title or f'{ylabel} per lab', fontweight='bold')

    else:
        fig, ax = plt.subplots(figsize=(max(14, 2 * len(condition_order)), 6))
        sns.boxplot(
            data  = data,
            x     = 'source name',
            y     = count_col,
            order = condition_order,
            ax    = ax,
            **_boxplot_kwargs,
        )
        ax.set(
            xlabel = 'Condition',
            ylabel = ylabel,
            title  = title or f'{ylabel} per condition',
        )
        ax.tick_params(axis='x', rotation=90)
        ax.get_legend().remove()
        axes = [ax]

    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
        _save_legend(axes[0], save_dir)

    if show:
        plt.show()
    else:
        plt.close(fig)