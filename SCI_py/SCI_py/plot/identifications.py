from __future__ import annotations
import os
import matplotlib.pyplot as plt
import seaborn as sns
from ..constants.plotting import CONDITION_ORDER, LAB_ORDER, LAB_COLORS, CONDITION_COLORS
from .utilities import _save, _save_legend, _save_legend_from_mapping

import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd
import warnings

def _metric_label(count_col: str) -> str:
    _METRIC_LABELS = {
        'n_proteins'  : 'Number of proteins',
        'n_precursors': 'Number of precursors',
        'n_proteins_withCont'  : 'Number of proteins',
        'n_precursors_withCont': 'Number of precursors',
        'pct_contaminant' : '% contaminant identifications',
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
        _save_legend_from_mapping(LAB_COLORS, save_dir, title='Lab')

    if show:
        plt.show()
    else:
        plt.close(fig)
        

def plot_identification_heatmap(
    adata,
    conditions:              list[str]    = CONDITION_ORDER,
    global_protein_order:    bool         = True,
    sort_conditions:         list[str] | None = None,
    max_samples_per_condition: int | None = None,
    conditions_to_sample:    list[str] | None = None,
    save_dir:                str | None   = None,
    show:                    bool         = True,
) -> dict[str, np.ndarray]:
    """
    Plot per-lab binary identification heatmaps.

    Parameters
    ----------
    adata:
        AnnData object with ``plate_idx`` and ``source name`` in ``obs``.
    conditions:
        Conditions to include and their display order.
    global_protein_order:
        If ``True``, proteins are ordered by global detection frequency.
        If ``False``, ordered per lab.
    sort_conditions:
        Conditions used to determine protein ordering.  ``None`` uses all.
    max_samples_per_condition:
        Maximum samples to display per condition (randomly subsampled).
    conditions_to_sample:
        Conditions to apply subsampling to.  ``None`` applies to all.
    save_dir:
        Directory for saved figures.  ``None`` disables saving.
    show:
        Whether to call ``plt.show()``.

    Returns
    -------
    dict mapping lab → binary identification matrix (samples × proteins).
    """
    condition_order = [c for c in CONDITION_ORDER if c in conditions]

    mask           = adata.obs['source name'].isin(condition_order)
    adata_sub      = adata[mask, :]
    X_full         = adata_sub.X.copy()
    obs_full       = adata_sub.obs.copy()
    global_order   = _protein_order(X_full, obs_full, sort_conditions) \
                     if global_protein_order else None

    binary_dfs = {}
    for lab in LAB_ORDER:
        obs_lab = obs_full[obs_full['plate_idx'] == lab]
        if obs_lab.empty:
            continue

        X_binary, obs_sorted = _prepare_lab(
            obs_lab, X_full, adata_sub,
            condition_order, global_order, global_protein_order,
            max_samples_per_condition, conditions_to_sample,
        )

        fig = plt.figure(figsize=(max(8, X_binary.shape[0] * 0.08 + 2), 8))
        gs  = gridspec.GridSpec(2, 1, height_ratios=[0.02, 1], hspace=0.01)
        _draw_heatmap(fig, gs[1], gs[0], X_binary, obs_sorted, condition_order, title=lab)

        if save_dir:
            _save(fig, f'{save_dir}/{lab}_identification_heatmap.png')
        plt.show() if show else plt.close(fig)

        binary_dfs[lab] = X_binary

    if save_dir:
        _save_heatmap_legend(condition_order, save_dir, show)

    return binary_dfs


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _protein_order(
    X:                np.ndarray,
    obs:              "pd.DataFrame",
    sort_conditions:  list[str] | None,
) -> np.ndarray:
    if sort_conditions is not None:
        mask = obs['source name'].isin(sort_conditions).values
        counts = np.sum(~np.isnan(X[mask, :]), axis=0)
    else:
        counts = np.sum(~np.isnan(X), axis=0)
    return np.argsort(counts)[::-1]


def _sample_obs(obs, condition_order, max_n, conditions_to_sample):
    if max_n is None:
        return obs
    to_sample = conditions_to_sample or condition_order
    parts = []
    for cond in condition_order:
        sub = obs[obs['source name'] == cond]
        if cond in to_sample and max_n < len(sub):
            sub = sub.sample(max_n, random_state=42)
        parts.append(sub)
    return pd.concat(parts)


def _prepare_lab(
    obs_lab, X_full, adata_sub,
    condition_order, global_order, global_protein_order,
    max_samples_per_condition, conditions_to_sample,
):
    obs = _sample_obs(obs_lab, condition_order, max_samples_per_condition, conditions_to_sample)
    idx = [adata_sub.obs.index.get_loc(s) for s in obs.index]
    X   = X_full[idx, :]

    valid_samples  = np.sum(~np.isnan(X), axis=1) > 0
    valid_features = np.sum(~np.isnan(X), axis=0) > 0
    X   = X[np.ix_(valid_samples, valid_features)]
    obs = obs.iloc[valid_samples]

    if global_protein_order:
        global_to_local      = np.full(X_full.shape[1], -1)
        global_to_local[valid_features] = np.arange(valid_features.sum())
        protein_order        = global_to_local[global_order[valid_features[global_order]]]
        protein_order        = protein_order[protein_order >= 0]
    else:
        protein_order = _protein_order(X, obs, None)

    obs = obs.copy()
    obs['n_identified']    = np.sum(~np.isnan(X), axis=1)
    obs['condition_order'] = obs['source name'].map({c: i for i, c in enumerate(condition_order)})
    obs_sorted             = obs.sort_values(['condition_order', 'n_identified'],
                                              ascending=[True, False])
    sort_idx  = [list(obs.index).index(s) for s in obs_sorted.index]
    X_binary  = (~np.isnan(X[np.ix_(sort_idx, protein_order)])).astype(float)

    return X_binary, obs_sorted


def _draw_heatmap(fig, gs_hmap, gs_cbar, X_binary, obs_sorted, condition_order, title=''):
    ax_cbar = fig.add_subplot(gs_cbar)
    ax_hmap = fig.add_subplot(gs_hmap)

    cond_idx  = {c: i for i, c in enumerate(condition_order)}
    cbar_data = obs_sorted['source name'].map(cond_idx).values.astype(float)
    cmap_cbar = ListedColormap([CONDITION_COLORS.get(c, 'grey') for c in condition_order])

    ax_cbar.imshow(cbar_data.reshape(1, -1), aspect='auto', interpolation='none',
                   cmap=cmap_cbar, vmin=0, vmax=len(condition_order) - 1)
    ax_cbar.set(xticks=[], yticks=[], title=title)
    ax_cbar.title.set_fontsize(14)
    ax_cbar.title.set_fontweight('bold')
    for spine in ax_cbar.spines.values():
        spine.set_visible(False)

    ax_hmap.imshow(X_binary.T, aspect='auto', interpolation='none',
                   cmap=ListedColormap(['#f0f0f0', '#2c7bb6']), vmin=0, vmax=1)
    ax_hmap.set(xlabel='Samples', ylabel='Proteins', xticks=[], yticks=[])


def _save_heatmap_legend(condition_order, save_dir, show):
    for handles, title, fname in [
        (
            [Patch(facecolor=CONDITION_COLORS.get(c, 'grey'), label=c) for c in condition_order],
            'Condition',
            'legend_heatmap_conditions.png',
        ),
        (
            [Patch(facecolor='#2c7bb6', label='Identified'),
             Patch(facecolor='#f0f0f0', edgecolor='grey', label='Missing')],
            'Identification',
            'legend_heatmap_identification.png',
        ),
    ]:
        fig, ax = plt.subplots(figsize=(3, len(handles) * 0.35 + 0.5))
        ax.legend(handles=handles, loc='center', frameon=False,
                  fontsize=9, title=title, title_fontsize=10)
        ax.axis('off')
        _save(fig, f'{save_dir}/{fname}')
        plt.show() if show else plt.close(fig)