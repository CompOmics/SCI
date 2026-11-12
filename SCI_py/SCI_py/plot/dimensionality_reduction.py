"""
PCA and 2D-embedding scatter grids, coloured by ``obs`` columns.

Used to look at sample structure (lab, condition) across datasets. The
embeddings are computed in ``SCI_py.metrics.dimensionality_reduction``.
"""
from __future__ import annotations

import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

from ..constants.plotting import CONDITION_COLORS, CONDITION_ORDER, LAB_COLORS, LAB_ORDER
from .export import save_figure

# Built-in (palette, category order) for obs columns with a project-wide convention.
_BUILTIN_PALETTES = {
    'source name': (CONDITION_COLORS, CONDITION_ORDER),
    'plate_idx':   (LAB_COLORS, LAB_ORDER),
}


def plot_pca_comparison(
    embeddings:        dict[str, dict],
    obs:               pd.DataFrame,
    color_by:          list[str],
    components:        tuple[int, int]   = (0, 1),
    palettes:          dict[str, dict]   | None = None,
    orders:            dict[str, list]   | None = None,
    point_size:        float             = 10,
    alpha:             float             = 0.5,
    panel_size:        tuple[float, float] = (5.2, 5),
    suptitle:          str | None        = None,
    save_dir:          str | None        = None,
    filename:          str               = 'pca_comparison.png',
    show:              bool              = True,
):
    """
    Grid of PCA scatter panels: one row per `color_by` variable, one column per
    embedding — typically used to compare a raw vs. batch-corrected embedding.

    Parameters
    ----------
    embeddings:
        Dict as returned by
        `SCI_py.metrics.dimensionality_reduction.compute_pca_embeddings`, i.e.
        name -> {'pcs': ..., 'explained_variance_ratio': ...}. Column order in the
        figure follows the order of this dict's keys.
    obs:
        Sample metadata, aligned row-for-row (same order, same length) with the
        samples used to fit `embeddings`.
    color_by:
        obs column name(s) to colour points by, one row per entry. Columns
        'source name' / 'plate_idx' use the project's built-in palette; anything
        else needs a matching entry in `palettes` (and optionally `orders`), or
        falls back to the default color cycle.
    components:
        Which zero-indexed principal components to plot as (x, y).
    palettes / orders:
        Optional dict mapping a `color_by` entry (without a built-in palette) to a
        {category: color} dict / a category ordering.
    suptitle:
        Optional figure-level title, e.g. to label which sample set is shown.
    save_dir / filename / show:
        Standard figure-saving controls (see `SCI_py.plot.export.save_figure`).
    """
    palettes = palettes or {}
    orders   = orders or {}
    pc_x, pc_y = components

    embedding_names = list(embeddings.keys())
    fig, axes = plt.subplots(
        len(color_by), len(embedding_names),
        figsize=(panel_size[0] * len(embedding_names), panel_size[1] * len(color_by)),
        squeeze=False,
    )

    for row, color_col in enumerate(color_by):
        labels = obs[color_col].to_numpy()

        if color_col in _BUILTIN_PALETTES:
            color_map, order = _BUILTIN_PALETTES[color_col]
        elif color_col in palettes:
            color_map = palettes[color_col]
            order     = orders.get(color_col, list(color_map.keys()))
        else:
            uniq      = pd.unique(labels)
            cycle     = plt.rcParams['axes.prop_cycle'].by_key()['color']
            color_map = dict(zip(uniq, cycle * (len(uniq) // len(cycle) + 1)))
            order     = list(uniq)

        categories = [c for c in order if c in set(labels)]

        for col, name in enumerate(embedding_names):
            ax  = axes[row, col]
            pcs = embeddings[name]['pcs']
            for cat in categories:
                mask = labels == cat
                ax.scatter(pcs[mask, pc_x], pcs[mask, pc_y], s=point_size, alpha=alpha,
                           color=color_map[cat], linewidths=0, label=cat)

            var = embeddings[name]['explained_variance_ratio']
            ax.set_xlabel(f"PC{pc_x + 1} ({var[pc_x]:.1%})")
            ax.set_ylabel(f"PC{pc_y + 1} ({var[pc_y]:.1%})")
            ax.set_title(f"{name}  —  colour: {color_col}")

        handles = [Line2D([0], [0], marker='o', linestyle='', color=color_map[c], label=c)
                   for c in categories]
        axes[row, -1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                             frameon=False, fontsize=8, title=color_col)

    fig.tight_layout()

    if suptitle:
        fig.suptitle(suptitle, fontsize=13, fontweight='bold')
        fig.subplots_adjust(top=0.90)

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))

    plt.show() if show else plt.close(fig)

    return fig, axes


_METHOD_LABELS = {'pca': 'PCA', 'umap': 'UMAP', 'tsne': 't-SNE'}


def plot_embedding_comparison(
    embeddings:        dict[str, dict],
    obs:               pd.DataFrame,
    color_by:          list[str],
    methods:           list[str]         = ('pca', 'umap', 'tsne'),
    palettes:          dict[str, dict]   | None = None,
    orders:            dict[str, list]   | None = None,
    point_size:        float             = 10,
    alpha:             float             = 0.5,
    panel_size:        tuple[float, float] = (5.2, 5),
    suptitle:          str | None        = None,
    save_dir:          str | None        = None,
    filename:          str               = 'embedding_comparison.png',
    show:              bool              = True,
):
    """
    Grid of 2D embedding scatter panels: one row per `color_by` variable, one column
    per (matrix name, method) combination — e.g. Raw-PCA, Raw-UMAP, Raw-t-SNE,
    Corrected-PCA, Corrected-UMAP, Corrected-t-SNE — for comparing dimensionality
    reduction methods across, typically, a raw vs. batch-corrected embedding.

    Parameters
    ----------
    embeddings:
        Dict as returned by
        `SCI_py.metrics.dimensionality_reduction.compute_embeddings`, i.e.
        name -> method -> {'coords': ..., ...}. Column order follows this dict's
        keys, then `methods`.
    obs:
        Sample metadata, aligned row-for-row (same order, same length) with the
        samples used to fit `embeddings`.
    color_by:
        obs column name(s) to colour points by, one row per entry. Columns
        'source name' / 'plate_idx' use the project's built-in palette; anything
        else needs a matching entry in `palettes` (and optionally `orders`), or
        falls back to the default color cycle.
    methods:
        Which methods to draw as columns, for matrices that have them (subset of
        ('pca', 'umap', 'tsne')).
    palettes / orders:
        Optional dict mapping a `color_by` entry (without a built-in palette) to a
        {category: color} dict / a category ordering.
    suptitle:
        Optional figure-level title, e.g. to label which sample set is shown.
    save_dir / filename / show:
        Standard figure-saving controls (see `SCI_py.plot.export.save_figure`).
    """
    palettes = palettes or {}
    orders   = orders or {}

    columns = [(name, method) for name in embeddings for method in methods if method in embeddings[name]]
    if not columns:
        raise ValueError(f"None of methods={methods} found in embeddings.")

    fig, axes = plt.subplots(
        len(color_by), len(columns),
        figsize=(panel_size[0] * len(columns), panel_size[1] * len(color_by)),
        squeeze=False,
    )

    for row, color_col in enumerate(color_by):
        labels = obs[color_col].to_numpy()

        if color_col in _BUILTIN_PALETTES:
            color_map, order = _BUILTIN_PALETTES[color_col]
        elif color_col in palettes:
            color_map = palettes[color_col]
            order     = orders.get(color_col, list(color_map.keys()))
        else:
            uniq      = pd.unique(labels)
            cycle     = plt.rcParams['axes.prop_cycle'].by_key()['color']
            color_map = dict(zip(uniq, cycle * (len(uniq) // len(cycle) + 1)))
            order     = list(uniq)

        categories = [c for c in order if c in set(labels)]

        for col, (name, method) in enumerate(columns):
            ax     = axes[row, col]
            entry  = embeddings[name][method]
            coords = entry['coords']
            for cat in categories:
                mask = labels == cat
                ax.scatter(coords[mask, 0], coords[mask, 1], s=point_size, alpha=alpha,
                           color=color_map[cat], linewidths=0, label=cat)

            label = _METHOD_LABELS.get(method, method)
            if method == 'pca':
                var = entry['explained_variance_ratio']
                ax.set_xlabel(f"PC1 ({var[0]:.1%})")
                ax.set_ylabel(f"PC2 ({var[1]:.1%})")
            else:
                ax.set_xlabel(f"{label}1")
                ax.set_ylabel(f"{label}2")
            ax.set_title(f"{name}  —  {label}  —  colour: {color_col}")

        handles = [Line2D([0], [0], marker='o', linestyle='', color=color_map[c], label=c)
                   for c in categories]
        axes[row, -1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                             frameon=False, fontsize=8, title=color_col)

    fig.tight_layout()

    if suptitle:
        fig.suptitle(suptitle, fontsize=13, fontweight='bold')
        fig.subplots_adjust(top=0.90)

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))

    plt.show() if show else plt.close(fig)

    return fig, axes
