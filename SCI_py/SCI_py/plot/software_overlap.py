"""
UpSet and error-boxplot panels for precursor overlap between search engines.

Part of the benchmark analysis (Notebooks/1_QCQuant/04_Overlap_between_search_engines).
The input tables are built in ``SCI_py.metrics.software_overlap``.
"""
import os
import warnings
from math import ceil

import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot
import string

from .export import save_figure

# ── Panel figure functions ────────────────────────────────────────────────────

def plot_upsetplot_panel(
    slices,
    labs,
    ordering_tools,
    nrow:      int        = 3,
    ncol:      int        = 3,
    save_dir:  str | None = None,
    filename:  str        = 'upsetplot_panel_3x3',
    show:      bool       = True,
):
    """3×3 UpSet panel — one subfigure per lab.

    Only shows degree-1 (unique), degree-(n-1) (all but one), and degree-n (all) subsets.
    Bar order within each degree follows ``ordering_tools``.
    Colours: steelblue = unique, #888888 = all-but-one, darkorange = all.

    Parameters
    ----------
    slices:
        Dict mapping lab name to ``{'ups': ..., 'boxplot': ...}`` as produced by the
        per-lab computation cell.
    labs:
        Ordered list of lab names to plot (at most nrow * ncol).
    ordering_tools:
        Tool display order for both category rows and bar sequence.
    """
    import matplotlib as mpl

    # Scale up the text drawn inside each UpSet subfigure (category labels,
    # intersection counts, totals) so the panel stays legible in print.
    old_rcparams = mpl.rcParams.copy()
    mpl.rcParams.update({
        'font.size':        15,
        'axes.labelsize':   16,
        'axes.titlesize':   18,
        'xtick.labelsize':  14,
        'ytick.labelsize':  14,
    })

    try:
        fig     = plt.figure(figsize=(24, 18))
        subfigs = fig.subfigures(nrow, ncol, hspace=0.06, wspace=0.04)

        panel_letters = string.ascii_lowercase

        for i, lab in enumerate(labs[:nrow * ncol]):
            subfig = subfigs.flat[i]
            ups    = slices[lab]['ups']
            tools  = list(ups.index.names)
            n_t    = len(tools)

            # Filter to desired degrees
            idx_df       = ups.index.to_frame()
            degrees      = idx_df.sum(axis=1)
            keep_degrees = {1, n_t} if n_t <= 2 else {1, n_t - 1, n_t}
            ups_filt     = ups[degrees.isin(keep_degrees)]

            # Build desired bar order
            tools_in_order = [t for t in ordering_tools if t in tools]
            desired = []
            for tool in tools_in_order:
                desired.append(tuple(t == tool for t in tools))
            if n_t > 2:
                for tool in reversed(tools_in_order):
                    desired.append(tuple(t != tool for t in tools))
            desired.append(tuple(True for _ in tools))
            desired     = [tup for tup in desired if tup in ups_filt.index]
            ups_ordered = ups_filt.loc[desired]

            upset = upsetplot.UpSet(
                ups_ordered,
                sort_by='input',
                sort_categories_by='input',
                show_counts='{:,}',
                totals_plot_elements=3,
                element_size=None,
            )
            for subset in ups_ordered.index:
                degree = sum(subset)
                if degree == 1:
                    color = 'steelblue'
                elif degree == n_t:
                    color = 'darkorange'
                else:
                    color = '#888888'
                upset.style_subsets(
                    present=[name for name, val in zip(tools, subset) if val],
                    absent =[name for name, val in zip(tools, subset) if not val],
                    facecolor=color,
                )
            axes_dict = upset.plot(fig=subfig)

            # Enlarge the intersection-count annotations and the y-axis label
            for txt in axes_dict['intersections'].texts:
                txt.set_fontsize(13)
            axes_dict['intersections'].set_ylabel(
                'precursors', fontsize=15,
            )

            panel = panel_letters[i] if i < len(panel_letters) else f'z{i - 25}'
            subfig.suptitle(
                f'({panel}) {lab}',
                fontsize=18, fontweight='bold', x=0.03, ha='left',
            )

        for i in range(len(labs), nrow * ncol):
            subfigs.flat[i].set_visible(False)

        fig.suptitle(
            'Precursor overlap across search engines — best settings per tool',
            fontsize=22, fontweight='bold', y=1.02,
        )
        legend_patches = [
            mpatches.Patch(facecolor='steelblue',  label='unique to one tool'),
            mpatches.Patch(facecolor='#888888',    label='all but one tool'),
            mpatches.Patch(facecolor='darkorange', label='all tools'),
        ]
        fig.legend(handles=legend_patches, loc='lower center', ncol=3,
                   fontsize=16, frameon=False, bbox_to_anchor=(0.5, -0.01))

        if save_dir:
            save_figure(fig, os.path.join(save_dir, filename))
        plt.show() if show else plt.close(fig)

    finally:
        mpl.rcParams.update(old_rcparams)


def plot_full_upsetplot_panel(
    slices,
    labs,
    ordering_tools,
    ncol:        int        = 2,
    max_upsets:  int | None = None,
    save_dir:    str | None = None,
    filename:    str        = 'upsetplot_panel_full',
    show:        bool       = True,
):
    """Full UpSet panel — every intersection per lab, two plots per row.

    Unlike :func:`plot_upsetplot_panel`, this shows all observed intersections
    (not only degree-1, degree-(n-1) and degree-n). Intersections are sorted by
    cardinality (largest first) so the panel stays readable, and bars are
    coloured by degree (steelblue = unique, #888888 = intermediate,
    darkorange = all tools).

    Parameters
    ----------
    slices:
        Dict mapping lab name to ``{'ups': ..., 'boxplot': ...}``.
    labs:
        Ordered list of lab names to plot (one panel each).
    ordering_tools:
        Tool display order for the category rows.
    ncol:
        Number of UpSet plots per row (default 2).
    max_upsets:
        If set, keep only the ``max_upsets`` largest intersections per lab.
        ``None`` keeps every intersection.
    """
    import matplotlib as mpl

    labs  = list(labs)
    nrow  = ceil(len(labs) / ncol)

    old_rcparams = mpl.rcParams.copy()
    mpl.rcParams.update({
        'font.size':        15,
        'axes.labelsize':   16,
        'axes.titlesize':   18,
        'xtick.labelsize':  14,
        'ytick.labelsize':  14,
    })

    try:
        fig     = plt.figure(figsize=(13 * ncol, 5.5 * nrow))
        subfigs = fig.subfigures(nrow, ncol, hspace=0.08, wspace=0.05)
        subfigs_flat = np.atleast_1d(subfigs).ravel()

        panel_letters = string.ascii_lowercase

        for i, lab in enumerate(labs):
            subfig = subfigs_flat[i]
            ups    = slices[lab]['ups']
            tools  = list(ups.index.names)
            n_t    = len(tools)

            # Sort every intersection by size (largest first), optionally capped
            sizes        = ups.groupby(level=list(range(n_t))).size()
            sorted_sizes = sizes.sort_values(ascending=False)
            if max_upsets is not None:
                sorted_sizes = sorted_sizes.head(max_upsets)
            ups_ordered = ups.loc[sorted_sizes.index]

            upset = upsetplot.UpSet(
                ups_ordered,
                sort_by='cardinality',
                sort_categories_by='input',
                show_counts=False,
                totals_plot_elements=3,
                element_size=None,
            )
            for subset in sorted_sizes.index:
                subset = subset if isinstance(subset, tuple) else (subset,)
                degree = sum(subset)
                if degree == 1:
                    color = 'steelblue'
                elif degree == n_t:
                    color = 'darkorange'
                else:
                    color = '#888888'
                upset.style_subsets(
                    present=[name for name, val in zip(tools, subset) if val],
                    absent =[name for name, val in zip(tools, subset) if not val],
                    facecolor=color,
                )
            axes_dict = upset.plot(fig=subfig)

            for txt in axes_dict['intersections'].texts:
                txt.set_fontsize(12)
            axes_dict['intersections'].set_ylabel('precursors', fontsize=15)

            panel = panel_letters[i] if i < len(panel_letters) else f'z{i - 25}'
            cap_note = f' (top {max_upsets})' if max_upsets is not None else ''
            subfig.suptitle(
                f'({panel}) {lab}{cap_note}',
                fontsize=18, fontweight='bold', x=0.03, ha='left',
            )

        for i in range(len(labs), len(subfigs_flat)):
            subfigs_flat[i].set_visible(False)

        fig.suptitle(
            'Precursor overlap across search engines — all intersections',
            fontsize=22, fontweight='bold', y=1.01,
        )
        legend_patches = [
            mpatches.Patch(facecolor='steelblue',  label='unique to one tool'),
            mpatches.Patch(facecolor='#888888',    label='shared by a subset'),
            mpatches.Patch(facecolor='darkorange', label='all tools'),
        ]
        fig.legend(handles=legend_patches, loc='lower center', ncol=3,
                   fontsize=16, frameon=False, bbox_to_anchor=(0.5, -0.005))

        if save_dir:
            save_figure(fig, os.path.join(save_dir, filename))
        plt.show() if show else plt.close(fig)

    finally:
        mpl.rcParams.update(old_rcparams)


def plot_boxplot_panel(
    slices,
    labs,
    ordering_tools,
    nrow:      int        = 3,
    ncol:      int        = 3,
    ylim:      tuple      = (-0.1, 3),
    save_dir:  str | None = None,
    filename:  str        = 'boxplot_panel_3x3',
    show:      bool       = True,
):
    """3×3 boxplot panel — absolute log2 FC error per tool and precursor subset.

    Two hue levels: 'unique to tool' and 'all tools'. Default seaborn palette.
    y-axis is capped at ``ylim[1]`` to keep the main distribution visible.

    Parameters
    ----------
    slices:
        Dict mapping lab name to ``{'ups': ..., 'boxplot': ...}``.
    labs:
        Ordered list of lab names to plot.
    ordering_tools:
        Controls the x-axis tool order within each panel.
    """
    fig, axes = plt.subplots(nrow, ncol, figsize=(18, 15), constrained_layout=True)

    panel_letters = string.ascii_lowercase

    for i, (lab, ax) in enumerate(zip(labs, axes.flat)):
        bpd        = slices[lab]['boxplot'].copy()
        tools_here = [t for t in ordering_tools if t in bpd['tool'].unique()]

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sns.boxplot(
                data       = bpd,
                x          = 'tool',
                y          = 'AbsQuantError',
                hue        = 'Quantified per software',
                order      = tools_here,
                hue_order  = ['unique to tool', 'all tools'],
                ax         = ax,
                width      = 0.65,
                linewidth  = 1.0,
                flierprops = dict(marker='.', markersize=1.5, alpha=0.3),
            )

        ax.set_ylim(*ylim)
        panel = panel_letters[i] if i < len(panel_letters) else f'z{i - 25}'
        ax.set_title(f'({panel}) {lab}', loc='left', fontweight='bold', fontsize=16)
        ax.set_xlabel('')
        ax.set_ylabel('|log2 FC error|' if i % ncol == 0 else '', fontsize=16)
        ax.tick_params(axis='x', rotation=45, labelsize=13)
        ax.tick_params(axis='y', labelsize=13)
        ax.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.35)
        ax.spines[['top', 'right']].set_visible(False)

        lgd = ax.get_legend()
        if i == 0:
            lgd.set_title('')
            for t in lgd.get_texts():
                t.set_fontsize(13)
        elif lgd:
            lgd.remove()

    for i in range(len(labs), nrow * ncol):
        axes.flat[i].set_visible(False)

    fig.suptitle(
        'Quantification error per search engine — best settings (MBR)',
        fontsize=22, fontweight='bold',
    )

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


# ── Supplemental figure functions ─────────────────────────────────────────────

def plot_upsetplot_panel_supplemental(
    slices,
    labs,
    ordering_tools,
    nrow:      int        = 3,
    ncol:      int        = 3,
    save_dir:  str | None = None,
    filename:  str        = 'upsetplot_panel_3x3',
    show:      bool       = True,
    sort_by_size: bool    = False,
    max_upsets:   int | None = None,
):
    """Publication-ready UpSet panel grid with fixed layout scaling."""
    
    # Scale up global text sizes temporarily so they remain highly legible 
    # even on a large canvas space
    import matplotlib as mpl
    old_rcParams = mpl.rcParams.copy()
    
    mpl.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 14,
        'axes.titlesize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })

    try:
        # Keep a larger canvas size so layout rows don't overlap vertically
        fig     = plt.figure(figsize=(24, 18)) 
        subfigs = fig.subfigures(nrow, ncol, hspace=0.06, wspace=0.04)
        
        abc_labels = string.ascii_lowercase

        for i, lab in enumerate(labs[:nrow * ncol]):
            subfig = subfigs.flat[i]
            ups    = slices[lab]['ups']
            tools  = list(ups.index.names)
            n_t    = len(tools)

            if sort_by_size:
                sizes = ups.groupby(level=list(range(n_t))).size()
                sorted_sizes = sizes.sort_values(ascending=False)
                if max_upsets is not None:
                    sorted_sizes = sorted_sizes.head(max_upsets)
                
                ups_ordered = ups.loc[sorted_sizes.index]
                sort_by_param = 'cardinality'
            else:
                idx_df       = ups.index.to_frame()
                degrees      = idx_df.sum(axis=1)
                keep_degrees = {1, n_t} if n_t <= 2 else {1, n_t - 1, n_t}
                ups_filt     = ups[degrees.isin(keep_degrees)]

                tools_in_order = [t for t in ordering_tools if t in tools]
                desired = []
                for tool in tools_in_order:
                    desired.append(tuple(t == tool for t in tools))
                if n_t > 2:
                    for tool in reversed(tools_in_order):
                        desired.append(tuple(t != tool for t in tools))
                desired.append(tuple(True for _ in tools))
                desired     = [tup for tup in desired if tup in ups_filt.index]
                ups_ordered = ups_filt.loc[desired]
                
                sort_by_param = 'input'

            # Leave element_size=None to completely avoid the get_figwidth SubFigure bug
            upset = upsetplot.UpSet(
                ups_ordered,
                sort_by=sort_by_param,
                sort_categories_by='input',
                show_counts=False,
                totals_plot_elements=3,
                element_size=None,  
                subset_size='count' if isinstance(ups_ordered, pd.DataFrame) else 'auto'
            )
            
            for subset in ups_ordered.index.unique():
                degree = sum(subset)
                if degree == 1:
                    color = 'steelblue'
                elif degree == n_t:
                    color = 'darkorange'
                else:
                    color = '#888888'
                    
                upset.style_subsets(
                    present=[name for name, val in zip(tools, subset) if val],
                    absent =[name for name, val in zip(tools, subset) if not val],
                    facecolor=color,
                )
                
            upset.plot(fig=subfig)
            
            # Label panel letter tokens clearly
            panel_letter = abc_labels[i] if i < len(abc_labels) else f"z{i-25}"
            subfig.suptitle(
                f"({panel_letter}) {lab}", 
                fontsize=14, 
                fontweight='bold', 
                x=0.03, 
                y=0.98,
                ha='left'
            )

        for i in range(len(labs), nrow * ncol):
            subfigs.flat[i].set_visible(False)

        fig.suptitle(
            'Precursor overlap across search engines — best settings per tool',
            fontsize=16, fontweight='bold', y=1.01,
        )
        
        legend_patches = [
            mpatches.Patch(facecolor='steelblue',  label='unique to one tool'),
            mpatches.Patch(facecolor='#888888',    label='all but one tool'),
            mpatches.Patch(facecolor='darkorange', label='all tools'),
        ]
        fig.legend(
            handles=legend_patches, loc='lower center', ncol=3,
            fontsize=13, frameon=False, bbox_to_anchor=(0.5, 0.01)
        )

        if save_dir:
            fig.savefig(os.path.join(save_dir, filename + '.png'), dpi=300, bbox_inches='tight')
        
        plt.show() if show else plt.close(fig)

    finally:
        # Always safely restore original font states for other graphs in your runtime session
        mpl.rcParams.update(old_rcParams)


def plot_boxplot_panel_supplemental(
    slices,
    labs,
    ordering_tools,
    nrow: int = 2,
    ncol: int = 4,
    ylim: tuple = (-0.1, 3),
    save_dir: str | None = None,
    filename: str = 'boxplot_panel',
    show: bool = True,
):
    """Publication-ready boxplot panel."""

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(20, 8),
        constrained_layout=True,
    )

    labels = string.ascii_lowercase

    for i, (lab, ax) in enumerate(zip(labs, axes.flat)):

        bpd = slices[lab]['boxplot'].copy()

        tools_here = [
            t for t in ordering_tools
            if t in bpd['tool'].unique()
        ]

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            sns.boxplot(
                data=bpd,
                x='tool',
                y='AbsQuantError',
                hue='Quantified per software',
                order=tools_here,
                hue_order=[
                    'unique to tool',
                    'all tools',
                ],
                ax=ax,
                width=0.65,
                linewidth=0.8,
                flierprops=dict(
                    marker='.',
                    markersize=1.5,
                    alpha=0.3,
                ),
            )

        ax.set_ylim(*ylim)

        ax.set_title(
            f'({labels[i]}) {lab}',
            loc='left',
            fontsize=11,
            fontweight='bold',
        )

        ax.set_xlabel('')

        if i % ncol == 0:
            ax.set_ylabel(
                '|log2 FC error|',
                fontsize=10,
            )
        else:
            ax.set_ylabel('')

        ax.tick_params(
            axis='x',
            rotation=45,
            labelsize=8,
        )

        ax.tick_params(
            axis='y',
            labelsize=8,
        )

        ax.grid(
            axis='y',
            linestyle='--',
            linewidth=0.5,
            alpha=0.35,
        )

        ax.spines[['top', 'right']].set_visible(False)

        lgd = ax.get_legend()

        if i == 0 and lgd:

            lgd.set_title('')

            for text in lgd.get_texts():
                text.set_fontsize(8)

        elif lgd:
            lgd.remove()

    for i in range(len(labs), nrow * ncol):
        axes.flat[i].set_visible(False)

    fig.suptitle(
        'Quantification error per search engine — best settings (MBR)',
        fontsize=14,
        fontweight='bold',
    )

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))

    plt.show() if show else plt.close(fig)