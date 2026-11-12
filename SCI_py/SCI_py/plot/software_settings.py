"""
Scatter plots comparing search-engine settings on the benchmark data.

Part of Notebooks/1_QCQuant/02_Settings_comparison: each point is one setting,
with mean absolute error (epsilon) on x and precursor count on y. Recommended
settings are selected with ``SCI_py.benchmark.filters.filter_recommended``.
"""
import os
import warnings

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from ..benchmark.filters import filter_recommended
from .export import save_figure

_MARKERS      = ['o', 's', 'D', '^', 'v', 'P', 'X', '*', 'h']
_INSTR_MARKER = {'Astral': 'o', 'timsTOF': 's'}


# ── Axis-level primitive ──────────────────────────────────────────────────────

def _scatter_on_ax(ax, df, color_column, color_map,
                   x_col='mean_abs_epsilon', shape_column=None,
                   symbol_map=None, size=55):
    """Draw benchmark scatter points on an existing Axes."""
    if shape_column is not None and symbol_map is None:
        uniq       = sorted(df[shape_column].dropna().unique())
        symbol_map = {v: _MARKERS[i % len(_MARKERS)] for i, v in enumerate(uniq)}

    for group in sorted(df[color_column].dropna().unique()):
        sub   = df[df[color_column] == group]
        color = color_map.get(group, 'grey')
        if shape_column is not None:
            for sv, ss in sub.groupby(shape_column):
                ax.scatter(ss[x_col], ss['nr_prec'],
                           c=color, marker=symbol_map.get(sv, 'o'),
                           s=size, alpha=0.8)
        else:
            ax.scatter(sub[x_col], sub['nr_prec'],
                       c=color, marker='o', s=size, alpha=0.8, label=group)

    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, color='lightgrey', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel('Mean absolute ε', fontsize=8)
    ax.set_ylabel('Precursors', fontsize=8)
    ax.tick_params(labelsize=7)


# ── Figure-level functions ────────────────────────────────────────────────────

def plot_scatter_all_tools(
    datapoints,
    color_map,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_all_tools',
    show:      bool       = True,
):
    """Section 4.1 — all tools, all datasets; colour = search engine."""
    fig, ax = plt.subplots(figsize=(8, 5))
    _scatter_on_ax(ax, datapoints, color_column='software_name', color_map=color_map)

    present = sorted(datapoints['software_name'].unique())
    handles = [plt.scatter([], [], c=color_map.get(e, 'grey'), s=55, alpha=0.8, label=e)
               for e in present]
    ax.legend(handles=handles, title='Search engine',
              bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize=8)
    ax.set_title('All tools — all datasets')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


def plot_scatter_per_tool(
    datapoints,
    lab_colors,
    lab_order,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_per_tool',
    show:      bool       = True,
):
    """Section 4.2 — one panel per search engine; colour = lab."""
    tools  = sorted(datapoints['software_name'].unique())
    ncols  = 3
    nrows  = -(-len(tools) // ncols)

    uniq_comments = sorted(datapoints['comments'].dropna().unique())
    symbol_map    = {v: _MARKERS[i % len(_MARKERS)] for i, v in enumerate(uniq_comments)}

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                             constrained_layout=True)
    axes_flat = axes.flatten()

    for i, tool in enumerate(tools):
        sub = datapoints[datapoints['software_name'] == tool]
        _scatter_on_ax(axes_flat[i], sub, color_column='lab', color_map=lab_colors,
                       shape_column='comments', symbol_map=symbol_map)
        axes_flat[i].set_title(tool, fontsize=10, fontweight='bold')

        present_shapes = sub['comments'].dropna().unique()
        handles_shape  = [
            plt.scatter([], [], c='grey', s=55, alpha=0.8,
                        marker=symbol_map[s], label=s)
            for s in sorted(present_shapes)
        ]
        if handles_shape:
            axes_flat[i].legend(handles=handles_shape, title='Comments',
                                fontsize=7, title_fontsize=7, loc='upper left')

    for j in range(len(tools), len(axes_flat)):
        axes_flat[j].set_visible(False)

    present_labs = [l for l in lab_order if l in datapoints['lab'].values]
    handles_lab  = [
        plt.scatter([], [], c=lab_colors.get(l, 'grey'), s=55, alpha=0.8, label=l)
        for l in present_labs
    ]
    fig.legend(handles=handles_lab, title='Lab', loc='lower right',
               bbox_to_anchor=(1.0, 0.0), fontsize=8)
    fig.suptitle('Per-tool benchmark — all settings, colour by lab', fontweight='bold')

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


def plot_scatter_instrument_split(
    datapoints,
    color_map,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_instrument_split',
    show:      bool       = True,
):
    """Section 4.3 — Astral vs timsTOF; colour = search engine."""
    dtps_astral  = datapoints[datapoints['lab'].str.contains('Astral')]
    dtps_timstof = datapoints[datapoints['lab'].str.contains('timsTOF')]

    fig, (ax_a, ax_t) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    _scatter_on_ax(ax_a, dtps_astral,  color_column='software_name',
                   color_map=color_map, shape_column='lab')
    _scatter_on_ax(ax_t, dtps_timstof, color_column='software_name',
                   color_map=color_map, shape_column='lab')

    ax_a.set_title('Astral',  fontweight='bold')
    ax_t.set_title('timsTOF', fontweight='bold')

    present = sorted(datapoints['software_name'].unique())
    handles = [plt.scatter([], [], c=color_map.get(e, 'grey'), s=55, alpha=0.8, label=e)
               for e in present]
    fig.legend(handles=handles, title='Search engine',
               bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0, fontsize=8)
    fig.suptitle('Astral vs timsTOF — all tools, all settings', fontweight='bold')

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


def plot_scatter_tool_instrument_grid(
    datapoints,
    lab_colors,
    lab_order,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_tool_instrument_grid',
    show:      bool       = True,
):
    """Section 4.4 — rows = tools, cols = Astral/timsTOF; colour = lab."""
    tools   = sorted(datapoints['software_name'].unique())
    n_tools = len(tools)

    uniq_comments = sorted(datapoints['comments'].dropna().unique())
    symbol_map    = {v: _MARKERS[i % len(_MARKERS)] for i, v in enumerate(uniq_comments)}

    fig, axes = plt.subplots(n_tools, 2, figsize=(10, 3.5 * n_tools),
                             constrained_layout=True)

    for row, tool in enumerate(tools):
        sub_a = datapoints[(datapoints['software_name'] == tool) &
                            datapoints['lab'].str.contains('Astral')]
        sub_t = datapoints[(datapoints['software_name'] == tool) &
                            datapoints['lab'].str.contains('timsTOF')]

        _scatter_on_ax(axes[row, 0], sub_a, color_column='lab', color_map=lab_colors,
                       shape_column='comments', symbol_map=symbol_map)
        _scatter_on_ax(axes[row, 1], sub_t, color_column='lab', color_map=lab_colors,
                       shape_column='comments', symbol_map=symbol_map)
        axes[row, 0].set_ylabel(tool, fontsize=9, fontweight='bold')

        present_shapes = sorted(set(
            sub_a['comments'].dropna().tolist() + sub_t['comments'].dropna().tolist()
        ))
        handles_shape = [
            plt.scatter([], [], c='grey', s=55, alpha=0.8,
                        marker=symbol_map[s], label=s)
            for s in present_shapes
        ]
        if handles_shape:
            leg = axes[row, 0].legend(
                handles=handles_shape, title=f'{tool} — comments',
                fontsize=7, title_fontsize=7, loc='upper left',
                bbox_to_anchor=(1.3, 1),
                bbox_transform=axes[row, 1].transAxes,
            )
            axes[row, 0].add_artist(leg)

    axes[0, 0].set_title('Astral',  fontweight='bold')
    axes[0, 1].set_title('timsTOF', fontweight='bold')

    present_labs = [l for l in lab_order if l in datapoints['lab'].values]
    handles_lab  = [
        plt.scatter([], [], c=lab_colors.get(l, 'grey'), s=55, alpha=0.8, label=l)
        for l in present_labs
    ]
    fig.legend(handles=handles_lab, title='Lab',
               bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0, fontsize=8)
    fig.suptitle('Per-tool × instrument — all settings, colour by lab', fontweight='bold')

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


def plot_scatter_recommended(
    datapoints,
    recommended,
    recommended_mbr,
    color_map,
    lab_order,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_recommended',
    show:      bool       = True,
):
    """Section 4.5 — two panels: no-MBR and MBR recommended settings."""
    dtps_rec     = filter_recommended(datapoints, recommended)
    dtps_rec_mbr = filter_recommended(datapoints, recommended_mbr)

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    _scatter_on_ax(ax_l, dtps_rec,     color_column='software_name',
                   color_map=color_map, shape_column='lab')
    _scatter_on_ax(ax_r, dtps_rec_mbr, color_column='software_name',
                   color_map=color_map, shape_column='lab')

    ax_l.set_title('Recommended (no MBR)', fontweight='bold')
    ax_r.set_title('MBR / 2nd-pass',       fontweight='bold')

    present    = sorted(datapoints['software_name'].unique())
    color_handles = [plt.scatter([], [], c=color_map.get(e, 'grey'),
                                 s=55, alpha=0.8, label=e) for e in present]
    color_leg = ax_r.legend(handles=color_handles, title='Search engine',
                             bbox_to_anchor=(1.02, 1), loc='upper left',
                             borderaxespad=0, fontsize=8)
    ax_r.add_artist(color_leg)

    present_labs  = [l for l in lab_order if l in datapoints['lab'].values]
    shape_handles = [
        mlines.Line2D([], [], marker=_MARKERS[i % len(_MARKERS)],
                      color='grey', linestyle='none', markersize=7, label=l)
        for i, l in enumerate(present_labs)
    ]
    ax_r.legend(handles=shape_handles, title='Lab',
                bbox_to_anchor=(1.02, 0.4), loc='upper left',
                borderaxespad=0, fontsize=8)

    fig.suptitle('Recommended settings — all tools and labs', fontweight='bold')

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)

    print(f'Recommended (no MBR): {len(dtps_rec)} points')
    print(f'MBR / 2nd-pass:       {len(dtps_rec_mbr)} points')


def plot_scatter_recommended_split(
    datapoints,
    recommended,
    recommended_mbr,
    color_map,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_recommended_split',
    show:      bool       = True,
):
    """Section 4.6 — recommended settings split by Astral/timsTOF; one figure per rec variant."""
    instrument_subsets = {
        'Astral':   datapoints[datapoints['lab'].str.contains('Astral')],
        'timsTOF': datapoints[datapoints['lab'].str.contains('timsTOF')],
    }

    for rec_label, rec_dict, suffix in [
        ('Recommended (no MBR)', recommended,     f'{filename}_no_mbr'),
        ('MBR / 2nd-pass',       recommended_mbr, f'{filename}_mbr'),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

        for ax, (instr, sub) in zip(axes, instrument_subsets.items()):
            filtered = filter_recommended(sub, rec_dict)
            _scatter_on_ax(ax, filtered, color_column='software_name',
                           color_map=color_map, shape_column='lab')
            ax.set_title(instr, fontweight='bold')

        present = sorted(datapoints['software_name'].unique())
        handles = [plt.scatter([], [], c=color_map.get(e, 'grey'),
                               s=55, alpha=0.8, label=e) for e in present]
        fig.legend(handles=handles, title='Search engine',
                   bbox_to_anchor=(1.01, 1), loc='upper left',
                   borderaxespad=0, fontsize=8)
        fig.suptitle(f'{rec_label} — Astral vs timsTOF', fontweight='bold')

        if save_dir:
            save_figure(fig, os.path.join(save_dir, suffix))
        plt.show() if show else plt.close(fig)


def plot_scatter_averaged(
    averaged,
    color_map,
    save_dir:  str | None = None,
    filename:  str        = 'scatter_averaged',
    show:      bool       = True,
):
    """Section 4.7 — instrument-averaged performance; colour = engine, shape = instrument."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for _, row in averaged.iterrows():
        ax.scatter(
            row['mean_abs_epsilon'], row['nr_prec'],
            c      = color_map.get(row['software_name'], 'grey'),
            marker = _INSTR_MARKER.get(row['instrument_type'], 'o'),
            s      = 160,
            alpha  = 0.9,
            zorder = 3,
        )

    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, color='lightgrey', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel('Mean absolute ε (averaged across labs)', fontsize=9)
    ax.set_ylabel('Precursors (averaged across labs)', fontsize=9)
    ax.set_title('Instrument-averaged performance — best settings per tool',
                 fontweight='bold')

    present      = sorted(averaged['software_name'].unique())
    color_handles = [
        plt.scatter([], [], c=color_map.get(e, 'grey'), s=80, alpha=0.9, label=e)
        for e in present
    ]
    color_leg = ax.legend(handles=color_handles, title='Search engine',
                           bbox_to_anchor=(1.02, 1), loc='upper left',
                           borderaxespad=0, fontsize=8)
    ax.add_artist(color_leg)

    shape_handles = [
        plt.scatter([], [], c='grey', marker=m, s=80, alpha=0.9, label=instr)
        for instr, m in _INSTR_MARKER.items()
    ]
    ax.legend(handles=shape_handles, title='Instrument',
              bbox_to_anchor=(1.02, 0.4), loc='upper left',
              borderaxespad=0, fontsize=8)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)
