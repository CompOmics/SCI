import os
import warnings
from math import ceil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from ..constants.plotting import SOFTWARE_COLORS
from .utilities import _save


# ── KDE of log2 fold changes ──────────────────────────────────────────────────

def plot_species_distributions(
    intermediates,
    species_config,
    color_map=None,
    normalize_kde=False,
    xlim_per_species=None,
    save_dir: str | None = None,
    filename: str = 'species_distributions',
    show: bool = True,
):
    """KDE distributions of log2 fold changes per lab and species.

    Parameters
    ----------
    intermediates:
        Long-format DataFrame with columns ``lab``, ``species``,
        ``search_engine``, ``log2_A_vs_B``.
    species_config:
        ``{species: {'A_vs_B': float, 'color': str}}``.
    color_map:
        ``{search_engine: color}``. Falls back to the mpl cycle.
    normalize_kde:
        If True, normalise each KDE to integrate to 1.
    xlim_per_species:
        ``{species: [xmin, xmax]}``.
    """
    lab_list = sorted(intermediates['lab'].unique())
    species_list = sorted(intermediates['species'].unique())
    software_list = sorted(intermediates['search_engine'].unique())

    if color_map is None:
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        color_map = {
            s: cycle[i % len(cycle)]
            for i, s in enumerate(software_list)
        }

    n_labs, n_species = len(lab_list), len(species_list)

    fig, axes = plt.subplots(
        n_labs,
        n_species,
        figsize=(4.8 * n_species, 2.4 * n_labs),
        constrained_layout=True,
        sharey='row',
    )

    if n_labs == 1:
        axes = axes[np.newaxis, :]
    if n_species == 1:
        axes = axes[:, np.newaxis]

    legend_handles = {}

    for row_idx, lab in enumerate(lab_list):
        for col_idx, species in enumerate(species_list):

            ax = axes[row_idx, col_idx]

            subset = intermediates[
                (intermediates['lab'] == lab)
                & (intermediates['species'] == species)
            ]

            for software in software_list:

                values = (
                    subset.loc[
                        subset['search_engine'] == software,
                        'log2_A_vs_B'
                    ]
                    .dropna()
                    .values
                )

                if len(values) < 2:
                    continue

                color = color_map.get(software, 'black')

                kde = stats.gaussian_kde(values)

                x_grid = np.linspace(
                    values.min() - 1,
                    values.max() + 1,
                    300,
                )

                y_kde = kde(x_grid)

                if not normalize_kde:
                    y_kde *= len(values)

                line, = ax.plot(
                    x_grid,
                    y_kde,
                    color=color,
                    linewidth=1.7,
                )

                legend_handles.setdefault(software, line)

            expected_log2 = np.log2(
                species_config[species]['A_vs_B']
            )

            ax.axvline(
                expected_log2,
                color=species_config[species]['color'],
                linestyle='--',
                linewidth=1.4,
                alpha=0.9,
            )

            # Light gridlines for easier visual comparison
            ax.grid(
                True,
                which='major',
                linestyle='--',
                linewidth=0.6,
                alpha=0.35,
            )

            if xlim_per_species and species in xlim_per_species:
                ax.set_xlim(xlim_per_species[species])

            if row_idx == 0:
                ax.set_title(
                    species,
                    fontsize=12,
                    fontweight='bold',
                    pad=10,
                )

            if col_idx == 0:
                ax.set_ylabel(
                    lab,
                    fontsize=22,
                )

            if row_idx == n_labs - 1:
                ax.set_xlabel(
                    'log2 FC (A vs B)',
                    fontsize=22,
                )

            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=14)

    fig.legend(
        handles=list(legend_handles.values()),
        labels=list(legend_handles.keys()),
        title='Search engine',
        bbox_to_anchor=(1.01, 1),
        loc='upper left',
        borderaxespad=0,
        fontsize=9,
        title_fontsize=10,
    )

    fig.suptitle(
        'Log2 fold change distributions per lab and species',
        fontsize=15,
        fontweight='bold',
        y=1.02,
    )

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))

    plt.show() if show else plt.close(fig)


# ── MA-plots per tool ─────────────────────────────────────────────────────────

def plot_ma_per_tool(
    intermediates,
    species_config,
    species_colors,
    lab,
    tools                  = None,
    ncols:           int   = 3,
    xlim:            tuple = (-4, 5),
    use_rank:        bool  = True,
    point_size:      int   = 4,
    alpha:           float = 0.25,
    max_points_per_species = 20_000,
    random_state:    int   = 0,
    save_dir:  str | None  = None,
    filename:  str         = 'ma_plots',
    show:      bool        = True,
):
    """MA-plot per tool for a single lab: x = measured log2 FC, coloured by species.

    Parameters
    ----------
    use_rank:
        If True, y = percentile rank (0–100) of mean log intensity within the tool/lab.
        If False, y = raw mean log intensity.
    max_points_per_species:
        Down-sample each species to at most this many points for speed.
    """
    df = intermediates[intermediates['lab'] == lab].copy()
    if df.empty:
        raise ValueError(f'No intermediates for lab {lab!r}.')
    df['abundance'] = 0.5 * (df['log_Intensity_mean_A'] + df['log_Intensity_mean_B'])

    if tools is None:
        tools = sorted(df['search_engine'].unique())
    species_list = [s for s in species_colors if s in df['species'].unique()]

    nrows = ceil(len(tools) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.6 * ncols, 4.0 * nrows),
        sharex=True, sharey=use_rank, squeeze=False,
    )
    axes_flat = axes.ravel()

    rng = np.random.default_rng(random_state)
    for ax, tool in zip(axes_flat, tools):
        sub = df[df['search_engine'] == tool].copy()
        if use_rank:
            sub['y'] = sub['abundance'].rank(pct=True) * 100
        else:
            sub['y'] = sub['abundance']

        for sp in species_list:
            d = sub[sub['species'] == sp].dropna(subset=['log2_A_vs_B', 'y'])
            if len(d) > max_points_per_species:
                d = d.iloc[rng.choice(len(d), max_points_per_species, replace=False)]
            ax.scatter(d['log2_A_vs_B'], d['y'], s=point_size, alpha=alpha,
                       color=species_colors[sp], label=sp, linewidths=0, rasterized=True)

        for sp in species_list:
            exp = np.log2(species_config[sp]['A_vs_B'])
            ax.axvline(exp, color=species_colors[sp], linestyle=':', linewidth=1.6, zorder=5)

        ax.set_title(f'{tool}  (n={len(sub):,})', fontweight='bold', fontsize=10)
        ax.set_xlim(xlim)
        ax.spines[['top', 'right']].set_visible(False)

    for ax in axes_flat[len(tools):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel('measured log2 FC (A vs B)')
    y_label = ('abundance rank (percentile, within tool)' if use_rank
               else 'mean log intensity (A, B)')
    for ax in axes[:, 0]:
        ax.set_ylabel(y_label)

    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=species_colors[s],
                   label=s, markersize=7)
        for s in species_list
    ]
    fig.legend(handles=handles, title='Species', bbox_to_anchor=(1.0, 0.5),
               loc='center left', frameon=False)
    fig.suptitle(
        f'MA-plots per tool - abundance vs measured fold change  (lab: {lab})',
        fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 0.93, 1])

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


# ── Pairwise tool comparison ──────────────────────────────────────────────────

def _abundance(df):
    return 0.5 * (df['log_Intensity_mean_A'] + df['log_Intensity_mean_B'])


def compare_tools_error(
    inter_a,
    inter_b,
    name_a,
    name_b,
    species_colors,
    lab,
    species_order  = ('HUMAN', 'YEAST', 'ECOLI'),
    match_on       = ('precursor ion',),
    n_bins:  int   = 12,
    color_a        = None,
    color_b        = None,
    min_per_bin    = 20,
    save_dir: str | None = None,
    filename: str        = 'pairwise_error',
    show:     bool       = True,
):
    """Compare quantification error of two tools on matching precursors within a single lab.

    Returns the merged per-precursor DataFrame. Prints per-species summaries.

    The x-axis uses an abundance rank (percentile over the matched precursor set) so that
    the very different raw intensity scales between tools are comparable.
    """
    match_on = list(match_on)
    col_a = color_a or SOFTWARE_COLORS.get(name_a, '#1f77b4')
    col_b = color_b or SOFTWARE_COLORS.get(name_b, '#d62728')

    inter_a = inter_a[inter_a['lab'] == lab]
    inter_b = inter_b[inter_b['lab'] == lab]
    if inter_a.empty or inter_b.empty:
        raise ValueError(f'No data for lab {lab!r} in one of the tools.')

    a = inter_a[match_on + ['epsilon', 'species',
                             'log_Intensity_mean_A', 'log_Intensity_mean_B']].copy()
    b = inter_b[match_on + ['epsilon',
                             'log_Intensity_mean_A', 'log_Intensity_mean_B']].copy()
    a['abund_a'] = _abundance(a)
    b['abund_b'] = _abundance(b)

    m = a[match_on + ['species', 'epsilon', 'abund_a']].merge(
        b[match_on + ['epsilon', 'abund_b']], on=match_on, suffixes=('_a', '_b'),
    )
    m = m.dropna(subset=['epsilon_a', 'epsilon_b', 'abund_a', 'abund_b'])
    m['abs_err_a'] = m['epsilon_a'].abs()
    m['abs_err_b'] = m['epsilon_b'].abs()
    rank_a = m['abund_a'].rank(pct=True) * 100
    rank_b = m['abund_b'].rank(pct=True) * 100
    m['abundance'] = 0.5 * (rank_a + rank_b)
    m['delta']     = m['abs_err_a'] - m['abs_err_b']

    species_list = [s for s in species_order if s in m['species'].unique()]
    n_sp = len(species_list)

    fig, axes = plt.subplots(2, n_sp, figsize=(5 * n_sp, 8),
                             squeeze=False, sharex='col')

    print(f'{name_a} vs {name_b}: {len(m):,} matched precursors')
    for col, sp in enumerate(species_list):
        d = m[m['species'] == sp].copy()
        d['bin'] = pd.qcut(d['abundance'], n_bins, duplicates='drop')
        g    = d.groupby('bin', observed=True)
        x    = g['abundance'].median().values
        keep = g.size().values >= min_per_bin

        mean_a, sem_a = g['abs_err_a'].mean().values, g['abs_err_a'].sem().values
        mean_b, sem_b = g['abs_err_b'].mean().values, g['abs_err_b'].sem().values
        med_delta     = g['delta'].median().values

        # Top: mean |epsilon| per tool
        ax = axes[0, col]
        ax.plot(x[keep], mean_a[keep], color=col_a, label=name_a, linewidth=2)
        ax.fill_between(x[keep], (mean_a - sem_a)[keep], (mean_a + sem_a)[keep],
                        color=col_a, alpha=0.2)
        ax.plot(x[keep], mean_b[keep], color=col_b, label=name_b, linewidth=2)
        ax.fill_between(x[keep], (mean_b - sem_b)[keep], (mean_b + sem_b)[keep],
                        color=col_b, alpha=0.2)
        ax.set_title(f'{sp}  (n={len(d):,})', fontweight='bold')
        if col == 0:
            ax.set_ylabel('mean |error|  (|measured - expected|)')
            ax.legend(frameon=False, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)

        # Bottom: paired delta
        ax = axes[1, col]
        xk, dk = x[keep], med_delta[keep]
        ax.plot(xk, dk, color='black', linewidth=1.5)
        ax.axhline(0, color='grey', linewidth=0.8)
        ax.fill_between(xk, dk, 0, where=dk < 0, color=col_a, alpha=0.45, interpolate=True)
        ax.fill_between(xk, dk, 0, where=dk > 0, color=col_b, alpha=0.45, interpolate=True)
        if col == 0:
            ax.set_ylabel(f'median(|err {name_a}| - |err {name_b}|)')
        ax.set_xlabel('abundance rank (percentile, matched precursors)')
        ax.spines[['top', 'right']].set_visible(False)

        frac_a = (d['delta'] < 0).mean() * 100
        med_a  = d['abs_err_a'].median()
        med_b  = d['abs_err_b'].median()
        print(f'  {sp:6s}: {name_a} more accurate in {frac_a:5.1f}% of matched precursors '
              f'| median |err|: {name_a}={med_a:.3f}  {name_b}={med_b:.3f}')

    fig.suptitle(
        f'Quantification error vs abundance - {name_a} vs {name_b}  (lab: {lab})\n'
        f'(lower is better; bottom panel below 0 = {name_a} wins, above 0 = {name_b} wins)',
        fontweight='bold',
    )
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)
    return m
