"""
Quantitative performance plots for search engines on the benchmark data.

Part of Notebooks/1_QCQuant/03_Performance_comparison: species log2 fold-change
distributions, MA plots, error versus abundance, and pairwise tool error
comparisons. There is no matching metrics module.
"""
import os
from math import ceil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from ..constants.plotting import SOFTWARE_COLORS
from .export import save_figure


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
        save_figure(fig, os.path.join(save_dir, filename))

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
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


# ── MA-plots per tool, pooled over a set of labs ──────────────────────────────

def plot_ma_multilab(
    intermediates,
    species_config,
    species_colors,
    labs,
    tools                  = None,
    ncols:           int   = 3,
    xlim:            tuple = (-4, 5),
    use_rank:        bool  = True,
    point_size:      int   = 4,
    alpha:           float = 0.25,
    max_points_per_species = 20_000,
    random_state:    int   = 0,
    group_name:  str | None = None,
    save_dir:  str | None  = None,
    filename:  str         = 'ma_plots_multilab',
    show:      bool        = True,
):
    """MA-plot per tool, pooling several labs into one panel.

    Same representation as :func:`plot_ma_per_tool` (x = measured log2 FC,
    coloured by species, expected per-species ratios as dotted vertical lines),
    but combining a list of labs (e.g. all Astral labs, or all timsTOF labs).

    Abundances are only comparable within a single lab, so the abundance rank is
    computed **separately within each (lab, tool)** and only then are the points
    pooled and plotted together. This keeps the per-lab error/rank calculation
    intact while showing the instrument group as a single cloud per tool.

    Parameters
    ----------
    labs:
        Iterable of lab names to pool into one figure.
    group_name:
        Optional label for the pooled group (e.g. ``'Astral'``), used in the title.
    use_rank:
        If True, y = percentile rank (0-100) of mean log intensity within each
        (lab, tool). If False, y = raw mean log intensity (only comparable within
        a lab, so ranking is recommended when pooling labs).
    max_points_per_species:
        Down-sample each species (after pooling) to at most this many points.
    """
    labs = list(labs)
    df = intermediates[intermediates['lab'].isin(labs)].copy()
    if df.empty:
        raise ValueError(f'No intermediates for labs {labs!r}.')
    df['abundance'] = 0.5 * (df['log_Intensity_mean_A'] + df['log_Intensity_mean_B'])

    if tools is None:
        tools = sorted(df['search_engine'].unique())
    species_list = [s for s in species_colors if s in df['species'].unique()]

    # Rank computed per (lab, tool) so intensity scales stay comparable, then pooled.
    if use_rank:
        df['y'] = (
            df.groupby(['lab', 'search_engine'])['abundance']
            .rank(pct=True) * 100
        )
    else:
        df['y'] = df['abundance']

    nrows = ceil(len(tools) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.6 * ncols, 4.0 * nrows),
        sharex=True, sharey=use_rank, squeeze=False,
    )
    axes_flat = axes.ravel()

    rng = np.random.default_rng(random_state)
    for ax, tool in zip(axes_flat, tools):
        sub = df[df['search_engine'] == tool]

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
    y_label = ('abundance rank (percentile, within lab/tool)' if use_rank
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

    group_txt = group_name if group_name else f'{len(labs)} labs'
    fig.suptitle(
        f'MA-plots per tool - abundance vs measured fold change  '
        f'({group_txt}: {", ".join(labs)})',
        fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 0.93, 1])

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)


# ── Multi-tool dynamic range, pooled over a set of labs ───────────────────────

def plot_dynamic_range_multi(
    intermediates,
    species_colors,
    labs,
    tools,
    color_map            = None,
    species_order        = ('HUMAN', 'YEAST', 'ECOLI'),
    match_on             = ('precursor ion',),
    n_bins:  int         = 12,
    min_per_bin:  int    = 20,
    group_name:  str | None = None,
    save_dir: str | None = None,
    filename: str        = 'dynamic_range_multi',
    show:     bool       = True,
):
    """Mean absolute error vs abundance rank for a set of tools, pooled over a set of labs.

    Generalises :func:`compare_tools_error` from two tools to an arbitrary list of
    tools. For each lab the precursors are restricted to those quantified by **all**
    requested tools (stringent overlap), the abundance rank is computed within that
    lab, and only then are the matched precursors pooled across labs. Ranks, bins and
    errors are therefore always computed per lab; pooling only affects what is drawn.

    The figure has, per species column:

    - **Top** - mean absolute error ``|epsilon|`` vs abundance rank, one line per
      tool (with SEM band). The lowest line is the most accurate tool at that
      abundance.
    - **Bottom** - a thin band coloured by the tool with the lowest mean error in
      each abundance bin (the winning tool across the dynamic range).

    Parameters
    ----------
    labs:
        Iterable of lab names. Ranks/bins/errors are computed separately per lab
        and pooled for plotting.
    tools:
        Tools to compare (typically the recommended set shown in the MA-plots).
    color_map:
        ``{tool: color}``. Falls back to ``SOFTWARE_COLORS`` then the mpl cycle.
    group_name:
        Optional label for the pooled group (e.g. ``'Astral'``), used in the title.
    """
    labs     = list(labs)
    tools    = list(tools)
    match_on = list(match_on)

    if color_map is None:
        color_map = SOFTWARE_COLORS
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    tool_color = {
        t: color_map.get(t, cycle[i % len(cycle)])
        for i, t in enumerate(tools)
    }

    err_cols = [f'abs_err__{t}' for t in tools]

    per_lab = []
    for lab in labs:
        inter_lab = intermediates[intermediates['lab'] == lab]

        merged = None
        for tool in tools:
            t = inter_lab[inter_lab['search_engine'] == tool][
                match_on + ['species', 'epsilon',
                            'log_Intensity_mean_A', 'log_Intensity_mean_B']
            ].copy()
            t['abund'] = _abundance(t)
            t = t.dropna(subset=['epsilon', 'abund'])
            keep_cols = match_on + ['epsilon', 'abund']
            if merged is None:
                keep_cols = match_on + ['species', 'epsilon', 'abund']
            t = t[keep_cols].rename(
                columns={'epsilon': f'eps__{tool}', 'abund': f'abund__{tool}'}
            )
            merged = t if merged is None else merged.merge(t, on=match_on, how='inner')

        if merged is None or merged.empty:
            continue

        # Ranks computed within this lab over the fully matched precursor set.
        rank_cols = []
        for tool in tools:
            rc = f'rank__{tool}'
            merged[rc] = merged[f'abund__{tool}'].rank(pct=True) * 100
            merged[f'abs_err__{tool}'] = merged[f'eps__{tool}'].abs()
            rank_cols.append(rc)
        merged['abundance'] = merged[rank_cols].mean(axis=1)
        merged['lab'] = lab
        per_lab.append(merged[['lab', 'species', 'abundance'] + err_cols])

    if not per_lab:
        raise ValueError(f'No precursors matched across all tools for labs {labs!r}.')

    m = pd.concat(per_lab, ignore_index=True)

    species_list = [s for s in species_order if s in m['species'].unique()]
    n_sp = len(species_list)

    fig, axes = plt.subplots(
        2, n_sp, figsize=(5 * n_sp, 7), squeeze=False, sharex='col', sharey='row',
        gridspec_kw={'height_ratios': [6, 1]},
    )

    group_txt = group_name if group_name else f'{len(labs)} labs'
    print(f'Dynamic range ({group_txt}): {len(m):,} precursors matched across all '
          f'{len(tools)} tools, pooled over {len(labs)} lab(s)')

    for col, sp in enumerate(species_list):
        d = m[m['species'] == sp].copy()
        d['bin'] = pd.qcut(d['abundance'], n_bins, duplicates='drop')
        g = d.groupby('bin', observed=True)
        x = g['abundance'].median().values
        keep = g.size().values >= min_per_bin

        ax = axes[0, col]
        mean_err = {}
        for tool in tools:
            mean_t = g[f'abs_err__{tool}'].mean().values
            sem_t  = g[f'abs_err__{tool}'].sem().values
            mean_err[tool] = mean_t
            c = tool_color[tool]
            ax.plot(x[keep], mean_t[keep], color=c, label=tool, linewidth=2)
            ax.fill_between(x[keep], (mean_t - sem_t)[keep], (mean_t + sem_t)[keep],
                            color=c, alpha=0.15)
        ax.set_title(f'{sp}  (n={len(d):,})', fontweight='bold')
        if col == 0:
            ax.set_ylabel('mean |error|  (|measured - expected|)')
            ax.legend(frameon=False, fontsize=8, ncol=1)
        ax.spines[['top', 'right']].set_visible(False)

        # Best-tool band: colour each bin by the tool with the lowest mean error.
        band = axes[1, col]
        intervals = list(g.groups.keys())
        err_matrix = np.vstack([mean_err[t] for t in tools])  # (n_tools, n_bins)
        best_idx = np.nanargmin(err_matrix, axis=0)
        for j, iv in enumerate(intervals):
            if not keep[j]:
                continue
            band.axvspan(iv.left, iv.right, color=tool_color[tools[best_idx[j]]],
                         alpha=0.85, linewidth=0)
        band.set_yticks([])
        band.set_ylabel('best', fontsize=8, rotation=0, ha='right', va='center')
        band.set_xlabel('abundance rank (percentile, matched precursors)')
        band.spines[['top', 'right', 'left']].set_visible(False)

    fig.suptitle(
        f'Quantification error vs abundance across tools  ({group_txt}: {", ".join(labs)})\n'
        f'(lower is better; bottom band = lowest-error tool per abundance bin)',
        fontweight='bold',
    )
    fig.tight_layout()

    if save_dir:
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)
    return m


# ── Pairwise tool comparison ──────────────────────────────────────────────────

def _abundance(df):
    """Return the mean of ``log_Intensity_mean_A`` and ``log_Intensity_mean_B`` per row."""
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
        save_figure(fig, os.path.join(save_dir, filename))
    plt.show() if show else plt.close(fig)
    return m
