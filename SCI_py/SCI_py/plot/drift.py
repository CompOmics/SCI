from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import linregress
from scipy.stats import spearmanr

from ..constants.plotting import CONDITION_ORDER, CONDITION_COLORS, LAB_COLORS, LAB_ORDER
from .utilities import _save, _save_legend

_DATE_COL = 'comment[Gradient time]_date'

def plot_features_by_day(
    adata,
    lab:        str,
    count_col:  str              = 'n_precursors',
    conditions: list[str]        = CONDITION_ORDER,
    normalize:  bool             = False,
    save_dir:   str | None       = None,
    filename:   str              = 'precursors_by_day.png',
    show:       bool             = True,
) -> None:
    """
    Plot per-condition feature counts by acquisition day for a single lab.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx``, ``source name``, ``count_col``, and
        ``comment[Gradient time]_date`` in ``obs``.
    lab:
        Lab identifier matching ``plate_idx``.
    count_col:
        Column in ``obs`` containing the count to plot.
    conditions:
        Conditions to include.
    normalize:
        If ``True``, normalise counts to % of max per condition.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    df = adata.obs[
        adata.obs['plate_idx'].eq(lab) &
        adata.obs['source name'].isin(conditions)
    ].copy()

    day_order = sorted(df[_DATE_COL].dropna().unique())
    if not day_order:
        print(f'No data for {lab}')
        return

    df['date_label'] = df[_DATE_COL].astype(str).str[:10]
    day_order_str    = [str(d)[:10] for d in day_order]

    if normalize:
        for cond in conditions:
            mask    = df['source name'] == cond
            max_val = df.loc[mask, count_col].max()
            if max_val > 0:
                df[count_col] = df[count_col].astype(float)
                df.loc[mask, count_col] = df.loc[mask, count_col] / max_val * 100

    ylabel = f'% of max {count_col}' if normalize else f'# {count_col}'

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    sns.boxplot(
        data      = df.sort_values(_DATE_COL),
        x         = 'source name', y=count_col,
        hue       = 'date_label', hue_order=day_order_str,
        order     = [c for c in CONDITION_ORDER if c in conditions],
        palette   = sns.color_palette('Blues', n_colors=len(day_order)),
        ax        = ax,
    )
    ax.set(xlabel='', ylabel=ylabel, title=lab)
    ax.title.set_fontweight('bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.get_legend().remove()

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')
        _save_legend(ax, save_dir, filename + '_legend')

    plt.show() if show else plt.close(fig)


def plot_feature_stability(
    adata,
    count_col:  str        = 'n_precursors',
    conditions: list[str]  = CONDITION_ORDER,
    save_dir:   str | None = None,
    filename:   str        = 'precursor_stability.png',
    show:       bool       = True,
) -> None:
    """
    Plot normalised feature counts over relative acquisition days per lab.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx``, ``source name``, ``count_col``, and
        ``comment[Gradient time]_date`` in ``obs``.
    count_col:
        Column in ``obs`` containing the count to plot.
    conditions:
        Conditions to include, each shown as a separate subplot.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    condition_order = [c for c in CONDITION_ORDER if c in conditions]
    combined        = _build_stability_df(adata, count_col, condition_order)

    if combined.empty:
        print('No data available.')
        return

    fig, axes = plt.subplots(
        1, len(condition_order),
        figsize=(6 * len(condition_order), 5),
        sharey=True, constrained_layout=True,
    )
    axes = [axes] if len(condition_order) == 1 else list(axes)

    for ax, cond in zip(axes, condition_order):
        cond_df = combined[combined['source name'] == cond]
        for lab in LAB_ORDER:
            lab_groups = cond_df[cond_df['plate_idx'] == lab].groupby('day_relative')[f'{count_col}_pct']
            if lab_groups.ngroups == 0:
                continue
            median = lab_groups.median()
            q25    = lab_groups.quantile(0.25)
            q75    = lab_groups.quantile(0.75)
            color  = LAB_COLORS.get(lab, 'black')
            ax.plot(median.index, median.values, color=color, linewidth=2, label=lab)
            ax.fill_between(median.index, q25.values, q75.values, color=color, alpha=0.15)

        ax.set(xlabel='Days since first acquisition',
               ylabel=f'% of max {count_col}' if ax == axes[0] else '',
               title=cond, xlim=(0, None))
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle(f'{count_col} stability over time', fontweight='bold')

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')
        _save_stability_legend(save_dir, filename)

    plt.show() if show else plt.close(fig)


def plot_drift(
    drift_df:           pd.DataFrame,
    lab:                str,
    count_col:          str,
    conditions:         list[str]  = CONDITION_ORDER,
    mode:               str        = 'conditions',
    line_per_condition: bool       = True,
    save_dir:           str | None = None,
    filename:           str        = 'drift.png',
    show:               bool       = True,
) -> None:
    """
    Plot intra-day drift for a single lab.

    Parameters
    ----------
    drift_df:
        Output of ``compute_drift``.
    lab:
        Lab identifier matching ``plate_idx``.
    count_col:
        Base count column name (without ``_pct``).
    conditions:
        Conditions to include.
    mode:
        ``'conditions'`` — one axis, one regression line per condition.
        ``'daily'``      — one subplot per condition, one line per day.
    line_per_condition:
        Only used when ``mode='conditions'``.  If ``False``, a single
        regression line is fitted across all conditions.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    pct_col         = f'{count_col}_pct'
    x_line          = np.linspace(0, 1, 100)
    condition_order = [c for c in CONDITION_ORDER if c in conditions]
    df              = drift_df[drift_df['plate_idx'] == lab].dropna(subset=[pct_col])

    if mode == 'conditions':
        fig, axes = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True)
        axes      = [axes]
        groups    = [(axes[0], condition_order)]
        palette   = CONDITION_COLORS
        group_col = 'source name'

    elif mode == 'daily':
        dates       = sorted(df['date'].unique())
        date_labels = [str(d)[:10] for d in dates]
        palette     = dict(zip(date_labels, sns.color_palette('viridis', n_colors=len(dates))))
        group_col   = 'date'
        fig, axes   = plt.subplots(1, len(condition_order),
                                    figsize=(7 * len(condition_order), 5),
                                    sharey=True, constrained_layout=True)
        axes  = [axes] if len(condition_order) == 1 else list(axes)
        groups = [(ax, [str(d)[:10] for d in sorted(
                      df[df['source name'] == cond]['date'].unique())])
                  for ax, cond in zip(axes, condition_order)]
    else:
        raise ValueError("mode must be 'conditions' or 'daily'")

    for ax, (cond, group_keys) in zip(axes, [(c, k) for c, (_, k) in
                                              zip(condition_order, groups)]):
        plot_df = df if mode == 'conditions' else df[df['source name'] == cond]

        if mode == 'conditions' and not line_per_condition:
            for c in condition_order:
                sub = plot_df[plot_df['source name'] == c]
                ax.scatter(sub['time_norm'], sub[pct_col],
                           color=CONDITION_COLORS.get(c, 'black'), alpha=0.3, s=15)
            _plot_regression(ax, plot_df['time_norm'].values,
                             plot_df[pct_col].values, 'black', 'All', x_line)
        else:
            for key in group_keys:
                sub   = (plot_df[plot_df['source name'] == key] if mode == 'conditions'
                         else plot_df[plot_df['date'].astype(str).str[:10] == key])
                color = palette.get(key, 'black')
                ax.scatter(sub['time_norm'], sub[pct_col], color=color, alpha=0.3, s=15)
                _plot_regression(ax, sub['time_norm'].values, sub[pct_col].values,
                                 color, key, x_line)

        title = f'{lab} — intra-day drift' if mode == 'conditions' else f'{lab} — {cond}'
        ax.set(xlabel='Normalised time within day [0–1]',
               ylabel=f'% of max {count_col}' if ax == axes[0] else '',
               title=title)
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8)

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')
        if mode == 'daily':
            _save_drift_daily_legend(palette, save_dir, filename)
        else:
            _save_legend(axes[0], save_dir, filename + '_legend')

    plt.show() if show else plt.close(fig)


def plot_drift_significance(
    stats_df:  pd.DataFrame,
    alpha:     float      = 0.05,
    save_dir:  str | None = None,
    filename:  str        = 'drift_significance.png',
    show:      bool       = True,
) -> None:
    """
    Plot a heatmap of drift regression p-values per lab and condition.

    Parameters
    ----------
    stats_df:
        Combined statistics DataFrame from ``fit_drift_stats`` with a
        ``p_value_bh`` column (BH-corrected) and ``adata_type`` column.
    alpha:
        Significance threshold shown on the colour scale.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    adata_types = stats_df['adata_type'].unique()
    fig, axes   = plt.subplots(1, len(adata_types),
                                figsize=(14, stats_df['lab'].nunique() * 0.5 + 2),
                                constrained_layout=True)
    axes = [axes] if len(adata_types) == 1 else list(axes)

    for ax, adata_type in zip(axes, adata_types):
        pivot = (stats_df[stats_df['adata_type'] == adata_type]
                 .pivot_table(index='lab', columns='condition',
                              values='p_value_bh', aggfunc='min'))
        annot = pivot.map(lambda p:
            f'{p:.3f} (***)' if p < 0.001 else
            f'{p:.3f} (**)' if p < 0.01  else
            f'{p:.3f} (*)'  if p < 0.05  else
            f'{p:.3f}'
        )
        sns.heatmap(
            pivot, ax=ax, annot=annot, fmt='',
            cmap='RdYlGn_r',
            norm=TwoSlopeNorm(vmin=0.00001, vcenter=alpha, vmax=0.5),
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'p-value (BH)', 'shrink': 0.6},
        )
        ax.set(title=f'{adata_type.capitalize()} drift significance',
               xlabel='', ylabel='Lab' if ax == axes[0] else '')
        ax.title.set_fontweight('bold')
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)

    fig.suptitle(f'Intra-day drift — regression significance (α={alpha})',
                 fontweight='bold')

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')

    plt.show() if show else plt.close(fig)


def plot_rt_drift_correlations(
    correlation_dict: dict,
    save_dir:         str | None = None,
    filename:         str        = 'rt_drift_correlations.png',
    show:             bool       = True,
):
    """
    Plot per-lab Spearman correlation distributions of RT drift.

    Parameters
    ----------
    correlation_dict:
        Dict mapping lab name → DataFrame with ``spearman_r`` column,
        as returned by ``compute_rt_drift_correlations``.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    n_labs = len(correlation_dict)
    ncols  = 2
    nrows  = (n_labs + 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.5 * nrows),
                             constrained_layout=True)
    axes = axes.flatten()

    for ax, (lab_name, correlations) in zip(axes, correlation_dict.items()):
        sns.boxplot(
            data       = correlations,
            x          = 'spearman_r',
            ax         = ax,
            color      = LAB_COLORS.get(lab_name, 'steelblue'),
            width      = 0.4,
            flierprops = dict(marker='o', markersize=2, alpha=0.4),
        )
        ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
        ax.set(title=lab_name, xlabel='Spearman ρ', ylabel='', xlim=(-1, 1))
        ax.spines[['top', 'right']].set_visible(False)

    for ax in axes[n_labs:]:
        ax.set_visible(False)

    fig.suptitle('RT drift correlations per lab', fontweight='bold')

    if save_dir:
        _save(fig, f'{save_dir}/{filename}')

    plt.show() if show else plt.close(fig)

def plot_rt_drift_precursor(
    adata,
    precursor:  str,
    rt_layer: str = 'EG.ApexRT',
    correlation_col = 'injection_order',
    color_by:   str       = 'source name',
    palette:    dict | None = None,
    title:      str | None  = None,
    save_dir:   str | None  = None,
    filename:   str | None  = None,
    show:       bool        = True,
):
    """
    Dotplot of apex RT vs injection order for a single precursor.

    Parameters
    ----------
    adata_precursor:
        Precursor-level AnnData.
    precursor:
        Precursor ID (must be in ``adata.var_names``).
    lab_name:
        Lab (plate_idx) to subset.
    conditions:
        Source names to include.
    color_by:
        obs column to colour dots by.
    palette:
        Dict mapping color_by values to colors. Falls back to color cycle if None.
    title:
        Figure title. Auto-generated if None.
    save_dir:
        Directory to save figure.
    filename:
        Output filename. Defaults to ``rt_drift_{precursor}.png``.
    show:
        Whether to call ``plt.show()``.
    """
    BUILTIN = {'source name': CONDITION_COLORS, 'plate_idx': LAB_COLORS}

    if precursor not in adata.var_names:
        raise ValueError(f"Precursor '{precursor}' not found in adata.var_names.")

    idx     = adata.var_names.get_loc(precursor)
    X       = adata.layers[f'{rt_layer}']
    rt_vals = (X[:, idx].toarray().flatten() if hasattr(X, 'toarray')
               else np.array(X[:, idx], dtype=float).flatten())
    rt_vals[rt_vals == 0] = np.nan

    obs             = adata.obs.copy()
    obs['_rt']      = rt_vals
    obs             = obs.dropna(subset=['_rt'])
    correlation_col_vals = obs[f'{correlation_col}'].values
    rt              = obs['_rt'].values

    _palette = (
        BUILTIN.get(color_by)
        or palette
        or dict(zip(
            obs[color_by].unique(),
            plt.rcParams['axes.prop_cycle'].by_key()['color'],
        ))
    )
    colors = obs[color_by].map(_palette).fillna('black')

    r, p = spearmanr(correlation_col_vals, rt)

    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.scatter(correlation_col_vals, rt, c=colors, s=18, alpha=0.7, linewidths=0)
    ax.set(
        xlabel = 'Injection order',
        ylabel = 'Apex RT (min)',
        title  = title or f'{precursor}  |  ρ = {r:.3f},  p = {p:.2e}  |  n = {len(rt)}',
    )
    ax.spines[['top', 'right']].set_visible(False)

    # legend
    for label, color in _palette.items():
        if label in obs[color_by].values:
            ax.scatter([], [], c=color, s=18, label=label)
    ax.legend(title=color_by, frameon=False, fontsize=8)

    if save_dir:
        fname = filename or f'rt_drift_{precursor}.png'
        _save(fig, f'{save_dir}/{fname}')

    plt.show() if show else plt.close(fig)
    
# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _plot_regression(ax, x, y, color, label, x_line):
    mask = ~np.isnan(y)
    x, y = x[mask], y[mask]
    if len(x) < 2:
        return
    slope, intercept, _, p, _ = linregress(x, y)
    ax.plot(x_line, intercept + slope * x_line, color=color, linewidth=2.5,
            label=f'{label} (slope={slope:.1f}%/inj, p={p:.3f})')


def _save_drift_daily_legend(palette, save_dir, filename):
    handles = [
        mlines.Line2D([], [], color=c, linewidth=2, marker='o', markersize=5, label=l)
        for l, c in palette.items()
    ]
    fig, ax = plt.subplots(figsize=(3, len(handles) * 0.35 + 0.5))
    ax.legend(handles=handles, title='Acquisition date', loc='center', frameon=False, fontsize=9)
    ax.axis('off')
    _save(fig, f'{save_dir}/{filename.replace(".png", "_legend.png")}')
    plt.close(fig)


def _build_stability_df(adata, count_col, condition_order) -> pd.DataFrame:
    dfs = []
    for lab in LAB_ORDER:
        df = adata.obs[
            adata.obs['plate_idx'].eq(lab) &
            adata.obs['source name'].isin(condition_order)
        ].copy().dropna(subset=[_DATE_COL])

        if df.empty:
            continue

        df[_DATE_COL]      = pd.to_datetime(df[_DATE_COL])
        df['day_relative'] = (df[_DATE_COL] - df[_DATE_COL].min()).dt.days.astype(int)

        for cond in condition_order:
            mask    = df['source name'] == cond
            max_val = df.loc[mask, count_col].max()
            if max_val > 0:
                df.loc[mask, f'{count_col}_pct'] = df.loc[mask, count_col] / max_val * 100

        dfs.append(df)

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def _save_stability_legend(save_dir, filename):
    handles = [
        mlines.Line2D([], [], color=LAB_COLORS.get(lab, 'black'), linewidth=2, label=lab)
        for lab in LAB_ORDER
    ]
    fig, ax = plt.subplots(figsize=(4, len(LAB_ORDER) * 0.4 + 0.5))
    ax.legend(handles=handles, title='Lab', loc='center', frameon=False, fontsize=9)
    ax.axis('off')
    _save(fig, f'{save_dir}/{filename + '_legend'}')
    plt.close(fig)