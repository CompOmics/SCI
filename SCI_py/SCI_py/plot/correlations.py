"""
Plots for the sample-correlation analysis.

Covers missingness masking (masking stats, identification heatmap, missingness
distributions), intra- and inter-lab correlations, instrument comparisons and
sample-pair scatters (Notebooks/3_Correlations/01_sample_correlations). The
computations live in ``SCI_py.metrics.correlations``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests
from itertools import combinations
from scipy.stats import linregress
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from ..constants.plotting import LAB_COLORS, LAB_ORDER
from .export import save_figure


def plot_masking_stats(
    stats_df: pd.DataFrame,
    save_dir: str | None = None,
    filename: str        = 'masking_stats.png',
    show:     bool       = True,
) -> None:
    """
    Bar chart of % values masked per lab, one subplot per condition.

    Parameters
    ----------
    stats_df:
        Output of ``create_reference_set``.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    conditions      = stats_df['condition'].unique()
    fig, axes = plt.subplots(1, len(conditions),
                              figsize=(6 * len(conditions), 5),
                              constrained_layout=True)
    axes = [axes] if len(conditions) == 1 else list(axes)

    for ax, cond in zip(axes, conditions):
        sub = stats_df[stats_df['condition'] == cond].sort_values('pct_masked', ascending=False)
        ax.barh(sub['lab'], sub['pct_masked'],
                color=[LAB_COLORS.get(l, '#2c7bb6') for l in sub['lab']], alpha=0.8)
        ax.axvline(sub['pct_masked'].mean(), color='grey', linestyle='--', linewidth=1,
                   label=f"mean: {sub['pct_masked'].mean():.1f}%")
        for _, row in sub.iterrows():
            ax.text(row['pct_masked'] + 0.2, row['lab'],
                    f"{int(row['total_masked']):,}", va='center', fontsize=8)
        ax.set(xlabel='% values masked', title=cond)
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=9)

    fig.suptitle('Masking applied to enforce equal missingness across labs', fontweight='bold')

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
    plt.show() if show else plt.close(fig)


def plot_masking_heatmap(
    masked_data: dict,
    condition:   str,
    save_dir:    str | None = None,
    filename:    str        = 'masking_heatmap.png',
    show:        bool       = True,
) -> None:
    """
    Binary identification heatmap after missingness equalisation.

    Parameters
    ----------
    masked_data:
        Output of ``create_reference_set``.
    condition:
        Condition to visualise.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    labs = [l for l in LAB_ORDER if l in masked_data.get(condition, {})]
    if not labs:
        print(f'No data for condition {condition}')
        return

    fig, axes = plt.subplots(1, len(labs), figsize=(4 * len(labs), 8),
                              sharey=True, constrained_layout=True)
    axes = [axes] if len(labs) == 1 else list(axes)
    cmap = ListedColormap(['#f0f0f0', '#2c7bb6'])

    for ax, lab in zip(axes, labs):
        X        = masked_data[condition][lab]
        row_ord  = np.argsort(np.sum(~np.isnan(X), axis=1))[::-1]
        col_ord  = np.argsort(np.sum(~np.isnan(X), axis=0))[::-1]
        X_binary = (~np.isnan(X[np.ix_(row_ord, col_ord)])).astype(float)

        ax.imshow(X_binary.T, aspect='auto', interpolation='none',
                  cmap=cmap, vmin=0, vmax=1)
        ax.set(title=lab, xlabel='Samples', xticks=[], yticks=[])
        ax.title.set_fontweight('bold')

    axes[0].set_ylabel('Proteins')
    fig.suptitle(f'Post-masking identification heatmap — {condition}', fontweight='bold')

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
        _save_binary_legend(save_dir)

    plt.show() if show else plt.close(fig)


def plot_missingness_distributions(
    masked_data:   dict,
    conditions:    list[str],
    level:         str        = 'protein',
    reference_lab: str | None = None,
    xlim_max:      float      = .5,
    save_dir:      str | None = None,
    filename:      str        = 'missingness_distributions.png',
    show:          bool       = True,
) -> pd.DataFrame:
    """
    Plot per-protein missingness distributions across labs after masking
    and test whether they differ using the Kolmogorov-Smirnov test.

    The KS test is run pairwise between all labs (or against a reference
    lab if specified).  Non-significant results validate that the masking
    step successfully equalised missingness across labs.

    Parameters
    ----------
    masked_data:
        Output of ``create_reference_set``.
    conditions:
        Conditions to visualise. Data is pooled across conditions.
    level:
        protein or sample.
    reference_lab:
        If provided, compare all labs against this lab only.
        If ``None``, all pairwise comparisons are made.
    xlim_max:
        Max value of the Per-protein missingness figure.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.

    Returns
    -------
    pd.DataFrame
        KS test results with columns ``lab_a``, ``lab_b``,
        ``ks_stat``, ``p_value``, ``p_value_bh``, ``sig``.
    """
    if level not in ('protein', 'sample'):
        raise ValueError("level must be 'protein' or 'sample'")

    axis = 0 if level == 'protein' else 1

    # Pool missingness across conditions per lab
    labs = [l for l in LAB_ORDER if any(l in masked_data.get(c, {}) for c in conditions)]

    missingness = {}
    for lab in labs:
        arrays = [
            1 - np.sum(~np.isnan(masked_data[c][lab]), axis=axis)
                / masked_data[c][lab].shape[axis]
            for c in conditions if lab in masked_data.get(c, {})
        ]
        missingness[lab] = np.concatenate(arrays)

    pairs = (
        [(reference_lab, lb) for lb in labs if lb != reference_lab]
        if reference_lab else list(combinations(labs, 2))
    )
    ks_records = []
    for la, lb in pairs:
        stat, p = ks_2samp(missingness[la], missingness[lb])
        ks_records.append({'lab_a': la, 'lab_b': lb, 'ks_stat': stat, 'p_value': p})

    ks_df = pd.DataFrame(ks_records)
    _, p_corrected, _, _ = multipletests(ks_df['p_value'], method='fdr_bh')
    ks_df['p_value_bh'] = p_corrected
    ks_df['sig']        = ks_df['p_value_bh'].apply(
        lambda p: '***' if p < 0.001 else '**' if p < 0.01
                  else '*' if p < 0.05 else 'ns'
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    ax = axes[0]
    for lab in labs:
        sns.kdeplot(missingness[lab], ax=ax, color=LAB_COLORS.get(lab, 'grey'),
                    linewidth=1.8, label=lab, fill=True, alpha=0.08)
    xlabel = ('Per-protein missingness fraction' if level == 'protein'
              else 'Per-sample missingness fraction')
    ax.set(xlabel=xlabel, ylabel='Density',
           title='Missingness distribution per lab',
           xlim=(0, xlim_max))
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    ks_matrix = pd.DataFrame(index=labs, columns=labs, dtype=float)
    for _, row in ks_df.iterrows():
        ks_matrix.loc[row['lab_a'], row['lab_b']] = row['ks_stat']
        ks_matrix.loc[row['lab_b'], row['lab_a']] = row['ks_stat']
    np.fill_diagonal(ks_matrix.values, 0.0)

    annot = ks_matrix.copy().astype(object)
    for _, row in ks_df.iterrows():
        label = f"{row['ks_stat']:.3f}\n({row['sig']})"
        annot.loc[row['lab_a'], row['lab_b']] = label
        annot.loc[row['lab_b'], row['lab_a']] = label
    np.fill_diagonal(annot.values, '0')

    sns.heatmap(
        ks_matrix.astype(float), ax=ax,
        annot=annot, fmt='',
        cmap='YlOrRd', vmin=0, vmax=0.5,
        linewidths=0.5, linecolor='white',
        cbar_kws={'label': 'KS statistic', 'shrink': 0.6},
    )
    ax.set(title='Pairwise KS test — FDR corrected (lower = more similar)')
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='y', rotation=0)

    cond_label = ', '.join(conditions)
    fig.suptitle(
        f'Missingness equalisation validation ({level}-level) — {cond_label}\n'
        f'ns = FDR-corrected p ≥ 0.05',
        fontweight='bold',
    )

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
    plt.show() if show else plt.close(fig)

    return ks_df


def plot_coobserved_diagnostic(
    correlation_results: dict,
    method:       str        = 'pearson',
    save_dir:     str | None = None,
    filename:     str        = 'coobserved_diagnostic.png',
    show:         bool       = True,
) -> None:
    """
    Diagnostic plot of co-observed protein counts vs correlation.

    Two panels per condition:

    1. KDE of co-observed protein counts per comparison type — validates
       that ``min_coobserved`` is appropriate and that intra/inter-lab
       pairs have similar co-observation depth.
    2. Scatter of co-observed count vs correlation with binned mean overlay
       — if correlation increases with n_coobserved, the threshold is too low.

    Parameters
    ----------
    correlation_results:
        Output of ``compute_pairwise_correlations``.
    method:
        Correlation method column to use (``'pearson'`` or ``'spearman'``).
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    n_coobs_col = f'n_coobserved{"_pearson" if method == "pearson" else ""}'
    conditions  = list(correlation_results.keys())

    _palette = {
        'intra_lab'                : '#2c7bb6',
        'inter_lab_same_instrument': '#abd9e9',
        'inter_lab_diff_instrument': '#d7191c',
    }

    fig, axes = plt.subplots(2, len(conditions),
                              figsize=(6 * len(conditions), 10),
                              constrained_layout=True)
    axes = axes.reshape(2, len(conditions))

    for col, cond in enumerate(conditions):
        df = correlation_results[cond].dropna(subset=[n_coobs_col, method])

        # --- Panel 1: co-observed distribution per comparison type ---
        ax = axes[0, col]
        for comp_type, color in _palette.items():
            sub = df[df['comparison_type'] == comp_type][n_coobs_col]
            if sub.empty:
                continue
            sns.kdeplot(sub, ax=ax, color=color, linewidth=2,
                        fill=True, alpha=0.2, label=comp_type)
            ax.axvline(sub.median(), color=color, linewidth=1.5, linestyle='--')

        ax.set(xlabel='# co-observed proteins', ylabel='Density',
               title=f'{cond} — co-observed distribution')
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8)

        # --- Panel 2: co-observed vs correlation ---
        ax = axes[1, col]
        for comp_type, color in _palette.items():
            sub = df[df['comparison_type'] == comp_type]
            if sub.empty:
                continue
            if len(sub) > 5000:
                sub = sub.sample(5000, random_state=42)

            ax.scatter(sub[n_coobs_col], sub[method],
                       color=color, alpha=0.1, s=5, label=comp_type)

            # Binned mean overlay
            sub         = sub.copy()
            sub['bin']  = pd.cut(sub[n_coobs_col], bins=20)
            mean_per_bin = sub.groupby('bin', observed=True)[method].mean()
            ax.plot([b.mid for b in mean_per_bin.index], mean_per_bin.values,
                    color=color, linewidth=2)

        ax.set(xlabel='# co-observed proteins',
               ylabel=f'{method.capitalize()} correlation',
               title=f'{cond} — correlation vs co-observed')
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8)

    fig.suptitle(f'Co-observed protein diagnostic ({method})', fontweight='bold')

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
    plt.show() if show else plt.close(fig)


def plot_intralab_correlation(
    correlation_results: dict,
    method:   str        = 'spearman',
    save_dir: str | None = None,
    filename_fmt: str    = '{cond}_{method}_intralab',
    ylim: tuple          = (0,1),
    show:     bool       = True,
) -> None:
    """
    Violin plot of intra-lab correlations per lab, one figure per condition.

    Returns
    -------
    dict mapping condition → pairwise Mann-Whitney statistics DataFrame.
    """
    all_stats = {}

    for cond, df in correlation_results.items():
        intra = df[df['same_lab']].dropna(subset=[method])
        lab_order = (intra.groupby('lab_i')[method]
                         .median()
                         .sort_values(ascending=False)
                         .index.tolist())

        fig, ax = plt.subplots(figsize=(len(lab_order) * 1.2 + 2, 5),
                                constrained_layout=True)
        sns.violinplot(data=intra, x='lab_i', y=method, order=lab_order,
                       hue='lab_i', hue_order=lab_order,
                       palette=LAB_COLORS, inner='quartile', legend=False, ax=ax)
        ax.set(xlabel='Lab', ylabel=f'{method.capitalize()} correlation',
               title=f'{cond} — intra-lab correlation', ylim=ylim)
        ax.title.set_fontweight('bold')
        ax.tick_params(axis='x', rotation=45)
        ax.spines[['top', 'right']].set_visible(False)

        if save_dir:
            save_figure(fig, f'{save_dir}/{filename_fmt.format(cond=cond, method=method)}')
        plt.show() if show else plt.close(fig)


def plot_interlab_correlation_heatmap(
    correlation_results: dict,
    method:       str        = 'spearman',
    vmin: float              = 0.0,
    vmax: float              = 1.0,
    save_dir:     str | None = None,
    filename_fmt: str        = '{cond}_{method}_interlab_heatmap.png',
    show:         bool       = True,
) -> None:
    """
    Clustered heatmap of median pairwise correlations between labs.
    Diagonal = intra-lab, off-diagonal = inter-lab.
    """
    for cond, df in correlation_results.items():
        df       = df.dropna(subset=[method])
        intra    = df[df['same_lab']]
        inter    = df[~df['same_lab']]
        all_labs = sorted(set(df['lab_i']) | set(df['lab_j']))

        matrix = pd.DataFrame(np.nan, index=all_labs, columns=all_labs)

        for lab in all_labs:
            sub = intra[intra['lab_i'] == lab][method]
            if not sub.empty:
                matrix.loc[lab, lab] = sub.median()

        for la, lb in combinations(all_labs, 2):
            pair = '__'.join(sorted([la, lb]))
            sub  = inter[inter['lab_pair'] == pair][method]
            if not sub.empty:
                matrix.loc[la, lb] = matrix.loc[lb, la] = sub.median()

        g = sns.clustermap(
            matrix.astype(float),
            cmap='RdYlGn', vmin=vmin, vmax=vmax,
            annot=True, fmt='.3f',
            figsize=(max(8, len(all_labs) * 0.9),) * 2,
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': f'Median {method} correlation', 'shrink': 0.5},
            annot_kws={'size': 8},
        )
        g.ax_row_dendrogram.set_visible(False)
        g.ax_row_dendrogram.set_position([0, 0, 0, 0])
        g.fig.suptitle(f'{cond} — median correlation heatmap ({method})',
                       fontweight='bold', y=1.02)

        if save_dir:
            # clustermap uses its own figure — save each format manually
            base = filename_fmt.format(cond=cond, method=method)
            base = os.path.splitext(base)[0]  # strip extension
            for fmt in ('png', 'svg', 'pdf'):
                out = os.path.join(save_dir, fmt, f'{base}.{fmt}')
                os.makedirs(os.path.dirname(out), exist_ok=True)
                g.fig.savefig(out, bbox_inches='tight', transparent=True, format=fmt)

        plt.show() if show else plt.close(g.fig)


def plot_instrument_correlation(
    correlation_results: dict,
    method:       str        = 'spearman',
    save_dir:     str | None = None,
    filename_fmt: str        = '{cond}_{method}_instrument.png',
    ylim: tuple              = (0,1),
    show:         bool       = True,
) -> None:
    """
    Violin plot comparing correlations by instrument type combination.
    """
    for cond, df in correlation_results.items():
        df = df.dropna(subset=[method]).copy()

        df['comparison'] = df.apply(
            lambda r: 'intra_lab' if r['same_lab']
            else f"inter_{r['instrument_i']}" if r['instrument_i'] == r['instrument_j']
            else f"inter_{r['instrument_i']}__{r['instrument_j']}",
            axis=1,
        )
        comp_order = ['intra_lab'] + sorted(
            c for c in df['comparison'].unique() if c != 'intra_lab'
        )
        palette = {'intra_lab': '#2c7bb6'}
        for color, key in zip(sns.color_palette('Set2', len(comp_order) - 1),
                               [c for c in comp_order if c != 'intra_lab']):
            palette[key] = color

        fig, ax = plt.subplots(figsize=(len(comp_order) * 1.5 + 2, 5),
                                constrained_layout=True)
        sns.violinplot(data=df, x='comparison', y=method, order=comp_order,
                       hue='comparison', hue_order=comp_order,
                       palette=palette, inner='quartile', legend=False, ax=ax)
        ax.set_xticks(range(len(comp_order)))
        ax.set_xticklabels(
            [c.replace('inter_', '').replace('__', ' vs ') for c in comp_order],
            rotation=30, ha='right',
        )
        ax.set(xlabel='', ylabel=f'{method.capitalize()} correlation',
               title=f'{cond} — correlation by instrument type', ylim=ylim)
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)

        if save_dir:
            save_figure(fig, f'{save_dir}/{filename_fmt.format(cond=cond, method=method)}')
        plt.show() if show else plt.close(fig)


def plot_sample_pair_correlations(
    correlation_results: dict,
    masked_data:  dict,
    condition:    str,
    lab_pair:     str | tuple,
    method:       str        = 'pearson',
    n_pairs:      int        = 6,
    corr_range:   tuple      = (None, None),
    random_seed:  int        = 42,
    save_dir:     str | None = None,
    filename:     str        = 'sample_pair_scatter.png',
    show:         bool       = True,
) -> None:
    """
    Scatter plots of log2 intensities for randomly sampled cell pairs.

    Parameters
    ----------
    correlation_results:
        Output of ``compute_pairwise_correlations``.
    masked_data:
        Output of ``create_reference_set``.
    condition:
        Condition to visualise.
    lab_pair:
        Lab pair as ``'LabA__LabB'`` string or ``('LabA', 'LabB')`` tuple.
    method:
        Correlation metric column to filter and annotate.
    n_pairs:
        Number of sample pairs to show.
    corr_range:
        ``(low, high)`` correlation bounds to filter pairs.
    random_seed:
        Seed for reproducible sampling.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    if isinstance(lab_pair, tuple):
        lab_pair = '__'.join(sorted(lab_pair))
    la, lb = lab_pair.split('__')

    df      = correlation_results[condition].dropna(subset=[method])
    pairs   = df[df['lab_pair'] == lab_pair].copy()
    low, high = corr_range
    if low  is not None: pairs = pairs[pairs[method] >= low]
    if high is not None: pairs = pairs[pairs[method] <= high]

    if pairs.empty:
        raise ValueError(f'No pairs found for {condition}, {lab_pair}, range={corr_range}')

    sampled  = pairs.sample(min(n_pairs, len(pairs)), random_state=random_seed)
    offsets  = _lab_offsets(masked_data, condition)
    X_a      = masked_data[condition][la]
    X_b      = masked_data[condition][lb]

    n_cols = min(3, n_pairs)
    n_rows = int(np.ceil(n_pairs / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 4, n_rows * 4),
                              constrained_layout=True)
    axes = np.array(axes).flatten()

    for ax_idx, (_, row) in enumerate(sampled.iterrows()):
        ax      = axes[ax_idx]
        i, j    = int(row['sample_i']), int(row['sample_j'])
        xi      = X_a[i - offsets[la], :]
        xj      = X_b[j - offsets[lb], :]
        obs     = ~np.isnan(xi) & ~np.isnan(xj)
        xi_obs  = xi[obs]
        xj_obs  = xj[obs]

        slope, intercept, *_ = linregress(xi_obs, xj_obs)
        x_line = np.linspace(xi_obs.min(), xi_obs.max(), 100)

        ax.scatter(xi_obs, xj_obs, alpha=0.3, s=8, color='#2c7bb6')
        ax.plot(x_line, intercept + slope * x_line, color='#d7191c', linewidth=1.5)
        ax.set(xlabel=f'{la} log2 intensity', ylabel=f'{lb} log2 intensity',
               title=f'r={row[method]:.3f}  n={int(row["n_coobserved"]):,}')
        ax.title.set_fontsize(9)
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)

    for ax in axes[len(sampled):]:
        ax.set_visible(False)

    fig.suptitle(f'{condition} — {la} vs {lb} ({method})', fontweight='bold')

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
    plt.show() if show else plt.close(fig)


def plot_intralab_pairwise_stats(
    stats_df:  pd.DataFrame,
    method:    str        = 'pearson',
    save_dir:  str | None = None,
    filename:  str        = 'intralab_pairwise_stats.png',
    show:      bool       = True,
) -> None:
    """
    Visualise pairwise Mann-Whitney statistics between lab intra-lab
    correlation distributions.

    Three panels:

    1. **Delta heatmap** — median correlation difference between each lab
       pair. Positive = row lab has higher intra-lab correlation than
       column lab. The dominant question: which labs are most reproducible?

    2. **CLES heatmap** — common language effect size. Values > 0.6
       indicate practically meaningful differences in reproducibility.

    3. **Delta distribution** — histogram of all pairwise deltas across
       conditions. Shows whether reproducibility differences are
       systematic or condition-specific.

    Parameters
    ----------
    stats_df:
        Output of the pairwise stats loop — must have columns
        ``group_a``, ``group_b``, ``delta``, ``cles``, ``sig``,
        ``condition``.
    method:
        Label for axis titles.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    conditions = stats_df['condition'].unique()
    labs       = sorted(set(stats_df['group_a']) | set(stats_df['group_b']))

    n_cond = len(conditions)
    fig, axes = plt.subplots(
        3, n_cond,
        figsize=(5 * n_cond, 14),
        constrained_layout=True,
    )
    axes = axes.reshape(3, n_cond)

    for col, cond in enumerate(conditions):
        sub = stats_df[stats_df['condition'] == cond]

        # --- Build symmetric matrices ---
        delta_matrix = pd.DataFrame(np.nan, index=labs, columns=labs)
        cles_matrix  = pd.DataFrame(0.5,   index=labs, columns=labs)

        for _, row in sub.iterrows():
            la, lb = row['group_a'], row['group_b']
            delta_matrix.loc[la, lb] =  row['delta']
            delta_matrix.loc[lb, la] = -row['delta']
            cles_matrix.loc[la, lb]  =  row['cles']
            cles_matrix.loc[lb, la]  =  1 - row['cles']

        np.fill_diagonal(delta_matrix.values, 0.0)
        np.fill_diagonal(cles_matrix.values,  0.5)

        # Annotation: delta value + significance
        delta_annot = delta_matrix.copy().astype(object)
        for _, row in sub.iterrows():
            la, lb = row['group_a'], row['group_b']
            label  = f"{row['delta']:+.3f}"
            delta_annot.loc[la, lb] =  label
            delta_annot.loc[lb, la] = f"{-row['delta']:+.3f}"
        np.fill_diagonal(delta_annot.values, '0')

        cles_annot = cles_matrix.copy().astype(object)
        for _, row in sub.iterrows():
            la, lb = row['group_a'], row['group_b']
            cles_annot.loc[la, lb] = f"{row['cles']:.2f}"
            cles_annot.loc[lb, la] = f"{1-row['cles']:.2f}"
        np.fill_diagonal(cles_annot.values, '0.50')

        # --- Panel 1: delta heatmap ---
        ax = axes[0, col]
        sns.heatmap(
            delta_matrix.astype(float), ax=ax,
            annot=delta_annot, fmt='',
            cmap='RdBu_r', center=0, vmin=-0.1, vmax=0.1,
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'Δ median correlation', 'shrink': 0.6},
            annot_kws={'size': 6, 'linespacing': .8}
        )
        ax.set(title=f'{cond} — Δ median {method} correlation')
        ax.title.set_fontweight('bold')
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)

        # --- Panel 2: CLES heatmap ---
        ax = axes[1, col]
        sns.heatmap(
            cles_matrix.astype(float), ax=ax,
            annot=cles_annot, fmt='',
            cmap='RdBu_r', center=0.5, vmin=0, vmax=1,
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'CLES', 'shrink': 0.6},
            annot_kws={'size': 6, 'linespacing': .8}
        )
        ax.set(title=f'{cond} — effect size (CLES)')
        ax.title.set_fontweight('bold')
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)

        # --- Panel 3: delta distribution ---
        ax = axes[2, col]
        sns.histplot(sub['delta'], ax=ax, bins=20,
                     color='#2c7bb6', alpha=0.7, kde=True)
        ax.axvline(0, color='grey', linewidth=1, linestyle='--')
        ax.axvline(sub['delta'].median(), color='#d7191c',
                   linewidth=1.5, linestyle='-',
                   label=f"median Δ = {sub['delta'].median():.3f}")
        ax.set(xlabel=f'Δ median {method} correlation',
               ylabel='Count',
               title=f'{cond} — distribution of pairwise Δ')
        ax.title.set_fontweight('bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8)

    fig.suptitle(
        f'Intra-lab reproducibility comparison — {method}\n'
        f'Δ = median(lab_a) − median(lab_b) | CLES > 0.6 = meaningful difference',
        fontweight='bold',
    )

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')
    plt.show() if show else plt.close(fig)


def _save_binary_legend(save_dir):
    """Save the Identified / Missing legend for the masking heatmap."""
    fig, ax = plt.subplots(figsize=(3, 1.5))
    ax.legend(
        handles=[Patch(facecolor='#2c7bb6', label='Identified'),
                 Patch(facecolor='#f0f0f0', edgecolor='grey', label='Missing / Masked')],
        loc='center', frameon=False, fontsize=9,
        title='Identification', title_fontsize=10,
    )
    ax.axis('off')
    save_figure(fig, f'{save_dir}/legend_masking_heatmap.png')
    plt.close(fig)
    
# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _lab_offsets(masked_data: dict, condition: str) -> dict[str, int]:
    """Compute cumulative sample offsets per lab in the concatenated matrix."""
    offsets, offset = {}, 0
    for lab in LAB_ORDER:
        if lab not in masked_data.get(condition, {}):
            continue
        offsets[lab] = offset
        offset += masked_data[condition][lab].shape[0]
    return offsets