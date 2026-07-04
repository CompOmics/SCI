"""
Plotting and pre-computation functions for multi-dataset DIA-NN version comparisons
and cross-tool (Spectronaut) analyses.

Public API
----------
Data helpers
    shared_obs
    compute_overlap_decomposition
    compute_binned_overlap

Part A — DIA-NN version comparison
    plot_id_counts_overview
    plot_id_counts_per_lab
    plot_protein_counts
    plot_version_improvement
    plot_feature_overlap
    plot_quantitative_correlation

Part B — DIA-NN 2.6 vs Spectronaut
    plot_tool_counts_per_lab
    plot_tool_overlap
    plot_paired_per_lab
    plot_overlap_decomposition_per_lab
    plot_overlap_overview
    plot_binned_overlap_boxplot
    plot_binned_overlap_profile

Part C — Completeness and reproducibility
    plot_completeness_multi_dataset
    plot_cv_multi_dataset
"""
from __future__ import annotations

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr

from ..constants.plotting import LAB_COLORS
from .utilities import _save


# ---------------------------------------------------------------------------
# Private data helpers
# ---------------------------------------------------------------------------

def _dense(X):
    return X.toarray() if sp.issparse(X) else X


def _leading_acc(s):
    return s.str.split(';').str[0].str.strip()


def _prec_feature_ids(adata):
    return np.asarray(adata.var.index)


def _prot_feature_ids(adata):
    if 'PG.ProteinAccessions' in adata.var.columns:
        return np.asarray(_leading_acc(adata.var['PG.ProteinAccessions']))
    return np.asarray(adata.var.index)


def _get_feat_fn(feature_level: str):
    if feature_level == 'precursor':
        return _prec_feature_ids
    elif feature_level == 'protein':
        return _prot_feature_ids
    raise ValueError(f"feature_level must be 'precursor' or 'protein', got {feature_level!r}")


def _count_ids(adata):
    X = _dense(adata.X)
    detected = np.nan_to_num(X, nan=0.0) > 0
    df = adata.obs[['plate_idx', 'source name']].copy()
    df['n_ids'] = detected.sum(axis=1)
    return df


def _detected_union(adata, feature_ids, samples):
    X = _dense(adata[samples].X)
    det = (np.nan_to_num(X, nan=0.0) > 0).any(axis=0)
    return set(np.asarray(feature_ids)[det])


def _per_sample_sets(adata, feature_ids, samples):
    sub = adata[samples]
    det = np.nan_to_num(_dense(sub.X), nan=0.0) > 0
    fid = np.asarray(feature_ids)
    return {o: set(fid[det[i]]) for i, o in enumerate(sub.obs.index)}


def _aligned_matrix(adata, samples, feature_ids_wanted, feature_id_array, layer='X_log2'):
    id_to_col = {}
    for col, fid in enumerate(feature_id_array):
        if fid not in id_to_col:
            id_to_col[fid] = col
    cols = [id_to_col[f] for f in feature_ids_wanted]
    pos = {o: i for i, o in enumerate(adata.obs.index)}
    rows = [pos[s] for s in samples]
    X = _dense(adata.layers[layer])
    return X[np.ix_(rows, cols)]


def _per_sample_pearson(Xa, Xb, min_coobs=10):
    rs = np.full(Xa.shape[0], np.nan)
    for i in range(Xa.shape[0]):
        m = ~np.isnan(Xa[i]) & ~np.isnan(Xb[i])
        if m.sum() >= min_coobs:
            rs[i] = pearsonr(Xa[i, m], Xb[i, m])[0]
    return rs


def _grouped_box(ax, frames, conditions, color_map, ylabel, fliers=False):
    names = list(frames)
    n = len(names)
    width = 0.8 / max(n, 1)
    for ci, cond in enumerate(conditions):
        for gi, name in enumerate(names):
            vals = frames[name].loc[frames[name]['source name'] == cond, 'n_ids'].values
            if len(vals) == 0:
                continue
            pos = ci + (gi - (n - 1) / 2) * width
            bp = ax.boxplot(vals, positions=[pos], widths=width * 0.9,
                            patch_artist=True, showfliers=fliers,
                            medianprops=dict(color='black'))
            for box in bp['boxes']:
                box.set_facecolor(color_map[name])
                box.set_alpha(0.85)
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=20, ha='right')
    ax.set_ylabel(ylabel)
    ax.spines[['top', 'right']].set_visible(False)
    handles = [mpatches.Patch(facecolor=color_map[nm], label=nm) for nm in names]
    ax.legend(handles=handles, frameon=False, fontsize=8)


# ---------------------------------------------------------------------------
# Public data helpers
# ---------------------------------------------------------------------------

def shared_obs(adatas, labs=None, conditions=None):
    """
    Return observation indices present in all AnnData objects.

    Parameters
    ----------
    adatas:
        Iterable of AnnData objects to intersect.
    labs:
        If given, restrict each AnnData to these ``plate_idx`` values before intersecting.
    conditions:
        If given, restrict each AnnData to these ``source name`` values.

    Returns
    -------
    Sorted list of shared observation index values.
    """
    sets = []
    for a in adatas:
        o = a.obs
        if labs is not None:
            o = o[o['plate_idx'].isin(labs)]
        if conditions is not None:
            o = o[o['source name'].isin(conditions)]
        sets.append(set(o.index))
    return sorted(set.intersection(*sets)) if sets else []


def compute_overlap_decomposition(
    store: dict,
    name_a: str,
    name_b: str,
    lab_order: list[str],
    conditions: list[str],
    feature_level: str = 'precursor',
) -> pd.DataFrame:
    """
    Compute per-sample feature overlap between two datasets.

    For each matched sample, the identified features are split into three groups:
    found only in ``name_a``, shared, and found only in ``name_b``.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData.
    name_a:
        First dataset name.
    name_b:
        Second dataset name.
    lab_order:
        Labs to iterate over.
    conditions:
        Conditions (``source name`` values) to include.
    feature_level:
        ``'precursor'`` or ``'protein'``.

    Returns
    -------
    DataFrame with columns ``sample``, ``lab``, ``only_A``, ``shared``, ``only_B``,
    ``total``, ``jaccard``.
    """
    feat_fn = _get_feat_fn(feature_level)
    feat_ids_a = feat_fn(store[name_a])
    feat_ids_b = feat_fn(store[name_b])

    records = []
    for lab in lab_order:
        samples = shared_obs([store[name_a], store[name_b]], labs=[lab], conditions=conditions)
        sets_a = _per_sample_sets(store[name_a], feat_ids_a, samples)
        sets_b = _per_sample_sets(store[name_b], feat_ids_b, samples)
        for s in samples:
            sa, sb = sets_a[s], sets_b[s]
            inter = len(sa & sb)
            records.append({'sample': s, 'lab': lab,
                            'only_A': len(sa) - inter,
                            'shared': inter,
                            'only_B': len(sb) - inter})

    ov = pd.DataFrame(records)
    ov['total'] = ov['only_A'] + ov['shared'] + ov['only_B']
    ov['jaccard'] = ov['shared'] / ov['total']
    return ov


def compute_binned_overlap(
    adata_a,
    adata_b,
    name_a: str,
    name_b: str,
    labs: list[str],
    conditions: list[str],
    n_bins: int = 10,
    feature_level: str = 'precursor',
) -> pd.DataFrame:
    """
    Compute per-sample overlap as a function of abundance rank bin.

    Each tool's detected features are ranked by their own intensity (ascending) and
    divided into ``n_bins`` equal-frequency bins. Within each bin the fraction shared
    with the other tool is computed. Tool-specific ranks make bins comparable despite
    different intensity scales.

    Parameters
    ----------
    adata_a:
        AnnData for the first tool.
    adata_b:
        AnnData for the second tool.
    name_a:
        Label for the first tool.
    name_b:
        Label for the second tool.
    labs:
        Labs to include.
    conditions:
        Conditions (``source name`` values) to include.
    n_bins:
        Number of equal-frequency rank bins.
    feature_level:
        ``'precursor'`` or ``'protein'``.

    Returns
    -------
    DataFrame with columns ``sample``, ``lab``, ``tool``, ``bin``, ``bin_center``,
    ``n_total``, ``n_shared``, ``pct_shared``.
    """
    feat_fn = _get_feat_fn(feature_level)
    samples = shared_obs([adata_a, adata_b], labs=labs, conditions=conditions)
    feat_a = feat_fn(adata_a)
    feat_b = feat_fn(adata_b)

    Xa = _dense(adata_a[samples].X)
    Xb = _dense(adata_b[samples].X)

    bin_centers = [(b + 0.5) / n_bins * 100 for b in range(n_bins)]
    records = []
    for i, s in enumerate(samples):
        lab = adata_a.obs.loc[s, 'plate_idx']
        xa, xb = Xa[i], Xb[i]
        det_a = np.where(np.nan_to_num(xa, nan=0.0) > 0)[0]
        det_b = np.where(np.nan_to_num(xb, nan=0.0) > 0)[0]
        if len(det_a) < n_bins * 2 or len(det_b) < n_bins * 2:
            continue
        names_a, intens_a = feat_a[det_a], xa[det_a]
        names_b, intens_b = feat_b[det_b], xb[det_b]
        set_b, set_a = set(names_b), set(names_a)
        for tool, names, intens, partner in [
            (name_a, names_a, intens_a, set_b),
            (name_b, names_b, intens_b, set_a),
        ]:
            order = np.argsort(intens)
            n = len(order)
            bids = (np.arange(n) / n * n_bins).astype(int).clip(0, n_bins - 1)
            sn = names[order]
            in_partner = np.array([nm in partner for nm in sn], dtype=bool)
            for b in range(n_bins):
                mask = bids == b
                n_tot = int(mask.sum())
                n_sh = int(in_partner[mask].sum())
                records.append({'sample': s, 'lab': lab, 'tool': tool,
                                'bin': b, 'bin_center': bin_centers[b],
                                'n_total': n_tot, 'n_shared': n_sh,
                                'pct_shared': n_sh / n_tot * 100})

    binned = pd.DataFrame(records)
    n_kept = binned['sample'].nunique()
    print(f'{len(samples)} matched samples -> {n_kept} with sufficient features | '
          f'{len(binned):,} rows')
    print(binned.groupby('tool')['pct_shared'].describe().round(1).to_string())
    return binned


# ---------------------------------------------------------------------------
# Part A — DIA-NN version comparison
# ---------------------------------------------------------------------------

def plot_id_counts_overview(
    store: dict,
    dataset_order: list[str],
    dataset_colors: dict[str, str],
    conditions: list[str],
    labs: list[str],
    title: str = 'Precursor identifications per condition',
    save_dir: str | None = None,
    filename: str = 'A1_precursor_counts_overview.png',
    show: bool = True,
) -> None:
    """
    Grouped boxplot of precursor identification counts per condition across datasets.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData (precursor level).
    dataset_order:
        Dataset names to include, in display order.
    dataset_colors:
        Dict mapping dataset name → hex color.
    conditions:
        Conditions (``source name`` values) to display.
    labs:
        Labs (``plate_idx`` values) to restrict to.
    title:
        Figure title.
    save_dir:
        Directory for saved figures. ``None`` disables saving.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    frames = {name: _count_ids(store[name]) for name in dataset_order}
    frames = {k: v[v['plate_idx'].isin(labs)] for k, v in frames.items()}

    fig, ax = plt.subplots(figsize=(11, 5))
    _grouped_box(ax, frames, conditions, dataset_colors, '# precursors / sample')
    ax.set_title(title, fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_id_counts_per_lab(
    store: dict,
    dataset_order: list[str],
    dataset_colors: dict[str, str],
    lab_order: list[str],
    labs: list[str],
    sc_conditions: list[str],
    save_dir: str | None = None,
    filename: str = 'A1_precursor_counts_per_lab.png',
    show: bool = True,
) -> None:
    """
    Lab-by-lab faceted boxplot of precursor identification counts (single-cell conditions).

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData (precursor level).
    dataset_order:
        Dataset names to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    lab_order:
        Display order for labs.
    labs:
        Labs to restrict to.
    sc_conditions:
        Single-cell conditions to include.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    frames = {name: _count_ids(store[name]) for name in dataset_order}
    frames = {k: v[v['plate_idx'].isin(labs)] for k, v in frames.items()}

    fig, axes = plt.subplots(1, len(lab_order), figsize=(3.6 * len(lab_order), 5),
                             sharey=True, squeeze=False)
    axes = axes[0]
    for ax, lab in zip(axes, lab_order):
        for gi, name in enumerate(dataset_order):
            df = frames[name]
            vals = df[(df['plate_idx'] == lab) & (df['source name'].isin(sc_conditions))]['n_ids'].values
            bp = ax.boxplot(vals, positions=[gi], widths=0.6, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black'))
            for b in bp['boxes']:
                b.set_facecolor(dataset_colors[name])
                b.set_alpha(0.85)
        ax.set_xticks(range(len(dataset_order)))
        ax.set_xticklabels([n.replace('DIA-NN ', '') for n in dataset_order],
                           rotation=90, fontsize=7)
        ax.set_title(lab, fontweight='bold', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('# precursors / sample (single-cell)')
    fig.suptitle('Precursor IDs lab-by-lab — DIA-NN versions', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_protein_counts(
    store: dict,
    dataset_order: list[str],
    dataset_colors: dict[str, str],
    lab_order: list[str],
    labs: list[str],
    conditions: list[str],
    sc_conditions: list[str],
    save_dir: str | None = None,
    filename: str = 'A2_protein_counts.png',
    show: bool = True,
) -> None:
    """
    Two-panel figure: protein IDs per condition (left) and per lab (right).

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData (protein level).
    dataset_order:
        Dataset names to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    lab_order:
        Display order for labs.
    labs:
        Labs to restrict to.
    conditions:
        Conditions to show on the left panel.
    sc_conditions:
        Single-cell conditions for the per-lab panel.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    frames_p = {name: _count_ids(store[name]) for name in dataset_order}
    frames_p = {k: v[v['plate_idx'].isin(labs)] for k, v in frames_p.items()}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _grouped_box(axes[0], frames_p, conditions, dataset_colors, '# proteins / sample')
    axes[0].set_title('Protein IDs per condition', fontweight='bold')

    for gi, name in enumerate(dataset_order):
        df = frames_p[name]
        for li, lab in enumerate(lab_order):
            vals = df[(df['plate_idx'] == lab) & (df['source name'].isin(sc_conditions))]['n_ids'].values
            pos = li + (gi - (len(dataset_order) - 1) / 2) * (0.8 / len(dataset_order))
            bp = axes[1].boxplot(vals, positions=[pos], widths=0.35, patch_artist=True,
                                 showfliers=False, medianprops=dict(color='black'))
            for b in bp['boxes']:
                b.set_facecolor(dataset_colors[name])
                b.set_alpha(0.85)
    axes[1].set_xticks(range(len(lab_order)))
    axes[1].set_xticklabels(lab_order, rotation=30, ha='right')
    axes[1].set_ylabel('# proteins / sample (single-cell)')
    axes[1].set_title('Protein IDs lab-by-lab', fontweight='bold')
    axes[1].spines[['top', 'right']].set_visible(False)
    axes[1].legend(handles=[mpatches.Patch(facecolor=dataset_colors[n], label=n)
                            for n in dataset_order], frameon=False, fontsize=8)
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_version_improvement(
    version_pairs: list[tuple],
    prec: dict,
    prot: dict,
    labs: list[str],
    lab_order: list[str],
    sc_conditions: list[str],
    save_dir: str | None = None,
    filename: str = 'A3_version_improvement.png',
    show: bool = True,
) -> None:
    """
    Paired per-sample improvement scatter and per-lab gain boxplots (DIA-NN 2.2 → 2.6).

    Parameters
    ----------
    version_pairs:
        List of tuples ``(label, store_key, name_a, name_b)`` where ``store_key``
        is ``'prec'`` or ``'prot'``.
    prec:
        Dict mapping dataset name → AnnData (precursor level).
    prot:
        Dict mapping dataset name → AnnData (protein level).
    labs:
        Labs to restrict to.
    lab_order:
        Display order for labs.
    sc_conditions:
        Single-cell conditions.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    stores = {'prec': prec, 'prot': prot}
    fig, axes = plt.subplots(2, len(version_pairs),
                             figsize=(5 * len(version_pairs), 9))

    for col, (label, store_key, name_a, name_b) in enumerate(version_pairs):
        store = stores[store_key]
        samples = shared_obs([store[name_a], store[name_b]], labs=labs, conditions=sc_conditions)
        ca = _count_ids(store[name_a]).loc[samples]
        cb = _count_ids(store[name_b]).loc[samples]
        sample_labs = ca['plate_idx'].values

        ax = axes[0, col]
        for lab in lab_order:
            m = sample_labs == lab
            ax.scatter(ca['n_ids'].values[m], cb['n_ids'].values[m], s=14, alpha=0.6,
                       color=LAB_COLORS.get(lab, 'grey'), label=lab)
        hi = max(ca['n_ids'].max(), cb['n_ids'].max())
        ax.plot([0, hi], [0, hi], color='grey', linestyle='--', linewidth=1, zorder=0)
        ax.set_xlabel(name_a)
        ax.set_ylabel(name_b)
        ax.set_title(label, fontweight='bold', fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)
        if col == 0:
            ax.legend(frameon=False, fontsize=7)

        ax = axes[1, col]
        delta = cb['n_ids'].values - ca['n_ids'].values
        for li, lab in enumerate(lab_order):
            m = sample_labs == lab
            bp = ax.boxplot(delta[m], positions=[li], widths=0.6, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black'))
            for b in bp['boxes']:
                b.set_facecolor(LAB_COLORS.get(lab, 'grey'))
                b.set_alpha(0.8)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.set_xticks(range(len(lab_order)))
        ax.set_xticklabels(lab_order, rotation=30, ha='right')
        ax.set_ylabel('n(2.6) - n(2.2)')
        ax.set_title(f'{label} — gain per lab', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        print(f'{label}: median gain = {np.median(delta):+.0f}  (n={len(samples)} samples)')

    fig.suptitle('DIA-NN 2.2 → 2.6 paired per-sample improvement', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_feature_overlap(
    store: dict,
    pairs: list[tuple],
    level: str,
    labs: list[str],
    conditions: list[str],
    dataset_colors: dict[str, str],
    feature_level: str = 'precursor',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Pooled feature overlap bars and per-sample subset coverage histograms.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData.
    pairs:
        List of tuples ``(label, name_a, name_b)``.
    level:
        Display label for y-axis (e.g. ``'precursors'``).
    labs:
        Labs to restrict to.
    conditions:
        Conditions to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    feature_level:
        ``'precursor'`` or ``'protein'``.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'A4_{level}_overlap.png'

    feat_fn = _get_feat_fn(feature_level)
    fig, axes = plt.subplots(len(pairs), 2, figsize=(11, 4.2 * len(pairs)), squeeze=False)

    for row, (label, name_a, name_b) in enumerate(pairs):
        samples = shared_obs([store[name_a], store[name_b]], labs=labs, conditions=conditions)
        set_a = _detected_union(store[name_a], feat_fn(store[name_a]), samples)
        set_b = _detected_union(store[name_b], feat_fn(store[name_b]), samples)
        only_a = len(set_a - set_b)
        only_b = len(set_b - set_a)
        both = len(set_a & set_b)

        ax = axes[row, 0]
        bars = ax.bar(['2.2 only', 'both', '2.6 only'], [only_a, both, only_b],
                      color=[dataset_colors[name_a], '#888888', dataset_colors[name_b]])
        ax.bar_label(bars, fmt='%d', fontsize=8)
        ax.set_ylabel(f'# {level}')
        ax.set_title(f'{label} — pooled overlap', fontweight='bold', fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)

        sets_a = _per_sample_sets(store[name_a], feat_fn(store[name_a]), samples)
        sets_b = _per_sample_sets(store[name_b], feat_fn(store[name_b]), samples)
        cov = []
        for s in samples:
            sa, sb = sets_a[s], sets_b[s]
            if not sa or not sb:
                continue
            inter = len(sa & sb)
            cov.append(inter / min(len(sa), len(sb)))

        ax = axes[row, 1]
        ax.hist(cov, bins=30, color=dataset_colors[name_b], alpha=0.7)
        ax.axvline(np.median(cov), color='black', linestyle='--',
                   label=f'median={np.median(cov):.3f}')
        ax.set_xlabel('Coverage of smaller version per sample')
        ax.set_ylabel('# samples')
        ax.set_title(f'{label} — subset coverage', fontsize=10)
        ax.set_xlim(0, 1.02)
        ax.legend(frameon=False, fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_quantitative_correlation(
    store: dict,
    pairs: list[tuple],
    level: str,
    labs: list[str],
    conditions: list[str],
    dataset_colors: dict[str, str],
    feature_level: str = 'precursor',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Per-sample Pearson correlation histograms for log2 intensities on shared features.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData.
    pairs:
        List of tuples ``(label, name_a, name_b)``.
    level:
        Feature level label for the axis (e.g. ``'precursors'``).
    labs:
        Labs to restrict to.
    conditions:
        Conditions to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    feature_level:
        ``'precursor'`` or ``'protein'``.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'correlation_{level}.png'

    feat_fn = _get_feat_fn(feature_level)
    fig, axes = plt.subplots(1, len(pairs), figsize=(6 * len(pairs), 5), squeeze=False)
    axes = axes[0]

    for ax, (label, name_a, name_b) in zip(axes, pairs):
        samples = shared_obs([store[name_a], store[name_b]], labs=labs, conditions=conditions)
        feat = sorted(set(feat_fn(store[name_a])) & set(feat_fn(store[name_b])))
        Xa = _aligned_matrix(store[name_a], samples, feat, feat_fn(store[name_a]))
        Xb = _aligned_matrix(store[name_b], samples, feat, feat_fn(store[name_b]))
        rs = _per_sample_pearson(Xa, Xb)
        sample_labs = store[name_a].obs.loc[samples, 'plate_idx'].values

        for lab in labs:
            m = (sample_labs == lab) & ~np.isnan(rs)
            if m.sum() == 0:
                continue
            ax.hist(rs[m], bins=25, alpha=0.55, color=LAB_COLORS.get(lab, 'grey'), label=lab)
        ax.axvline(np.nanmedian(rs), color='black', linestyle='--',
                   label=f'median={np.nanmedian(rs):.3f}')
        ax.set_xlabel('Pearson r (log2)')
        ax.set_ylabel('# samples')
        ax.set_title(f'{label}\n{len(feat)} shared {level}', fontweight='bold', fontsize=10)
        ax.set_xlim(0, 1)
        ax.legend(frameon=False, fontsize=7)
        ax.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Part B — DIA-NN 2.6 vs Spectronaut
# ---------------------------------------------------------------------------

def plot_tool_counts_per_lab(
    store: dict,
    pair: list[str],
    dataset_colors: dict[str, str],
    lab_order: list[str],
    labs: list[str],
    sc_conditions: list[str],
    level: str = 'precursors',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Per-lab faceted boxplots comparing identification counts between two tools.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData.
    pair:
        Two dataset names to compare.
    dataset_colors:
        Dict mapping dataset name → hex color.
    lab_order:
        Display order for labs.
    labs:
        Labs to restrict to.
    sc_conditions:
        Single-cell conditions to include.
    level:
        Level label used in y-axis and title (e.g. ``'precursors'``).
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated from ``level`` if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'B1_{level}_counts_per_lab.png'

    frames_b = {name: _count_ids(store[name]) for name in pair}
    frames_b = {k: v[v['plate_idx'].isin(labs)] for k, v in frames_b.items()}
    short_names = [n.replace('DIA-NN 2.6 (2nd pass)', 'DIA-NN 2.6') for n in pair]

    fig, axes = plt.subplots(1, len(lab_order), figsize=(3.4 * len(lab_order), 4.5),
                             sharey=True, squeeze=False)
    axes = axes[0]
    for ax, lab in zip(axes, lab_order):
        for gi, name in enumerate(pair):
            df = frames_b[name]
            vals = df[(df['plate_idx'] == lab) & (df['source name'].isin(sc_conditions))]['n_ids'].values
            bp = ax.boxplot(vals, positions=[gi], widths=0.6, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black'))
            for b in bp['boxes']:
                b.set_facecolor(dataset_colors[name])
                b.set_alpha(0.85)
        ax.set_xticks(range(len(pair)))
        ax.set_xticklabels(short_names, rotation=90, fontsize=7)
        ax.set_title(lab, fontweight='bold', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel(f'# {level} / sample (single-cell)')
    fig.suptitle(f'{level.capitalize()} lab-by-lab — {pair[0]} vs {pair[1]}', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_tool_overlap(
    store: dict,
    name_a: str,
    name_b: str,
    level: str,
    labs: list[str],
    conditions: list[str],
    dataset_colors: dict[str, str],
    feature_level: str = 'precursor',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Pooled overlap bars and directional coverage histograms for two tools.

    Parameters
    ----------
    store:
        Dict mapping dataset name → AnnData.
    name_a:
        First dataset name.
    name_b:
        Second dataset name.
    level:
        Display label for y-axis (e.g. ``'precursors'``).
    labs:
        Labs to restrict to.
    conditions:
        Conditions to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    feature_level:
        ``'precursor'`` or ``'protein'``.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'B2_{level}_overlap.png'

    feat_fn = _get_feat_fn(feature_level)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    samples = shared_obs([store[name_a], store[name_b]], labs=labs, conditions=conditions)
    set_a = _detected_union(store[name_a], feat_fn(store[name_a]), samples)
    set_b = _detected_union(store[name_b], feat_fn(store[name_b]), samples)

    bars = axes[0].bar(
        [f'{name_a}\nonly', 'both', f'{name_b}\nonly'],
        [len(set_a - set_b), len(set_a & set_b), len(set_b - set_a)],
        color=[dataset_colors[name_a], '#888888', dataset_colors[name_b]])
    axes[0].bar_label(bars, fmt='%d', fontsize=8)
    axes[0].set_ylabel(f'# {level}')
    axes[0].set_title('Pooled overlap', fontweight='bold', fontsize=10)
    axes[0].spines[['top', 'right']].set_visible(False)

    sets_a = _per_sample_sets(store[name_a], feat_fn(store[name_a]), samples)
    sets_b = _per_sample_sets(store[name_b], feat_fn(store[name_b]), samples)
    cov_a, cov_b = [], []
    for s in samples:
        sa, sb = sets_a[s], sets_b[s]
        if not sa or not sb:
            continue
        inter = len(sa & sb)
        cov_a.append(inter / len(sa))
        cov_b.append(inter / len(sb))
    axes[1].hist(cov_a, bins=30, alpha=0.6, color=dataset_colors[name_a],
                 label=f'{name_a} in {name_b}')
    axes[1].hist(cov_b, bins=30, alpha=0.6, color=dataset_colors[name_b],
                 label=f'{name_b} in {name_a}')
    axes[1].set_xlabel('Per-sample coverage fraction')
    axes[1].set_ylabel('# samples')
    axes[1].set_title('Directional coverage', fontsize=10)
    axes[1].set_xlim(0, 1.02)
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].spines[['top', 'right']].set_visible(False)
    fig.suptitle(f'{level.capitalize()} — {name_a} vs {name_b}', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_paired_per_lab(
    levels: list[tuple],
    name_a: str,
    name_b: str,
    lab_order: list[str],
    sc_conditions: list[str],
    save_dir: str | None = None,
    filename: str = 'B4_paired_per_lab.png',
    show: bool = True,
) -> None:
    """
    Paired per-sample identification scatter per lab and per level.

    Parameters
    ----------
    levels:
        List of tuples ``(level_label, store)`` where ``store`` is a dict
        mapping dataset name → AnnData.
    name_a:
        First dataset name (x-axis).
    name_b:
        Second dataset name (y-axis).
    lab_order:
        Display order for labs.
    sc_conditions:
        Single-cell conditions to include.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    short = {name_a: name_a.split(' (')[0], name_b: name_b.split(' (')[0]}
    ncol = len(lab_order)
    nrow = len(levels)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 3.5 * nrow), squeeze=False)

    for r, (level, store) in enumerate(levels):
        for c, lab in enumerate(lab_order):
            ax = axes[r][c]
            samples = shared_obs([store[name_a], store[name_b]], labs=[lab],
                                 conditions=sc_conditions)
            if len(samples) == 0:
                ax.set_visible(False)
                continue
            na = _count_ids(store[name_a]).loc[samples, 'n_ids'].values
            nb = _count_ids(store[name_b]).loc[samples, 'n_ids'].values
            ax.scatter(na, nb, s=12, alpha=0.6, color=LAB_COLORS.get(lab, 'grey'))
            hi = max(na.max(), nb.max())
            ax.plot([0, hi], [0, hi], color='grey', linestyle='--', linewidth=1, zorder=0)
            med = np.median(nb - na)
            ax.set_title(f'{lab}\n{level}: Δ={med:+.0f}', fontsize=8, fontweight='bold')
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=6)
            if r == nrow - 1:
                ax.set_xlabel(f'# {short[name_a]}', fontsize=7)
            if c == 0:
                ax.set_ylabel(f'{level}\n# {short[name_b]}', fontsize=8)

    fig.suptitle(f'Paired per-sample identifications per lab — '
                 f'{short[name_a]} (x) vs {short[name_b]} (y)', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_overlap_decomposition_per_lab(
    ov: pd.DataFrame,
    name_a: str,
    name_b: str,
    lab_order: list[str],
    dataset_colors: dict[str, str],
    save_dir: str | None = None,
    filename: str = 'B5_overlap_per_sample.png',
    show: bool = True,
) -> None:
    """
    Per-sample stacked feature overlap bars, faceted by lab.

    Parameters
    ----------
    ov:
        DataFrame with columns ``lab``, ``only_A``, ``shared``, ``only_B``, ``total``,
        as returned by ``compute_overlap_decomposition``.
    name_a:
        Label for tool A (used in legend).
    name_b:
        Label for tool B (used in legend).
    lab_order:
        Display order for labs.
    dataset_colors:
        Dict mapping dataset name → hex color.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    cA = dataset_colors[name_a]
    cB = dataset_colors[name_b]

    fig, axes = plt.subplots(1, len(lab_order), figsize=(3.6 * len(lab_order), 4.5),
                             sharey=True, squeeze=False)
    axes = axes[0]
    for ax, lab in zip(axes, lab_order):
        sub = ov[ov['lab'] == lab].sort_values('total', ascending=False).reset_index(drop=True)
        x = np.arange(len(sub))
        ax.bar(x, sub['only_A'], width=1.0, color=cA, label=f'{name_a} only')
        ax.bar(x, sub['shared'], width=1.0, bottom=sub['only_A'],
               color='#999999', label='shared')
        ax.bar(x, sub['only_B'], width=1.0, bottom=sub['only_A'] + sub['shared'],
               color=cB, label=f'{name_b} only')
        ax.set_title(lab, fontweight='bold', fontsize=9)
        ax.set_xlabel('samples (sorted by total)')
        ax.set_xticks([])
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('# precursors')
    axes[-1].legend(frameon=False, fontsize=7, loc='upper right')
    fig.suptitle(f'Per-sample precursor overlap — {name_a} vs {name_b} (per lab)',
                 fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_overlap_overview(
    ov: pd.DataFrame,
    name_a: str,
    name_b: str,
    lab_order: list[str],
    dataset_colors: dict[str, str],
    save_dir: str | None = None,
    filename: str = 'B5_overlap_overview.png',
    show: bool = True,
) -> None:
    """
    Mean per-sample overlap across all samples and per lab (counts and fractions of union).

    Parameters
    ----------
    ov:
        DataFrame with columns ``lab``, ``only_A``, ``shared``, ``only_B``,
        as returned by ``compute_overlap_decomposition``.
    name_a:
        Label for tool A.
    name_b:
        Label for tool B.
    lab_order:
        Display order for labs.
    dataset_colors:
        Dict mapping dataset name → hex color.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    cA = dataset_colors[name_a]
    cB = dataset_colors[name_b]
    groups = ['All'] + lab_order

    def _grp(g):
        return ov if g == 'All' else ov[ov['lab'] == g]

    mean_counts = {g: _grp(g)[['only_A', 'shared', 'only_B']].mean() for g in groups}
    x = np.arange(len(groups))
    oa = np.array([mean_counts[g]['only_A'] for g in groups])
    sh = np.array([mean_counts[g]['shared'] for g in groups])
    ob = np.array([mean_counts[g]['only_B'] for g in groups])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].bar(x, oa, color=cA, label=f'{name_a} only')
    axes[0].bar(x, sh, bottom=oa, color='#999999', label='shared')
    axes[0].bar(x, ob, bottom=oa + sh, color=cB, label=f'{name_b} only')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(groups, rotation=30, ha='right')
    axes[0].set_ylabel('Mean # precursors / sample')
    axes[0].set_title('Mean overlap (counts)', fontweight='bold')
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].spines[['top', 'right']].set_visible(False)

    tot = oa + sh + ob
    axes[1].bar(x, oa / tot * 100, color=cA, label=f'{name_a} only')
    axes[1].bar(x, sh / tot * 100, bottom=oa / tot * 100, color='#999999', label='shared')
    axes[1].bar(x, ob / tot * 100, bottom=(oa + sh) / tot * 100, color=cB,
                label=f'{name_b} only')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(groups, rotation=30, ha='right')
    axes[1].set_ylabel('% of union')
    axes[1].set_title('Mean overlap (fraction of union)', fontweight='bold')
    axes[1].set_ylim(0, 100)
    axes[1].spines[['top', 'right']].set_visible(False)

    fig.suptitle(f'Precursor overlap overview — {name_a} vs {name_b}', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_binned_overlap_boxplot(
    binned_df: pd.DataFrame,
    name_a: str,
    name_b: str,
    dataset_colors: dict[str, str],
    n_bins: int = 10,
    level: str = 'precursors',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Boxplot of overlap percentage per abundance rank bin for each tool.

    Parameters
    ----------
    binned_df:
        DataFrame from ``compute_binned_overlap`` with columns ``tool``, ``bin``,
        ``pct_shared``, ``sample``.
    name_a:
        First dataset name.
    name_b:
        Second dataset name.
    dataset_colors:
        Dict mapping dataset name → hex color.
    n_bins:
        Number of abundance rank bins used when computing ``binned_df``.
    level:
        Feature level label for titles (e.g. ``'precursors'``).
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'B6_binned_overlap_boxplot_{level}.png'

    short = {name_a: name_a.split(' (')[0], name_b: name_b.split(' (')[0]}
    bin_positions = [(b + 0.5) / n_bins * 100 for b in range(n_bins)]
    bin_labels = [f'{b * 100 // n_bins}–{(b + 1) * 100 // n_bins}%' for b in range(n_bins)]
    box_width = 100 / n_bins * 0.65

    fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=True)

    for ax, tool, partner in zip(axes, [name_a, name_b], [name_b, name_a]):
        col = dataset_colors[tool]
        sub = binned_df[binned_df['tool'] == tool]
        groups = [sub[sub['bin'] == b]['pct_shared'].dropna().values for b in range(n_bins)]
        bp = ax.boxplot(groups, positions=bin_positions, widths=box_width,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=1.5))
        for box in bp['boxes']:
            box.set_facecolor(col)
            box.set_alpha(0.75)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 105)
        ax.set_xticks(bin_positions)
        ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('abundance rank bin (percentile within tool, 0 = lowest)')
        ax.set_ylabel('% shared with other tool')
        ax.set_title(
            f'{short[tool]}\n% of its {level} per rank bin also detected by '
            f'{short[partner]}  (n={sub["sample"].nunique()} samples)',
            fontweight='bold', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle(
        f'{level.capitalize()} overlap per abundance rank bin — {name_a} vs {name_b}\n'
        'Each box = distribution across all matched samples. '
        'Bins are equal-frequency within each tool.',
        fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_binned_overlap_profile(
    binned_df: pd.DataFrame,
    name_a: str,
    name_b: str,
    dataset_colors: dict[str, str],
    example_sample: str,
    n_bins: int = 10,
    level: str = 'precursors',
    save_dir: str | None = None,
    filename: str | None = None,
    show: bool = True,
) -> None:
    """
    Per-rank profile: single example sample (stacked bar) and all-sample mean ± SD.

    Parameters
    ----------
    binned_df:
        DataFrame from ``compute_binned_overlap``.
    name_a:
        First dataset name.
    name_b:
        Second dataset name.
    dataset_colors:
        Dict mapping dataset name → hex color.
    example_sample:
        Observation index to use for the example panel.
    n_bins:
        Number of abundance rank bins used when computing ``binned_df``.
    level:
        Feature level label for titles.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename. Auto-generated if ``None``.
    show:
        Whether to call ``plt.show()``.
    """
    if filename is None:
        filename = f'B6_binned_overlap_profile_{level}.png'

    short = {name_a: name_a.split(' (')[0], name_b: name_b.split(' (')[0]}
    bin_x = [(b + 0.5) / n_bins * 100 for b in range(n_bins)]
    bwidth = 100 / n_bins * 0.8

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey='row')

    for row, (tool, partner) in enumerate([(name_a, name_b), (name_b, name_a)]):
        col = dataset_colors[tool]
        sub = binned_df[binned_df['tool'] == tool]

        ax = axes[row, 0]
        ex = sub[sub['sample'] == example_sample].sort_values('bin')
        if not ex.empty:
            pct_sh = ex['pct_shared'].values
            ax.bar(bin_x, pct_sh, width=bwidth, color=col, alpha=0.85,
                   label=f'shared with {short[partner]}')
            ax.bar(bin_x, 100 - pct_sh, width=bwidth, color='#cccccc', alpha=0.85,
                   bottom=pct_sh, label=f'unique to {short[tool]}')
        else:
            ax.text(0.5, 0.5, 'example sample not found', ha='center', va='center',
                    transform=ax.transAxes, color='grey')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 110)
        sample_label = example_sample.split('|')[0] if '|' in example_sample else example_sample
        ax.set_title(f'{short[tool]} — example sample\n{sample_label}',
                     fontsize=9, fontweight='bold')
        ax.set_ylabel(f'% of {level}')
        if row == 1:
            ax.set_xlabel('abundance rank (percentile within tool, 0 = lowest)')
        ax.legend(frameon=False, fontsize=8, loc='lower right')
        ax.spines[['top', 'right']].set_visible(False)

        ax = axes[row, 1]
        for s in sub['sample'].unique():
            ss = sub[sub['sample'] == s].sort_values('bin')
            if len(ss) == n_bins:
                ax.plot(bin_x, ss['pct_shared'].values,
                        color='grey', alpha=0.04, linewidth=0.8)
        mean_v = sub.groupby('bin')['pct_shared'].mean().sort_index().values
        std_v = sub.groupby('bin')['pct_shared'].std().sort_index().values
        ax.plot(bin_x, mean_v, color=col, linewidth=2.5, label='mean')
        ax.fill_between(bin_x,
                        np.clip(mean_v - std_v, 0, 100),
                        np.clip(mean_v + std_v, 0, 100),
                        color=col, alpha=0.3, label='±1 SD')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 110)
        ax.set_title(f'{short[tool]} — all samples\n'
                     f'% shared with {short[partner]}  (n={sub["sample"].nunique()})',
                     fontsize=9, fontweight='bold')
        if row == 1:
            ax.set_xlabel('abundance rank (percentile within tool, 0 = lowest)')
        ax.legend(frameon=False, fontsize=8, loc='lower right')
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylabel('')

    fig.suptitle(
        f'% {level} shared with other tool as a function of abundance rank\n'
        f'{name_a} vs {name_b} — top: {short[name_a]} perspective, '
        f'bottom: {short[name_b]} perspective',
        fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Part C — Completeness and reproducibility
# ---------------------------------------------------------------------------

def plot_completeness_multi_dataset(
    comp: dict,
    dataset_order: list[str],
    dataset_colors: dict[str, str],
    lab_order: list[str],
    save_dir: str | None = None,
    filename: str = 'C1_completeness_per_lab.png',
    show: bool = True,
) -> None:
    """
    Completeness curves per lab, one line per dataset.

    Parameters
    ----------
    comp:
        Dict mapping dataset name → ``(thresholds, retention_curves, _)`` as returned
        by ``compute_completeness``.
    dataset_order:
        Dataset names to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    lab_order:
        Display order for labs.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    fig, axes = plt.subplots(1, len(lab_order), figsize=(3.6 * len(lab_order), 4.5),
                             sharey=True, squeeze=False)
    axes = axes[0]
    for ax, lab in zip(axes, lab_order):
        for name in dataset_order:
            thr, ret, _ = comp[name]
            if lab in ret:
                ax.plot(thr, ret[lab], color=dataset_colors[name], linewidth=1.6,
                        label=name, alpha=0.9)
        ax.axhline(50, color='grey', linewidth=0.7, linestyle=':')
        ax.set(xlim=(0, 100), ylim=(0, 105))
        ax.set_xlabel('Min. % of samples')
        ax.set_title(lab, fontweight='bold', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('% precursors retained')
    axes[-1].legend(frameon=False, fontsize=7, loc='upper right')
    fig.suptitle('Precursor completeness lab-by-lab (single-cell)', fontweight='bold')
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_cv_multi_dataset(
    cv_dfs: dict,
    dataset_order: list[str],
    dataset_colors: dict[str, str],
    lab_order: list[str],
    labs: list[str],
    save_dir: str | None = None,
    filename: str = 'C2_cv_per_lab.png',
    show: bool = True,
) -> None:
    """
    Grouped bar chart of median precursor CV per lab across datasets.

    Parameters
    ----------
    cv_dfs:
        Dict mapping dataset name → DataFrame with columns ``lab`` and ``cv``,
        as returned by ``compute_cv_matrix``.
    dataset_order:
        Dataset names to include.
    dataset_colors:
        Dict mapping dataset name → hex color.
    lab_order:
        Display order for labs.
    labs:
        Labs to include (typically the set available across all datasets).
    save_dir:
        Directory for saved figures.
    filename:
        Output filename base.
    show:
        Whether to call ``plt.show()``.
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    n = len(dataset_order)
    width = 0.8 / n
    for gi, name in enumerate(dataset_order):
        df = cv_dfs[name]
        meds = []
        for lab in lab_order:
            sub = df[(df['lab'] == lab) & np.isfinite(df['cv']) & (df['cv'] > 0)]
            meds.append(sub['cv'].median() if len(sub) else np.nan)
        positions = np.arange(len(lab_order)) + (gi - (n - 1) / 2) * width
        ax.bar(positions, meds, width=width * 0.9, color=dataset_colors[name],
               label=name, alpha=0.9)
    ax.set_xticks(range(len(lab_order)))
    ax.set_xticklabels(lab_order, rotation=30, ha='right')
    ax.set_ylabel('Median precursor CV (%)')
    ax.set_title(
        'Median CV per lab — DIA-NN versions comparable; Spectronaut on different scale',
        fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))
    if show:
        plt.show()
    else:
        plt.close(fig)
