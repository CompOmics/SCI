"""Computations for the DIA-NN version and DIA-NN vs Spectronaut comparisons.

Works on dicts of AnnData objects (one per dataset/version). Plots for these
results live in :mod:`SCI_py.plot.software_versions`.

- :func:`shared_obs` -- samples present in every dataset
- :func:`compute_overlap_decomposition` -- per-sample shared / unique feature counts
- :func:`compute_binned_overlap` -- feature overlap as a function of intensity bin
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr


# ---------------------------------------------------------------------------
# Private data helpers (also used by SCI_py.plot.software_versions)
# ---------------------------------------------------------------------------


def _dense(X):
    """Return `X` as a dense array if it is sparse."""
    return X.toarray() if sp.issparse(X) else X


def _leading_acc(s):
    """Return the first ';'-separated accession of each entry in Series `s`."""
    return s.str.split(';').str[0].str.strip()


def _prec_feature_ids(adata):
    """Return precursor feature ids (the ``var`` index) as an array."""
    return np.asarray(adata.var.index)


def _prot_feature_ids(adata):
    """Return the leading ``PG.ProteinAccessions`` per feature, or the ``var`` index if absent."""
    if 'PG.ProteinAccessions' in adata.var.columns:
        return np.asarray(_leading_acc(adata.var['PG.ProteinAccessions']))
    return np.asarray(adata.var.index)


def _get_feat_fn(feature_level: str):
    """Return the feature-id function for 'precursor' or 'protein'."""
    if feature_level == 'precursor':
        return _prec_feature_ids
    elif feature_level == 'protein':
        return _prot_feature_ids
    raise ValueError(f"feature_level must be 'precursor' or 'protein', got {feature_level!r}")


def _count_ids(adata):
    """Return ``obs`` lab/condition with ``n_ids``: features > 0 per sample in ``X``."""
    X = _dense(adata.X)
    detected = np.nan_to_num(X, nan=0.0) > 0
    df = adata.obs[['plate_idx', 'source name']].copy()
    df['n_ids'] = detected.sum(axis=1)
    return df


def _detected_union(adata, feature_ids, samples):
    """Return the set of feature ids detected (> 0) in any of `samples`."""
    X = _dense(adata[samples].X)
    det = (np.nan_to_num(X, nan=0.0) > 0).any(axis=0)
    return set(np.asarray(feature_ids)[det])


def _per_sample_sets(adata, feature_ids, samples):
    """Return a dict of sample -> set of feature ids detected (> 0) in it."""
    sub = adata[samples]
    det = np.nan_to_num(_dense(sub.X), nan=0.0) > 0
    fid = np.asarray(feature_ids)
    return {o: set(fid[det[i]]) for i, o in enumerate(sub.obs.index)}


def _aligned_matrix(adata, samples, feature_ids_wanted, feature_id_array, layer='X_log2'):
    """Return the `layer` sub-matrix for `samples` x `feature_ids_wanted` (first column per id)."""
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
    """Return the row-wise Pearson r of `Xa` vs `Xb` (NaN if < `min_coobs` co-observed values)."""
    rs = np.full(Xa.shape[0], np.nan)
    for i in range(Xa.shape[0]):
        m = ~np.isnan(Xa[i]) & ~np.isnan(Xb[i])
        if m.sum() >= min_coobs:
            rs[i] = pearsonr(Xa[i, m], Xb[i, m])[0]
    return rs


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
