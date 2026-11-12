"""Coefficients of variation (CV) per feature, lab and condition.

Plots for these results live in :mod:`SCI_py.plot.cvs`.
"""
import numpy as np
import pandas as pd
import warnings

from ..constants.labs import LABS

def compute_cv_matrix(adata, conditions, layer=''):
    """
    Compute the CV of every feature within each lab x condition group.

    CV = sample standard deviation (ddof=1) / mean x 100, ignoring NaNs. Only
    features with more than one non-NaN value in the group are kept. Values are
    used as stored, so pass a linear-scale (not log) layer for a standard CV.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx`` and ``source name`` in ``obs``.
    conditions:
        ``source name`` values to include.
    layer:
        Layer to use. Empty string uses ``X``.

    Returns
    -------
    pd.DataFrame
        Long format, one row per lab x condition x feature, with columns
        ``lab``, ``condition``, ``feature_id``, ``cv`` (in percent), ``n``
        (number of non-NaN values) and ``lab|condition``.
    """
    cvs_list = []
    samples = adata.obs[adata.obs['source name'].isin(conditions)].index

    for (plate, cond), group_samples in adata.obs.loc[samples].groupby(['plate_idx', 'source name']):
        if layer != '':
            X_sub = adata[group_samples.index, :].layers[layer]
        else:
            X_sub = adata[group_samples.index, :].X
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            col_means = np.nanmean(X_sub, axis=0)
            col_stds  = np.nanstd(X_sub, axis=0, ddof=1)
            n_valid   = np.sum(~np.isnan(X_sub), axis=0)

        mask = n_valid > 1
        cvs = (col_stds[mask] / col_means[mask]) * 100
        feature_ids = adata.var.index[mask]

        for feature_id, cv, n in zip(feature_ids, cvs, n_valid[mask]):
            cvs_list.append({
                'lab':           plate,
                'condition':     cond,
                'feature_id':    feature_id,
                'cv':            cv,
                'n':             n,
                'lab|condition': plate + '|' + cond
            })
    return pd.DataFrame(cvs_list)

def filter_common_proteome(cvs, condition_order):
    """Keep only features quantified in all labs, per condition."""
    n_labs = len(LABS)
    common = (
        cvs.groupby(['condition', 'feature_id'])['lab']
           .nunique()
           .reset_index(name='n_labs')
           .query('n_labs == @n_labs')
           .groupby('condition')['feature_id']
           .apply(set)
    )
    return cvs[cvs.apply(
        lambda row: row['feature_id'] in common.get(row['condition'], set()),
        axis=1
    )]
