import numpy as np
import pandas as pd
import warnings

def compute_cv_matrix(adata, conditions, layer=''):
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