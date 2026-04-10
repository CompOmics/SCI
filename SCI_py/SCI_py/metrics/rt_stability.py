import pandas as pd

def compute_rt_proportions(
    adata,
    rt_layer:     str,
    bins:         int = 20,
    normalize_rt: bool = False,
) -> pd.DataFrame:
    """
    Compute per-sample retention time proportion distributions.

    Parameters
    ----------
    adata:
        AnnData object. Retention times are taken from ``adata.layers[rt_layer]``,
        with shape (cells × features). ``adata.obs`` must contain
        ``source name`` and ``plate_idx``.
    rt_layer:
        Key in ``adata.layers`` containing retention time values.
    bins:
        Number of RT bins.
    normalize_rt:
        If ``True``, normalise RT per sample to ``[0, 1]``.

    Returns
    -------
    pd.DataFrame
        Long-format table with columns ``sample_id``, ``source name``,
        ``plate_idx``, ``bin_mid``, and ``proportion``.
    """
    import numpy as np

    X     = adata.layers[rt_layer]
    if hasattr(X, 'toarray'):
        X = X.toarray()

    obs   = adata.obs[['source name', 'plate_idx']].copy()

    rows = []
    for i, sample_id in enumerate(adata.obs_names):
        rt_vals = X[i, :]
        rt_vals = rt_vals[~np.isnan(rt_vals)]
        if len(rt_vals) == 0:
            continue

        if normalize_rt:
            rt_min, rt_max = rt_vals.min(), rt_vals.max()
            if rt_max > rt_min:
                rt_vals = (rt_vals - rt_min) / (rt_max - rt_min)

        counts, edges = np.histogram(rt_vals, bins=bins)
        bin_mids      = (edges[:-1] + edges[1:]) / 2
        proportions   = counts / counts.sum()

        source = obs.at[sample_id, 'source name']
        plate  = obs.at[sample_id, 'plate_idx']

        for mid, prop in zip(bin_mids, proportions):
            rows.append({
                'sample_id':   sample_id,
                'source name': source,
                'plate_idx':   plate,
                'bin_mid':     mid,
                'proportion':  prop,
            })

    return pd.DataFrame(rows)