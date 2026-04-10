import pandas as pd
import numpy as np

def compute_rt_norm_layer(
    adata,
    rt_layer:     str = 'EG.ApexRT'
) -> None:
    """
    Normalise retention times per sample to [0, 1] and store as a new layer.

    Normalisation is based on the earliest and latest eluting precursor
    within each sample. Samples where all RT values are identical are set to NaN.

    Parameters
    ----------
    adata:
        Precursor-level AnnData with ``rt_layer`` in ``adata.layers``.
    rt_layer:
        Source layer containing raw retention times.
    """
    X = adata.layers[rt_layer]
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.array(X, dtype=float)
    X[X == 0] = np.nan

    rt_min = np.nanmin(X, axis=1, keepdims=True)
    rt_max = np.nanmax(X, axis=1, keepdims=True)
    span   = rt_max - rt_min
    span[span == 0] = np.nan

    return (X - rt_min) / span
    
def compute_rt_proportions(
    adata,
    rt_layer: str,
    bins:     int = 20,
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
        Use ``'EG.ApexRT_norm'`` for sample-normalised RTs
        (see ``compute_rt_norm_layer``).
    bins:
        Number of RT bins.

    Returns
    -------
    pd.DataFrame
        Long-format table with columns ``sample_id``, ``source name``,
        ``plate_idx``, ``bin_mid``, and ``proportion``.
    """
    X = adata.layers[rt_layer]
    if hasattr(X, 'toarray'):
        X = X.toarray()

    obs  = adata.obs[['source name', 'plate_idx']].copy()
    rows = []

    for i, sample_id in enumerate(adata.obs_names):
        rt_vals = X[i, :]
        rt_vals = rt_vals[~np.isnan(rt_vals)]
        if len(rt_vals) == 0:
            continue

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