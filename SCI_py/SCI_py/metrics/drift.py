from __future__ import annotations

import pandas as pd

from scipy.stats import linregress
from tqdm import tqdm
import numpy as np
from scipy.stats import spearmanr
from ..constants.maps import LABS

_DATETIME_COL = 'comment[Gradient time]_datetime'

def compute_time_norm(
    adata,
    inplace: bool = False,
):
    """
    Normalise acquisition datetime to [0, 1] within each lab × day group.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx`` and ``comment[Gradient time]_datetime`` in ``obs``.
    inplace:
        If ``True``, add ``time_norm`` directly to ``adata.obs``.
        If ``False``, return a Series indexed by ``adata.obs_names``.

    Returns
    -------
    pd.Series | None
        Series of normalised times if ``inplace=False``, else ``None``.
    """
    df = adata.obs[[_DATETIME_COL, 'plate_idx']].copy()
    df[_DATETIME_COL] = pd.to_datetime(df[_DATETIME_COL])
    df['date']        = df[_DATETIME_COL].dt.floor('D')

    grp   = df.groupby(['plate_idx', 'date'])[_DATETIME_COL]
    t_min = grp.transform('min')
    t_max = grp.transform('max')
    span  = (t_max - t_min).dt.total_seconds()

    time_norm = (df[_DATETIME_COL] - t_min).dt.total_seconds() / span
    time_norm = time_norm.where(span > 0, 0.5)
    time_norm.name = 'time_norm'

    if inplace:
        adata.obs['time_norm'] = time_norm
    else:
        return time_norm
    
def compute_drift(
    adata,
    conditions: list[str],
    count_col:  str,
    labs:       list[str] | None = None,
) -> pd.DataFrame:
    """
    Compute normalised intra-day drift data for all labs.

    Parameters
    ----------
    adata:
        AnnData with ``plate_idx``, ``source name``, ``count_col``,
        and ``comment[Gradient time]_datetime`` in ``obs``.
    conditions:
        Source names to include.
    count_col:
        Column in ``obs`` containing the count to normalise.
    labs:
        Labs to include.  ``None`` uses ``LAB_ORDER``.

    Returns
    -------
    pd.DataFrame
        Long-format table with normalised time and count columns,
        ready for ``plot_drift`` and ``plot_drift_daily``.
    """
    labs = labs or LABS
    dfs  = []

    for lab in labs:
        df = adata.obs[
            adata.obs['plate_idx'].eq(lab) &
            adata.obs['source name'].isin(conditions)
        ].copy().dropna(subset=[_DATETIME_COL])

        if df.empty:
            continue
        
        df[_DATETIME_COL] = pd.to_datetime(df[_DATETIME_COL])
        df['date'] = df[_DATETIME_COL].dt.floor('D')

        # Normalise time within each day to [0, 1]
        for date in df['date'].unique():
            mask        = df['date'] == date
            t_min, t_max = df.loc[mask, _DATETIME_COL].min(), df.loc[mask, _DATETIME_COL].max()
            df.loc[mask, 'time_norm'] = (
                (df.loc[mask, _DATETIME_COL] - t_min) / (t_max - t_min)
                if t_max > t_min else 0.5
            )

        # Normalise counts to % of max per condition per day
        for cond in conditions:
            for date in df['date'].unique():
                mask    = df['source name'].eq(cond) & df['date'].eq(date)
                max_val = df.loc[mask, count_col].max()
                if max_val > 0:
                    df.loc[mask, f'{count_col}_pct'] = df.loc[mask, count_col] / max_val * 100

        df['plate_idx'] = lab
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

def compute_rt_drift_correlations(
    adata,
    correlation_col: str = 'injection_order',
    rt_layer:        str = 'EG.ApexRT',
    min_n:           int = 10,
) -> pd.DataFrame:
    """
    Compute per-precursor Spearman correlation between apex RT and a given column.

    Parameters
    ----------
    adata:
        Precursor-level AnnData with ``rt_layer`` in layers and
        ``correlation_col`` in obs.
    correlation_col:
        Column in ``obs`` to correlate RT against.
    rt_layer:
        Layer to use for retention times. Use ``'EG.ApexRT_norm'`` for
        sample-normalised RTs (see ``compute_rt_norm_layer``).
    min_n:
        Minimum number of non-NaN identifications required per precursor.

    Returns
    -------
    pd.DataFrame
        One row per precursor with columns:
        ``precursor``, ``n_identified``, ``spearman_r``, ``spearman_p``,
        ``rt_min``, ``rt_max``, ``rt_range``.
    """
    X = adata.layers[rt_layer]
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.array(X, dtype=float)
    X[X == 0] = np.nan

    correlation_col_vals = adata.obs[correlation_col].values
    precursor_ids        = adata.var_names

    rows = []
    for i, precursor in tqdm(enumerate(precursor_ids), desc='Calculating correlations...'):
        rt_vals = X[:, i]
        valid   = ~np.isnan(rt_vals)
        n       = valid.sum()

        if n < min_n:
            continue

        r, p = spearmanr(correlation_col_vals[valid], rt_vals[valid])
        rows.append({
            'precursor':    precursor,
            'n_identified': n,
            'spearman_r':   r,
            'spearman_p':   p,
            'rt_min':       rt_vals[valid].min(),
            'rt_max':       rt_vals[valid].max(),
            'rt_range':     rt_vals[valid].max() - rt_vals[valid].min(),
        })

    return pd.DataFrame(rows).sort_values('spearman_r', ascending=False).reset_index(drop=True)


def fit_drift_stats(
    drift_df:  pd.DataFrame,
    count_col: str,
    mode:      str = 'conditions',
) -> pd.DataFrame:
    """
    Fit linear regressions on drift data and return statistics.

    Parameters
    ----------
    drift_df:
        Output of ``compute_drift``.
    count_col:
        Base count column name (without ``_pct``).
    mode:
        ``'conditions'`` — one regression per lab × condition.
        ``'daily'``      — one regression per lab × condition × day.

    Returns
    -------
    pd.DataFrame
        Regression statistics with columns ``lab``, ``condition``,
        ``date``, ``slope``, ``intercept``, ``r``, ``p_value``, ``se``.
    """
    pct_col  = f'{count_col}_pct'
    records  = []
    group_by = ['plate_idx', 'source name'] if mode == 'conditions' else \
               ['plate_idx', 'source name', 'date']

    for keys, group in drift_df.dropna(subset=[pct_col]).groupby(group_by):
        keys   = keys if isinstance(keys, tuple) else (keys,)
        x, y   = group['time_norm'].values, group[pct_col].values
        if len(x) < 2:
            continue
        slope, intercept, r, p, se = linregress(x, y)
        record = dict(zip(['lab', 'condition', 'date'], list(keys) + [None, None]))
        record.update(slope=slope, intercept=intercept, r=r, p_value=p, se=se)
        records.append(record)

    return pd.DataFrame(records)