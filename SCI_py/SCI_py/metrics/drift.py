from __future__ import annotations

import pandas as pd

from scipy.stats import linregress
from ..constants.maps import LABS

_DATETIME_COL = 'comment[Gradient time]_datetime'

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