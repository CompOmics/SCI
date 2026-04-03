import pandas as pd

def compute_rt_proportions(
    report,
    sample_ids:   "pd.Index | list",
    rt_col:       str,
    conditions:   list[str] | None = None,
    bins:         int              = 20,
    normalize_rt: bool             = False,
) -> pd.DataFrame:
    """
    Compute per-sample retention time proportion distributions.

    Parameters
    ----------
    report:
        Long-format report DataFrame with columns ``sample_id``,
        ``source name``, ``plate_idx``, and ``rt_col``.
    sample_ids:
        Sample IDs to include.
    rt_col:
        Column containing retention time values.
    conditions:
        Source names to include.  ``None`` includes all.
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
    mask = report['sample_id'].isin(sample_ids)
    if conditions is not None:
        mask &= report['source name'].isin(conditions)

    df = report[mask].copy()

    if normalize_rt:
        rt_min = df.groupby('sample_id')[rt_col].transform('min')
        rt_max = df.groupby('sample_id')[rt_col].transform('max')
        df[rt_col] = (df[rt_col] - rt_min) / (rt_max - rt_min)

    df['bin_edges'] = pd.cut(df[rt_col], bins=bins)

    mapping = df[['sample_id', 'source name', 'plate_idx']].drop_duplicates().set_index('sample_id')

    props = (
        df.groupby('sample_id')['bin_edges']
          .value_counts(normalize=True)
          .rename('proportion')
          .reset_index()
    )
    props['bin_mid']     = props['bin_edges'].apply(lambda x: x.mid)
    props['source name'] = props['sample_id'].map(mapping['source name'])
    props['plate_idx']   = props['sample_id'].map(mapping['plate_idx'])

    return props