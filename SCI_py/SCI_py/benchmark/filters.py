"""Row filters for the ProteoBench results table (one row per datapoint)."""


def filter_recommended(df, rec_dict):
    """Keep rows where comments matches the value in rec_dict for that software_name."""
    mask = df.apply(
        lambda r: rec_dict.get(r['software_name']) == r['comments'], axis=1
    )
    return df[mask]


def apply_version_filter(df, version_filter):
    """Keep rows matching version_filter; None means accept all versions."""
    def _match(row):
        """Return True if the row's software_version equals the target for its software_name (or no target is set)."""
        target = version_filter.get(row['software_name'])
        if target is None:
            return True
        return str(row['software_version']).strip() == str(target).strip()
    return df[df.apply(_match, axis=1)]
