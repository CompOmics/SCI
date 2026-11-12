"""Clean, impute and type the columns of the merged metadata table.

The entry point is :func:`parse_metadata_columns`, which runs the steps below
in order on the merged SDRF + cellenONE + docx-annotation table:

1. :func:`parse_temporal_columns` -- split date/time columns
2. :func:`compute_injection_order` -- order runs per plate by acquisition time
3. :func:`impute_missing_values` -- fill known gaps per lab
4. :func:`convert_column_types` -- cast columns to their final dtypes

Typical usage
-------------
>>> sdrf        = read_SDRF(ROOT_SDRF_IN)
>>> isolated_df = read_cellenONE(ROOT_CELLENONE_IN)
>>> metadata    = sdrf.merge(isolated_df,
...                 on=['comment[X position]', 'comment[Y position]',
...                     'comment[Plate Number]'],
...                 how='left')
>>> metadata    = add_docx_annotations(metadata)
>>> metadata    = parse_metadata_columns(metadata)
"""

from __future__ import annotations

import pandas as pd

from ..constants.labs import LAB_NAME_MAPPING
from ..constants.samples import (
    _X_TO_COL,
    N_ROWS,
    SCALAR_IMPUTATIONS,
    SOURCE_NAME_MAPS,
)

__all__ = [
    "parse_metadata_columns",
    "parse_temporal_columns",
    "compute_injection_order",
    "impute_missing_values",
    "convert_column_types",
    "xy_to_rc",
]


# ---------------------------------------------------------------------------
# Plate geometry
# ---------------------------------------------------------------------------

def xy_to_rc(x: str, y: int) -> str:
    """
    Convert physical plate coordinates (X, Y) to table notation (R, C).

    Parameters
    ----------
    x:
        Column letter, e.g. ``"A"`` or ``"a"``.  Must be in A–P.
    y:
        Row number 1–24 (Y = 1 corresponds to the bottom of the plate).

    Returns
    -------
    str
        Table label such as ``"R24-C1"``.

    Raises
    ------
    ValueError
        If ``x`` or ``y`` are outside the valid range.
    """
    x = x.upper()
    if x not in _X_TO_COL:
        raise ValueError(f"X must be a letter A–P, got '{x}'")
    if not (1 <= y <= N_ROWS):
        raise ValueError(f"Y must be 1–{N_ROWS}, got {y}")

    col_idx = _X_TO_COL[x]
    row_idx = N_ROWS + 1 - y          # Y=1 → R24,  Y=24 → R1
    return f"R{row_idx}-C{col_idx}"


# ---------------------------------------------------------------------------
# Temporal parsing
# ---------------------------------------------------------------------------

def parse_temporal_columns(sdrf: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """
    Parse a raw datetime string column into three derived columns.

    The source column is expected to follow the ``%Y%m%d-%H%M%S`` format
    (inexact matching is enabled to tolerate minor variations).

    Derived columns
    ---------------
    ``{column_name}_date``
        Date component, floored to midnight (``pd.Timestamp``).
    ``{column_name}_time``
        Decimal hours as a float (e.g. ``14.5`` = 14:30:00).
    ``{column_name}_datetime``
        Full ``pd.Timestamp``.

    Parameters
    ----------
    sdrf:
        Metadata table containing ``column_name``.
    column_name:
        Name of the raw datetime column.

    Returns
    -------
    pd.DataFrame
        Input table with three additional columns appended.
    """
    parsed = pd.to_datetime(sdrf[column_name], format="%Y%m%d-%H%M%S", exact=False)
    sdrf[f"{column_name}_date"]     = parsed.dt.floor("D")
    sdrf[f"{column_name}_time"]     = (
        parsed.dt.hour
        + parsed.dt.minute / 60
        + parsed.dt.second / 3600
    )
    sdrf[f"{column_name}_datetime"] = parsed
    return sdrf


# ---------------------------------------------------------------------------
# Injection order
# ---------------------------------------------------------------------------

def compute_injection_order(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a normalised injection order in ``[0, 1]`` within each plate.

    The order is derived from ``comment[Gradient time]_datetime``.
    Missing datetimes are imputed with the plate-level median injection order;
    plates with a single observed datetime are set to ``0.5``.

    Parameters
    ----------
    metadata:
        Table containing ``plate_id`` and ``comment[Gradient time]_datetime``.

    Returns
    -------
    pd.DataFrame
        Input table with an added ``injection_order`` column.
    """
    df = metadata.copy()
    df["_datetime"] = pd.to_datetime(
        df["comment[Gradient time]_datetime"], errors="coerce"
    )

    def _plate_order(group: pd.DataFrame) -> pd.Series:
        """Scale acquisition-time ranks within one plate to 0-1; missing times get the plate median (0.5 if <= 1 time)."""
        dt         = group["_datetime"]
        n_observed = dt.notna().sum()
        if n_observed <= 1:
            return pd.Series(0.5, index=group.index)
        ranks  = dt.rank(method="average", na_option="keep")
        order  = (ranks - 1) / (n_observed - 1)
        median = order[dt.notna()].median()
        return order.where(dt.notna(), median)

    df["injection_order"] = (
        df.groupby("plate_id", group_keys=False).apply(_plate_order)
    ).fillna(0.5)

    return df.drop(columns=["_datetime"])


# ---------------------------------------------------------------------------
# Missing value imputation
# ---------------------------------------------------------------------------

def _impute_gradient_time_KI_Astral_2(df: pd.DataFrame) -> pd.DataFrame:
    """
    Impute missing ``comment[Gradient time]`` columns for plate ``KI_Astral_2``.

    Strategy
    --------
    1. **Date** – copied from ``comment[thawing timepoint]_date``.
    2. **Time** – exact copy from ``KI_Astral_1``, matched on ``sample_name``.
    3. **Datetime** – reconstructed from the imputed date and time.
    """
    ki2_mask = df["plate_id"] == "KI_Astral_2"
    ki1_mask = df["plate_id"] == "KI_Astral_1"

    # 1 — date
    df.loc[ki2_mask, "comment[Gradient time]_date"] = (
        df.loc[ki2_mask, "comment[thawing timepoint]_date"].values
    )

    # 2 — time matched by sample_name from KI_Astral_1
    ki1_time = (
        df[ki1_mask]
        .dropna(subset=["comment[Gradient time]_time"])
        .set_index("sample_name")["comment[Gradient time]_time"]
    )
    df.loc[ki2_mask, "comment[Gradient time]_time"] = (
        df.loc[ki2_mask, "sample_name"].map(ki1_time).values
    )

    # 3 — reconstruct datetime
    def _to_timedelta(t: float) -> pd.Timedelta:
        """Convert a time in decimal hours to a Timedelta (hours, minutes, seconds)."""
        h = int(t)
        m = int((t - h) * 60)
        s = int(((t - h) * 60 - m) * 60)
        return pd.Timedelta(hours=h, minutes=m, seconds=s)

    for idx in df[ki2_mask].index:
        date = df.at[idx, "comment[Gradient time]_date"]
        time = df.at[idx, "comment[Gradient time]_time"]
        if pd.notna(date) and pd.notna(time):
            df.at[idx, "comment[Gradient time]_datetime"] = (
                pd.Timestamp(date) + _to_timedelta(time)
            )

    return df


def impute_missing_values(sdrf: pd.DataFrame) -> pd.DataFrame:
    """
    Impute missing values in the merged metadata table.

    Two strategies are applied in order:

    1. **Scalar imputation** – fill with fixed sentinel values defined in
       ``SCALAR_IMPUTATIONS`` (``-1`` for numeric fields, SDRF controlled
       vocabulary strings for text fields).
    2. **Source-name lookup** – overwrite SDRF characteristic columns
       (organism, organism part, cell type, disease, cell line) from
       ``SOURCE_NAME_MAPS``.  These columns are fully overwritten rather
       than only filling NaN, as they are deterministically derived from
       ``source name``.

    Plate-specific gradient time imputation (``_impute_gradient_time_KI_Astral_2``)
    is handled earlier in ``parse_metadata_columns``, before injection order
    computation.

    Parameters
    ----------
    sdrf:
        Merged metadata table (SDRF + cellenONE + docx annotations).

    Returns
    -------
    pd.DataFrame
        Table with missing values filled.
    """
    df = sdrf.copy()

    # 1 — scalar
    for col, value in SCALAR_IMPUTATIONS.items():
        if col in df.columns:
            df[col] = df[col].fillna(value)

    # 2 — source-name derived columns (full overwrite)
    for col, mapping in SOURCE_NAME_MAPS.items():
        if col in df.columns:
            df[col] = df["source name"].map(mapping).fillna("not applicable")

    return df


# ---------------------------------------------------------------------------
# Type coercion
# ---------------------------------------------------------------------------

def convert_column_types(sdrf: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce expected dtypes on parsed metadata columns.

    Parameters
    ----------
    sdrf:
        Metadata table after parsing and imputation.

    Returns
    -------
    pd.DataFrame
        Table with corrected column types.
    """
    int_cols = [
        "comment[X position]",
        "comment[Y position]",
    ]
    str_cols = [
        "comment[Plate Number]",
        "comment[Gradient time]",
        "comment[thawing timepoint]",
        "comment[plate arrival]",
        "comment[plate arrival]_date",
        "comment[plate arrival]_datetime",
        "comment[Gradient time]_date",
        "comment[Gradient time]_datetime",
        "comment[thawing timepoint]_date",
        "comment[thawing timepoint]_datetime"
    ]
    bool_cols = [
        "Reinjection",
        "Doublet",
        "Dividing",
        "Dying",
        "CellInBlank",
    ]

    for col in int_cols:
        if col in sdrf.columns:
            sdrf[col] = sdrf[col].astype(int)
    for col in str_cols:
        if col in sdrf.columns:
            sdrf[col] = sdrf[col].astype(str)
    for col in bool_cols:
        if col in sdrf.columns:
            sdrf[col] = sdrf[col].astype(bool)

    return sdrf


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def parse_metadata_columns(sdrf: pd.DataFrame) -> pd.DataFrame:
    """
    Full metadata parsing pipeline for merged SDRF + cellenONE tables.

    Steps
    -----
    1. Normalise ``sample_name`` to lower-case.
    2. Parse temporal columns (plate arrival, gradient time, thawing timepoint).
    3. Derive ``instrument_type``, ``sample_id``, and ``plate_idx``.
    4. Normalise ``source name`` labels (Cell A → K562, Cell B → HEK293T).
    5. Impute plate-specific gradient times (``KI_Astral_2`` from ``KI_Astral_1``).
    6. Compute normalised injection order within each plate.
    7. Impute missing values (scalar sentinels, source-name lookups).
    8. Enforce column types.
    9. Drop the duplicate ``comment[plate arrival].1`` column introduced by
       some SDRF files.

    Parameters
    ----------
    sdrf:
        Merged metadata table (output of ``SCI_read`` functions).

    Returns
    -------
    pd.DataFrame
        Fully parsed and annotated metadata table.
    """
    sdrf = sdrf.copy()

    # 1 — sample name normalisation
    sdrf["sample_name"] = sdrf["sample_name"].str.lower()

    # 2 — temporal columns
    for col in [
        "comment[plate arrival]",
        "comment[Gradient time]",
        "comment[thawing timepoint]",
    ]:
        sdrf = parse_temporal_columns(sdrf, column_name=col)

    # 3 — derived identifiers
    sdrf["instrument_type"] = sdrf["plate_id"].str.split("_").str[1]
    sdrf["sample_id"]       = sdrf["sample_name"] + "|" + sdrf["plate_id"]
    sdrf["plate_idx"]       = sdrf["plate_id"].map(LAB_NAME_MAPPING)

    # 4 — source name normalisation
    sdrf["source name"] = (
        sdrf["source name"]
        .str.replace("Cell A", "K562",    regex=False)
        .str.replace("Cell B", "HEK293T", regex=False)
    )
    sdrf = sdrf.rename(
        columns={
            'Date        ': 'Date',
            'Time    ': 'Time'
        }
    )

    # 5 — plate-specific gradient time imputation (must precede injection order)
    sdrf = _impute_gradient_time_KI_Astral_2(sdrf)

    # 6 — injection order
    sdrf = compute_injection_order(sdrf)

    # 7 — missing value imputation
    sdrf = impute_missing_values(sdrf)

    # 7 — type coercion
    sdrf = convert_column_types(sdrf)

    # 8 — drop duplicate column introduced by some SDRF files
    if "comment[plate arrival].1" in sdrf.columns:
        sdrf = sdrf.drop(columns=["comment[plate arrival].1"])

    return sdrf
