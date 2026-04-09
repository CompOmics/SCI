"""
SCI_py.read
===========
Reading and ingestion functions for single-cell proteomics metadata sources.

Three data sources are supported:

* **SDRF** – per-plate ``SDRF.xlsx`` files following the MAGE-TAB/SDRF format
  for proteomics (one file per plate, stored in sub-directories named by
  plate ID).
* **cellenONE** – isolation reports (``Reordered_*_isolated.xls``) exported
  by the cellenONE instrument.
* **Docx annotations** – manual QC annotations (dying, dividing, doublet,
  blank wells) extracted from plate preparation documentation.

Typical usage
-------------
>>> sdrf        = SCI_read.read_SDRF(ROOT_SDRF_IN)
>>> isolated_df = SCI_read.read_cellenONE(ROOT_CELLENONE_IN)
>>> metadata    = sdrf.merge(isolated_df,
...                 on=['comment[X position]', 'comment[Y position]',
...                     'comment[Plate Number]'],
...                 how='left')
>>> metadata    = SCI_read.add_docx_annotations(metadata)
"""

from __future__ import annotations

import os
from glob import glob

import pandas as pd

from ..constants.metadata import LAB_TO_PLATENR, PLATE_ANNOTATIONS
from .parse import FILENAME_CONVERTERS, xy_to_rc

__all__ = [
    "read_SDRF",
    "read_cellenONE",
    "add_docx_annotations",
]


# ---------------------------------------------------------------------------
# SDRF reader
# ---------------------------------------------------------------------------

def read_SDRF(root: str) -> pd.DataFrame:
    """
    Read and merge all per-plate ``SDRF.xlsx`` files found under ``root``.

    Expected directory layout::

        root/
          Wyss_Astral_1/
            SDRF.xlsx
          Broad_timsTOF_1/
            SDRF.xlsx
          ...

    In addition to concatenating the individual SDRF tables, this function:

    * Derives ``sample_name`` from ``comment[data file]`` using
      plate-specific filename normalisation rules.
    * Maps ``comment[Plate Number]`` from the plate ID.
    * Converts letter-based X positions (A–P) to 1-based integers.
    * Computes ``WellName`` in ``R{row}-C{col}`` notation.
    * Annotates and deduplicates reinjected samples (Wyss, Broad).
    * Corrects a known BYU source-name labelling error.

    Parameters
    ----------
    root:
        Root directory containing one sub-directory per plate.

    Returns
    -------
    pd.DataFrame
        Merged SDRF table with additional derived columns.
    """
    sdrf_list = []
    for sdrf_path in glob(os.path.join(root, "*", "SDRF.xlsx")):
        df       = pd.read_excel(sdrf_path)
        plate_id = os.path.basename(os.path.dirname(sdrf_path))
        df["plate_id"]       = plate_id
        df["institute_name"] = plate_id.split("_")[0]
        sdrf_list.append(df)

    if not sdrf_list:
        raise FileNotFoundError(f"No SDRF.xlsx files found under '{root}'")

    sdrf = pd.concat(sdrf_list, axis=0, ignore_index=True)

    # --- Derived columns ---
    sdrf["sample_name"] = sdrf.apply(
        lambda row: FILENAME_CONVERTERS[row["plate_id"]](
            row["comment[data file]"]
        ).replace("control", "Control"),
        axis=1,
    )
    sdrf["comment[Plate Number]"] = sdrf["plate_id"].map(LAB_TO_PLATENR)
    sdrf["WellName"] = sdrf.apply(
        lambda row: xy_to_rc(row["comment[X position]"], row["comment[Y position]"]),
        axis=1,
    )

    # Convert letter X positions (A–P) to 1-based integers
    _x_to_int = {letter: i + 1 for i, letter in enumerate("ABCDEFGHIJKLMNOP")}
    sdrf["comment[X position]"] = sdrf["comment[X position]"].map(_x_to_int)

    # --- Reinjection annotation ---
    sdrf["Reinjection"] = False
    sdrf = _annotate_wyss_reinjections(sdrf, root)
    sdrf = _annotate_broad_reinjections(sdrf)

    # Deduplicate: keep reinjection over original when both exist
    sdrf = (
        sdrf.sort_values("Reinjection", ascending=False)
            .drop_duplicates(subset=["sample_name", "plate_id"], keep="first")
    )

    # --- BYU source-name correction ---
    # BYU swapped QC and blank well labels during pipetting
    _byu_swap = {
        "Site specific blanks": "Site specific QC",
        "Site specific QC"    : "Site specific blanks",
    }
    byu_mask = sdrf["plate_id"] == "BYU_timsTOF_1"
    sdrf.loc[byu_mask, "source name"] = (
        sdrf.loc[byu_mask, "source name"].map(
            lambda x: _byu_swap.get(x, x)
        )
    )

    return sdrf


def _annotate_wyss_reinjections(sdrf: pd.DataFrame, root: str) -> pd.DataFrame:
    """
    Mark Wyss_Astral_1 reinjections using the CSI_Plate4run3 log file.

    The run3 file lists the second acquisition of cells that were originally
    run in run1; these are flagged as reinjections.
    """
    run3_path = os.path.join(root, "Wyss_Astral_1", "CSI_Plate4run3.csv")
    if not os.path.exists(run3_path):
        return sdrf

    run3 = pd.read_csv(run3_path, skiprows=1)
    wyss_reinjections = set(
        run3["File Name"].apply(lambda x: "_".join(x.split("_")[6:-1]))
    )

    mask = (sdrf["plate_id"] == "Wyss_Astral_1") & (
        sdrf["sample_name"].isin(wyss_reinjections)
    )
    sdrf.loc[mask, "Reinjection"] = True
    return sdrf


def _annotate_broad_reinjections(sdrf: pd.DataFrame) -> pd.DataFrame:
    """
    Mark Broad_timsTOF_1 reinjections identified by ``'redo'`` in the filename.
    """
    mask = (sdrf["plate_id"] == "Broad_timsTOF_1") & (
        sdrf["comment[data file]"].str.contains("redo", na=False)
    )
    sdrf.loc[mask, "Reinjection"] = True
    return sdrf


# ---------------------------------------------------------------------------
# cellenONE reader
# ---------------------------------------------------------------------------

def read_cellenONE(folder: str) -> pd.DataFrame:
    """
    Read and merge all cellenONE isolation reports found under ``folder``.

    Expected directory layout::

        folder/
          Plate_1/
            CellType_A/
              Reordered_*_isolated.xls
            CellType_B/
              Reordered_*_isolated.xls
          Plate_2/
            ...

    Parameters
    ----------
    folder:
        Root directory containing one sub-directory per plate.

    Returns
    -------
    pd.DataFrame
        Merged cellenONE table with added ``comment[Plate Number]``,
        ``comment[Sample Type]``, and ``image_path`` columns, and
        X/Y position columns renamed to the SDRF convention.
    """
    records = []

    for plate_path in glob(os.path.join(folder, "*")):
        plate_nr = int(os.path.basename(plate_path).split("_")[1])

        for sample_type_folder in glob(os.path.join(plate_path, "*")):
            sample_type = os.path.basename(sample_type_folder).split("_")[0]
            hits        = glob(os.path.join(sample_type_folder, "Reordered_*_isolated.xls"))

            if not hits:
                print(f"Not found — Plate {plate_nr}, Sample type '{sample_type}'")
                continue

            isolated_path = hits[0]
            df            = pd.read_csv(isolated_path, sep="\t")
            df["comment[Plate Number]"] = plate_nr
            df["comment[Sample Type]"]  = sample_type
            df["image_path"]            = df["ImageFile"].apply(
                lambda x: os.path.join(os.path.dirname(isolated_path), x.split('"')[1])
            )
            records.append(df)

    if not records:
        raise FileNotFoundError(f"No cellenONE isolation files found under '{folder}'")

    merged = pd.concat(records, ignore_index=True)

    # Drop unnamed index columns and rename position columns to SDRF convention
    merged = merged.loc[:, ~merged.columns.str.startswith("Unnamed")]
    merged['YPos'] = merged['YPos'].apply(lambda x: 25-x)
    merged = merged.rename(columns={
        "XPos": "comment[X position]",
        "YPos": "comment[Y position]",
    })

    return merged


# ---------------------------------------------------------------------------
# Docx annotation
# ---------------------------------------------------------------------------

def add_docx_annotations(sdrf: pd.DataFrame) -> pd.DataFrame:
    """
    Apply manual per-well QC annotations from plate preparation documentation.

    Annotations are defined in ``_PLATE_ANNOTATIONS`` and cover four boolean
    flags: ``Doublet``, ``Dividing``, ``Dying``, and ``CellInBlank``.  Wells
    not listed in the annotation table remain ``False``.

    Plates with only free-text notes (plates 1, 11, 12) are not annotated
    at the well level; their notes are recorded in ``_PLATE_NOTES`` for
    reference.

    Parameters
    ----------
    sdrf:
        Merged SDRF table containing ``comment[Plate Number]`` and
        ``WellName`` columns.

    Returns
    -------
    pd.DataFrame
        Input table with added boolean columns
        ``Doublet``, ``Dividing``, ``Dying``, ``CellInBlank``.
    """
    qc_flags = ["Doublet", "Dividing", "Dying", "CellInBlank"]
    for flag in qc_flags:
        sdrf[flag] = False

    valid_plates = set(sdrf["comment[Plate Number]"].unique())

    for plate_nr, flags in PLATE_ANNOTATIONS.items():
        if plate_nr not in valid_plates:
            continue
        plate_mask = sdrf["comment[Plate Number]"] == plate_nr
        for flag, wells in flags.items():
            well_mask = sdrf["WellName"].isin(wells)
            sdrf.loc[plate_mask & well_mask, flag] = True

    return sdrf