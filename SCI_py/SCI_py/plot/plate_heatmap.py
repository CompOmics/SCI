from __future__ import annotations

import numpy as np
import anndata
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

from ..constants.plotting import CONDITION_COLORS, LAB_ORDER
from .utilities import _save


# ---------------------------------------------------------------------------
# Well colour map
# ---------------------------------------------------------------------------

def get_well_color(row: int, col: int) -> str:
    """
    Map a well position to its section colour based on plate layout.

    Parameters
    ----------
    row:
        Y position (1–16).
    col:
        X position (1–24).

    Returns
    -------
    str
        Hex colour string.
    """
    if row in (1, 16):                                    return '#6DBF8A'
    if row in (2, 15) and col in (1, 24):                 return '#6DBF8A'
    if col in (1, 24) and 3 <= row <= 14:                 return '#E8735A'
    if col in (2, 3)  and 3 <= row <= 14:                 return '#5B8DB8'
    if col in (22, 23) and 3 <= row <= 14:                return '#A8CCE8'
    if row == 2  and 2  <= col <= 12:                     return '#F9D776'
    if row == 15 and 13 <= col <= 23:                     return '#F9D776'
    if row == 2  and 13 <= col <= 21:                     return '#E8A598'
    if row == 15 and 4  <= col <= 12:                     return '#E8A598'
    if 3 <= row <= 14 and 4  <= col <= 12:                return '#9B72CF'
    if 3 <= row <= 14 and 13 <= col <= 21:                return '#F4A7B9'
    return '#ffffff'


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_plate_heatmap(
    adata,
    count_col:    str,
    extra_adata:  "anndata.AnnData | None" = None,
    label:        str                      = 'Features',
    save_dir:     str | None               = None,
    filename_fmt: str                      = '{lab}_plate_heatmap.png',
    show:         bool                     = True,
) -> None:
    """
    Plot per-well identification counts as a plate heatmap for each lab.

    Parameters
    ----------
    adata:
        AnnData object with ``plate_idx``, ``comment[X position]``,
        ``comment[Y position]``, ``source name``, and ``count_col`` in ``obs``.
    count_col:
        Column in ``obs`` containing the per-well count to display
        (e.g. ``'n_proteins'``).
    extra_adata:
        Optional second AnnData whose ``obs`` contains additional wells
        (e.g. QC wells) to overlay on the heatmap.  Must have the same
        columns as ``adata.obs``.
    label:
        Feature type label used in the plot title (e.g. ``'Proteins'``).
    save_dir:
        Directory for saved figures.  ``None`` disables saving.
    filename_fmt:
        Filename template.  ``{lab}`` is replaced with the lab name.
    show:
        Whether to call ``plt.show()``.
    """
    obs = _prepare_obs(adata)

    for lab in LAB_ORDER:
        plate_obs = obs[obs['plate_idx'] == lab]
        if plate_obs.empty:
            continue

        # Merge extra wells if provided
        if extra_adata is not None:
            extra_obs   = _prepare_obs(extra_adata)
            extra_plate = extra_obs[extra_obs['plate_idx'] == lab]
            plate_obs   = plate_obs._append(extra_plate).drop_duplicates(
                subset=['comment[X position]', 'comment[Y position]']
            )

        pivot = plate_obs.pivot(
            index   = 'comment[Y position]',
            columns = 'comment[X position]',
            values  = count_col,
        )
        condition_map = plate_obs.set_index(
            ['comment[Y position]', 'comment[X position]']
        )['source name'].to_dict()

        fig, ax = _draw_plate(pivot, condition_map, lab, label)

        if save_dir:
            _save(fig, f'{save_dir}/{filename_fmt.format(lab=lab)}')
        plt.show() if show else plt.close(fig)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _prepare_obs(adata) -> "pd.DataFrame":
    obs = adata.obs.copy()
    obs['comment[X position]'] = obs['comment[X position]'].astype(int)
    obs['comment[Y position]'] = obs['comment[Y position]'].astype(int)
    return obs


def _draw_plate(pivot, condition_map, lab, label):
    row_labels = pivot.index.tolist()
    col_labels = pivot.columns.tolist()
    n_rows, n_cols = len(row_labels), len(col_labels)
    values = pivot.values.astype(float)
    vmin, vmax = np.nanmin(values), np.nanmax(values)

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    for i, row in enumerate(row_labels):
        for j, col in enumerate(col_labels):
            val       = values[i, j]
            condition = condition_map.get((row, col))
            base_color = (CONDITION_COLORS.get(condition, '#ffffff')
                          if condition else get_well_color(col, row))

            rgb   = mcolors.to_rgb(base_color)
            alpha = (0.3 + 0.7 * (val - vmin) / (vmax - vmin)
                     if not np.isnan(val) else 0.1)

            ax.add_patch(plt.Rectangle(
                [j, n_rows - i - 1], 1, 1, color=(*rgb, alpha), zorder=1
            ))
            if not np.isnan(val):
                ax.text(j + 0.5, n_rows - i - 0.5, f'{int(val):,}',
                        ha='center', va='center', fontsize=6.5,
                        color='#333333', zorder=2)

    for boundary in [1, 2, 4, 13, 22, 24]:
        if boundary in row_labels:
            ax.axhline(n_rows - row_labels.index(boundary),
                       color='white', linewidth=2, zorder=3)
    for boundary in [1, 2, 3, 15, 16]:
        if boundary in col_labels:
            ax.axvline(col_labels.index(boundary),
                       color='white', linewidth=2, zorder=3)

    ax.set(xlim=(0, n_cols), ylim=(0, n_rows),
           xlabel='Column', ylabel='Row',
           title=f'{lab}  ·  {label} per well')
    ax.set_xticks([j + 0.5 for j in range(n_cols)])
    ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticks([n_rows - i - 0.5 for i in range(n_rows)])
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.title.set_fontsize(14)
    ax.title.set_fontweight('bold')

    ax.legend(
        handles=[mpatches.Patch(color=c, label=l)
                 for l, c in CONDITION_COLORS.items()],
        bbox_to_anchor=(1.01, 1), loc='upper left',
        fontsize=9, frameon=False,
    )

    return fig, ax