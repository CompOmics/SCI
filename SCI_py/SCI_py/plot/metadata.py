import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import os
from .constants import CONDITION_COLORS, LAB_COLORS


def plot_cell_types_per_plate(ax, metadata, plate_order):
    ct_counts = (metadata.groupby(['plate_idx', 'source name'])
                 .size().unstack(fill_value=0)
                 .reindex(plate_order, fill_value=0))

    bottom = np.zeros(len(plate_order))
    for ct in ct_counts.columns:
        color = CONDITION_COLORS.get(ct, '#aaaaaa')
        ax.bar(plate_order, ct_counts[ct], bottom=bottom,
               color=color, label=ct, width=0.6)
        bottom += ct_counts[ct].values

    ax.set_title('Cell types per plate', fontweight='bold')
    ax.set_xlabel('Plate ID')
    ax.set_ylabel('Number of cells')
    ax.set_xticklabels(plate_order, rotation=45, ha='right')
    ax.legend(title='Source name', frameon=False, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def plot_qc_flags_per_plate(ax, metadata, plate_order):
    qc_cols   = ['Reinjection', 'Doublet', 'Dividing', 'CellInBlank']
    qc_colors = ['#e07b54', '#5b8db8', '#76b876', '#b87bb8']

    qc_counts = (metadata.groupby('plate_idx')[qc_cols]
                 .sum()
                 .reindex(plate_order, fill_value=0))

    x      = np.arange(len(plate_order))
    n_cols = len(qc_cols)
    width  = 0.18

    for i, (col, color) in enumerate(zip(qc_cols, qc_colors)):
        ax.bar(x + i * width - (n_cols - 1) * width / 2,
               qc_counts[col], width=width, color=color, label=col)

    ax.set_title('QC flags per plate', fontweight='bold')
    ax.set_xlabel('Plate ID')
    ax.set_ylabel('Count')
    ax.set_xticks(x)
    ax.set_xticklabels(plate_order, rotation=45, ha='right')
    ax.legend(frameon=False, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def plot_acquisition_timeline(ax, metadata, plate_order):
    for plate in plate_order:
        sub = (metadata[metadata['plate_idx'] == plate]
               .dropna(subset=['comment[Gradient time]_datetime'])
               .copy())
        sub['comment[Gradient time]_datetime'] = pd.to_datetime(
            sub['comment[Gradient time]_datetime'])
        sub = sub.sort_values('comment[Gradient time]_datetime')
        color = LAB_COLORS.get(plate, '#aaaaaa')
        ax.scatter(sub['comment[Gradient time]_datetime'],
                   [plate] * len(sub),
                   color=color, s=8, alpha=0.6, label=plate)

    ax.set_title('Acquisition timeline per plate', fontweight='bold')
    ax.set_xlabel('Acquisition datetime')
    ax.set_ylabel('Plate ID')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', rotation=35)


def plot_plates_per_instrument(ax, metadata):
    instrument_counts = (metadata.groupby('instrument_type')['plate_idx']
                         .nunique()
                         .sort_values(ascending=False))

    colors = ['#6aaed6', '#f0b84a', '#5bbfaa', '#f0956e', '#dda0cc']
    ax.bar(instrument_counts.index, instrument_counts.values,
           color=colors[:len(instrument_counts)], width=0.5)

    ax.set_title('Plates per instrument type', fontweight='bold')
    ax.set_xlabel('Instrument type')
    ax.set_ylabel('Number of plates')
    ax.set_xticklabels(instrument_counts.index, rotation=35, ha='right')
    ax.spines[['top', 'right']].set_visible(False)

    for i, v in enumerate(instrument_counts.values):
        ax.text(i, v + 0.05, str(v), ha='center', fontweight='bold', fontsize=10)


def plot_plate_arrival(ax, metadata, plate_order):
    arrival = (metadata.groupby('plate_idx')['comment[plate arrival]_datetime']
               .first()
               .dropna()
               .apply(pd.to_datetime)
               .reindex(plate_order)
               .dropna()
               .sort_values())

    colors = [LAB_COLORS.get(p, '#aaaaaa') for p in arrival.index]
    ax.barh(range(len(arrival)),
            [(d - arrival.min()).days for d in arrival.values],
            color=colors, height=0.6)
    ax.set_yticks(range(len(arrival)))
    ax.set_yticklabels(arrival.index)
    ax.set_title('Plate arrival (days from first arrival)', fontweight='bold')
    ax.set_xlabel('Days from first plate arrival')
    ax.spines[['top', 'right']].set_visible(False)

    for i, (plate, date) in enumerate(arrival.items()):
        ax.text((date - arrival.min()).days + 0.3, i,
                date.strftime('%Y-%m-%d'), va='center', fontsize=8)


def plot_metadata_overview(metadata, save_path=None):
    plate_order = sorted(metadata['plate_idx'].unique())

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2, 0])
    ax5 = fig.add_subplot(gs[2, 1])

    plot_cell_types_per_plate(ax1,  metadata, plate_order)
    plot_qc_flags_per_plate(ax2,    metadata, plate_order)
    plot_acquisition_timeline(ax3,  metadata, plate_order)
    plot_plates_per_instrument(ax4, metadata)
    plot_plate_arrival(ax5,         metadata, plate_order)

    fig.suptitle('Single-Cell Proteomics — Metadata Overview',
                 fontsize=16, fontweight='bold', y=1.01)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', transparent=True)
        print(f"Saved to {save_path}")

    plt.show()
    return fig


def plot_diameter_per_plate_and_condition(ax, metadata, plate_order):
    """Violin + boxplot of cell diameter per plate, split by source name."""
    
    # Ensure Diameter is numeric
    metadata = metadata.copy()
    metadata['Diameter'] = pd.to_numeric(metadata['Diameter'], errors='coerce')
    
    conditions = ['K562', 'HEK293T']
    n_cond     = len(conditions)
    offsets    = np.linspace(-0.2, 0.2, n_cond)

    for j, cond in enumerate(conditions):
        sub   = metadata[metadata['source name'] == cond]
        color = CONDITION_COLORS.get(cond, '#aaaaaa')

        valid_plates = [(i, p) for i, p in enumerate(plate_order)
                        if len(sub[sub['plate_idx'] == p]['Diameter'].dropna()) > 1]

        if not valid_plates:
            continue

        indices, plates = zip(*valid_plates)
        groups  = [sub[sub['plate_idx'] == p]['Diameter'].dropna().values
                   for p in plates]
        x_pos   = np.array(indices, dtype=float) + offsets[j]

        parts = ax.violinplot(
            groups,
            positions   = x_pos,
            widths      = 0.35,
            showmedians = False,
            showextrema = False
        )
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.6)

        bp = ax.boxplot(
            groups,
            positions    = x_pos,
            widths       = 0.12,
            patch_artist = True,
            showfliers   = False,
            medianprops  = dict(color='black', linewidth=1.5)
        )
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.8)

        ax.plot([], [], color=color, linewidth=6, alpha=0.7, label=cond)

    ax.set_title('Cell diameter per plate', fontweight='bold')
    ax.set_xlabel('Plate ID')
    ax.set_ylabel('Diameter (µm)')
    ax.set_xticks(np.arange(len(plate_order)))
    ax.set_xticklabels(plate_order, rotation=45, ha='right')
    ax.legend(title='Source name', frameon=False, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def plot_injection_order_per_plate(ax, metadata, plate_order):
    """
    Distribution of injection order per plate per day, colored by source name.
    Plates with only one unique injection order value are skipped.
    Each plate x day combination is shown as a separate row.
    """
    # Add acquisition date
    df = metadata.copy()
    df['acq_date'] = pd.to_datetime(
        df['comment[Gradient time]_datetime'], errors='coerce').dt.date

    # Build ordered list of (plate, date) groups — skip plates with 1 unique order
    groups = []
    for plate in plate_order:
        sub = df[df['plate_idx'] == plate]

        # Skip plate if injection order has only one unique value
        n_unique = sub['injection_order'].dropna().nunique()
        if n_unique <= 1:
            continue

        # Get sorted unique dates for this plate
        dates = sorted(sub['acq_date'].dropna().unique())
        for date in dates:
            groups.append((plate, date))

    if not groups:
        ax.text(0.5, 0.5, 'No data to display',
                ha='center', va='center', transform=ax.transAxes)
        return

    conditions = sorted(df['source name'].dropna().unique())

    for row_idx, (plate, date) in enumerate(groups):
        sub = df[(df['plate_idx'] == plate) & (df['acq_date'] == date)]

        for cond in conditions:
            csub  = sub[sub['source name'] == cond]['injection_order'].dropna()
            color = CONDITION_COLORS.get(cond, '#aaaaaa')
            ax.scatter(csub, [row_idx] * len(csub),
                       color=color, s=8, alpha=0.6)

    # Y axis labels: "plate | date"
    y_labels = [f"{plate}\n{str(date)}" for plate, date in groups]
    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(y_labels, fontsize=7)

    # Separate plates visually with horizontal lines
    prev_plate = None
    for row_idx, (plate, date) in enumerate(groups):
        if prev_plate is not None and plate != prev_plate:
            ax.axhline(row_idx - 0.5, color='grey',
                       linewidth=0.8, linestyle='--', alpha=0.4)
        prev_plate = plate

    # Legend
    for cond in conditions:
        ax.scatter([], [], color=CONDITION_COLORS.get(cond, '#aaaaaa'),
                   s=20, label=cond)

    ax.set_title('Injection order per plate × day by cell type',
                 fontweight='bold')
    ax.set_xlabel('Injection order [0–1]')
    ax.set_ylabel('Plate × acquisition day')
    ax.axvline(0.5, color='grey', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.legend(title='Source name', frameon=False, fontsize=8, loc='upper right')
    ax.spines[['top', 'right']].set_visible(False)

    # Adjust figure height dynamically based on number of rows
    fig = ax.get_figure()
    fig.set_figheight(max(6, len(groups) * 0.5))


def plot_plate_position(ax, metadata, plate_order):
    """
    Scatter of X vs Y position colored by source name.
    Shows physical plate layout and cell type arrangement.
    Uses a grid of subplots — one per plate.
    """
    n_plates = len(plate_order)
    n_cols   = min(3, n_plates)
    n_rows   = int(np.ceil(n_plates / n_cols))

    # Clear the passed ax and use it as a container
    ax.set_visible(False)
    fig = ax.get_figure()

    # Find position of ax in figure to place subgrid
    pos      = ax.get_position()
    sub_axes = []
    for i, plate in enumerate(plate_order):
        row  = i // n_cols
        col  = i  % n_cols
        left = pos.x0 + col * pos.width  / n_cols
        bot  = pos.y0 + (n_rows - row - 1) * pos.height / n_rows
        sub_ax = fig.add_axes([left, bot,
                               pos.width / n_cols * 0.85,
                               pos.height / n_rows * 0.85])
        sub   = metadata[metadata['plate_idx'] == plate]
        for cond in sorted(sub['source name'].dropna().unique()):
            csub  = sub[sub['source name'] == cond]
            color = CONDITION_COLORS.get(cond, '#aaaaaa')
            sub_ax.scatter(csub['comment[X position]'],
                           csub['comment[Y position]'],
                           color=color, s=15, alpha=0.7)
        sub_ax.set_title(plate, fontsize=7, fontweight='bold')
        sub_ax.set_xlabel('X', fontsize=6)
        sub_ax.set_ylabel('Y', fontsize=6)
        sub_ax.tick_params(labelsize=6)
        sub_ax.spines[['top', 'right']].set_visible(False)
        sub_axes.append(sub_ax)

    # Shared legend
    conditions = sorted(metadata['source name'].dropna().unique())
    handles    = [Line2D([0], [0], marker='o', color='w',
                             markerfacecolor=CONDITION_COLORS.get(c, '#aaaaaa'),
                             markersize=6, label=c)
                  for c in conditions]
    fig.legend(handles=handles, title='Source name', frameon=False,
               fontsize=7, loc='lower right')


def plot_acquisition_date_heatmap(ax, metadata, plate_order):
    """
    Heatmap of acquisition dates per plate.
    Rows = plates, columns = unique acquisition dates,
    values = number of cells acquired on that date.
    """
    df = metadata.copy()
    df['acq_date'] = pd.to_datetime(
        df['comment[Gradient time]_datetime'], errors='coerce').dt.date

    pivot = (df.groupby(['plate_idx', 'acq_date'])
               .size()
               .unstack(fill_value=0)
               .reindex(plate_order, fill_value=0))

    pivot.columns = [str(c) for c in pivot.columns]

    im = ax.imshow(pivot.values, aspect='auto', cmap='Blues',
                   interpolation='nearest')

    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(np.arange(len(plate_order)))
    ax.set_yticklabels(plate_order, fontsize=8)
    ax.set_title('Cells acquired per plate per day', fontweight='bold')
    ax.set_xlabel('Acquisition date')
    ax.set_ylabel('Plate ID')

    # Annotate cells with counts > 0
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if val > 0:
                ax.text(j, i, str(val), ha='center', va='center',
                        fontsize=6,
                        color='white' if val > pivot.values.max() * 0.6
                        else 'black')

    plt.colorbar(im, ax=ax, label='Number of cells', shrink=0.8)


def plot_cell_qc_overview(metadata, save_path=None):
    """
    Master function — assembles all cell-level QC figures.
    """
    plate_order = sorted(metadata['plate_idx'].unique())

    fig = plt.figure(figsize=(22, 26))
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])   # plate position uses full row
    ax4 = fig.add_subplot(gs[2, :])   # heatmap uses full row

    plot_diameter_per_plate_and_condition(ax1, metadata, plate_order)
    plot_injection_order_per_plate(ax2,       metadata, plate_order)
    plot_plate_position(ax3,                  metadata, plate_order)
    plot_acquisition_date_heatmap(ax4,        metadata, plate_order)

    fig.suptitle('Cell-level QC — metadata overview',
                 fontsize=16, fontweight='bold', y=1.01)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', transparent=True)
        print(f"Saved to {save_path}")

    plt.show()
    return fig