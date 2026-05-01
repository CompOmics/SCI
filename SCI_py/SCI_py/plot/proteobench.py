import pandas as pd
from ..constants.plotting import (
    LAB_COLORS,
    LAB_ORDER,
    CONDITION_COLORS,
    CONDITION_ORDER,
    SOFTWARE_COLORS,
    SOFTWARE_ORDER
)

def plot_benchmark_scatter(
    df:            pd.DataFrame,
    color_column:  str          = 'lab',
    color_map:     dict | None  = None,
    shape_column:  str | None   = None,
    symbol_map:    dict | None  = None,
    tool_filter:   str | None   = None,
    metric:        str          = 'mean',
    size:          int          = 70,
    title:         str | None   = None,
    save_dir:      str | None   = None,
    filename:      str          = 'benchmark_scatter.png',
    show:          bool         = True,
) -> None:
    """
    Scatter plot of precursor count vs absolute epsilon.

    Parameters
    ----------
    df:
        Long-format DataFrame with columns ``software_name``, ``nr_prec``,
        ``mean_abs_epsilon`` or ``median_abs_epsilon``, ``comments``,
        and ``color_column``.
    color_column:
        Column to use for point colors.
    color_map:
        Dict mapping ``color_column`` values to colors.
        Falls back to matplotlib color cycle if ``None``.
    shape_column:
        Column to use for point markers. ``None`` uses circles for all.
        Pass ``'comment_type'`` to split on comment suffix.
        Requires ``tool_filter`` to be set — separate legends are only
        supported for single-tool plots.
    symbol_map:
        Dict mapping ``shape_column`` values to matplotlib markers.
        Falls back to a default cycle if ``None``.
    tool_filter:
        If provided, filter to a single tool (``software_name``).
    metric:
        ``'mean'`` or ``'median'``.
    size:
        Marker size.
    title:
        Figure title. Auto-generated if ``None``.
    save_dir:
        Directory to save figures. ``None`` disables saving.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    if shape_column is not None and tool_filter is None:
        n_colors = df[color_column].nunique()
        n_shapes = df[shape_column].nunique() if shape_column else 0
        if n_colors * n_shapes > 20:
            warnings.warn(
                f"Large number of color × shape combinations ({n_colors} × {n_shapes}). "
                "The legend may be unreadable. Consider using tool_filter to subset the data."
            )

    x_col = f'{metric}_abs_epsilon'
    if x_col not in df.columns:
        raise ValueError(f"Column '{x_col}' not found. Use metric='mean' or metric='median'.")

    df = df.copy()

    if tool_filter is not None:
        df = df[df['software_name'] == tool_filter]
        if df.empty:
            raise ValueError(f"No rows found for tool_filter='{tool_filter}'.")

    # parse comment_type if needed
    if shape_column == 'comment_type' and 'comment_type' not in df.columns:
        df['comment_type'] = df['comments'].apply(lambda x: x.split('||')[-1])

    # ── color mapping ─────────────────────────────────────────────────────────
    if color_map is None:
        categories = df[color_column].unique()
        cycle      = plt.rcParams['axes.prop_cycle'].by_key()['color']
        color_map  = {cat: cycle[i % len(cycle)] for i, cat in enumerate(categories)}

    # ── symbol mapping ────────────────────────────────────────────────────────
    _default_markers = ['o', 's', 'D', '^', 'v', 'P', 'X', '*', 'h']
    if shape_column is not None and symbol_map is None:
        unique_shapes = df[shape_column].unique()
        symbol_map    = {val: _default_markers[i % len(_default_markers)]
                         for i, val in enumerate(unique_shapes)}

    fig, ax = plt.subplots(figsize=(8, 5))

    # ── plot points ───────────────────────────────────────────────────────────
    for group in df[color_column].unique():
        sub   = df[df[color_column] == group]
        color = color_map.get(group, 'black')

        if shape_column is not None:
            for shape_val, sub_shape in sub.groupby(shape_column):
                ax.scatter(
                    sub_shape[x_col],
                    sub_shape['nr_prec'],
                    c      = color,
                    marker = symbol_map.get(shape_val, 'o'),
                    s      = size,
                    alpha  = 0.8,
                )
        else:
            ax.scatter(
                sub[x_col],
                sub['nr_prec'],
                c      = color,
                marker = 'o',
                s      = size,
                alpha  = 0.8,
                label  = group,
            )

    # ── legends ───────────────────────────────────────────────────────────────
    if shape_column is not None:
        # color legend
        color_handles = [
            plt.scatter([], [], c=color_map.get(g, 'black'), s=size, alpha=0.8, label=g)
            for g in df[color_column].unique()
        ]
        color_legend = ax.legend(
            handles        = color_handles,
            title          = color_column,
            bbox_to_anchor = (1.02, 1),
            loc            = 'upper left',
            borderaxespad  = 0,
            fontsize       = 8,
        )
        ax.add_artist(color_legend)

        # shape legend
        shape_handles = [
            plt.scatter([], [], c='black', marker=symbol_map.get(s, 'o'),
                        s=size, alpha=0.8, label=s)
            for s in df[shape_column].unique()
        ]
        ax.legend(
            handles        = shape_handles,
            title          = shape_column,
            bbox_to_anchor = (1.02, 0.4),
            loc            = 'upper left',
            borderaxespad  = 0,
            fontsize       = 8,
        )
    else:
        ax.legend(
            title          = color_column,
            bbox_to_anchor = (1.02, 1),
            loc            = 'upper left',
            borderaxespad  = 0,
            fontsize       = 8,
        )

    ax.set(
        xlabel = f'{metric.capitalize()} absolute ε',
        ylabel = 'Number of precursors',
        title  = title or (tool_filter or 'All tools'),
    )
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, color='lightgrey', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        fig.tight_layout()

    if save_dir:
        _save(fig, os.path.join(save_dir, filename))

    plt.show() if show else plt.close(fig)