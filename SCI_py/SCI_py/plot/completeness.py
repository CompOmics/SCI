"""
Completeness plots: feature- and sample-level retention curves per lab.

Part of the dataset statistics / QC-metrics analysis
(Notebooks/2_DatasetStatistics/02_QC_Metrics). The curves are computed in
``SCI_py.metrics.completeness``.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from ..constants.plotting import LAB_COLORS
from .export import save_figure


def plot_completeness(
    thresholds:       np.ndarray,
    retention_curves: dict,
    sample_curves:    dict,
    combined_feature: np.ndarray,
    combined_sample:  np.ndarray,
    label:            str        = 'Features',
    conditions:       list[str] | None = None,
    save_dir:         str | None = None,
    filename:         str        = 'completeness.png',
    show:             bool       = True,
):
    """
    Plot feature- and sample-level completeness curves per lab.

    Parameters
    ----------
    thresholds:
        Array of thresholds in percent (0–100), as returned by
        ``compute_completeness``.
    retention_curves:
        Dict mapping lab → % features retained, as returned by
        ``compute_completeness``.
    sample_curves:
        Dict mapping lab → % samples retained, as returned by
        ``compute_completeness``.
    combined_feature:
        Combined feature-level curve across all labs, as returned by
        ``compute_combined_curve``.
    combined_sample:
        Combined sample-level curve across all labs, as returned by
        ``compute_combined_curve``.
    label:
        Feature type label for axis titles (e.g. ``'Proteins'``).
    conditions:
        Used only for the figure title.
    save_dir:
        Directory for saved figures.
    filename:
        Output filename.
    show:
        Whether to call ``plt.show()``.
    """
    panels = [
        (retention_curves, combined_feature,
         f'% {label} retained',
         f'Min. % of samples a {label[:-1].lower()} must be observed in',
         f'{label}-level completeness'),
        (sample_curves, combined_sample,
         '% samples retained',
         f'Min. % of {label.lower()} a sample must have observed',
         'Sample-level completeness'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    for ax, (curves, combined, ylabel, xlabel, title) in zip(axes, panels):
        seen = set()
        for lab, values in curves.items():
            color = LAB_COLORS.get(lab, 'grey')
            ax.plot(thresholds, values, color=color, linewidth=1.8,
                    label=lab if lab not in seen else '_nolegend_')
            seen.add(lab)

        ax.plot(thresholds, combined, color='black', linewidth=1.8,
                linestyle='--', label='All', zorder=10)

        ax.set(xlabel=xlabel, ylabel=ylabel, title=title,
               xlim=(0, 100), ylim=(0, 105))
        ax.axhline(50, color='grey', linewidth=0.8, linestyle=':', alpha=0.6)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8, loc='lower left')

    cond_label = ', '.join(conditions) if conditions else 'all conditions'
    fig.suptitle(f'Completeness — {label} ({cond_label})', fontweight='bold')

    if save_dir:
        save_figure(fig, f'{save_dir}/{filename}')

    plt.show() if show else plt.close(fig)