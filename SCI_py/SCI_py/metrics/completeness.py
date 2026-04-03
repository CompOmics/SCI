from __future__ import annotations

import numpy as np
from ..constants.maps import LABS


def compute_completeness(
    adata,
    conditions: list[str] | None = None,
    layer:      str              = '',
) -> tuple[np.ndarray, dict, dict]:
    """
    Compute feature- and sample-level completeness curves per lab.

    Parameters
    ----------
    adata:
        AnnData object with ``plate_idx`` in ``obs``.
    conditions:
        Source names to include.  ``None`` includes all.
    layer:
        Layer to use.  Empty string uses ``X``.

    Returns
    -------
    thresholds:
        Array of completeness thresholds in percent (0–100).
    retention_curves:
        Dict mapping lab → % features retained at each threshold.
    sample_curves:
        Dict mapping lab → % samples retained at each threshold.
    """
    obs = adata.obs.copy()
    if conditions is not None:
        obs = obs[obs['source name'].isin(conditions)]

    thresholds       = np.arange(0, 1.01, 0.01)
    retention_curves = {}
    sample_curves    = {}

    for lab in LABS:
        idx = obs[obs['plate_idx'] == lab].index
        if len(idx) == 0:
            continue

        X = adata[idx, :].layers[layer] if layer else adata[idx, :].X
        valid_samples  = np.sum(~np.isnan(X), axis=1) > 0
        valid_features = np.sum(~np.isnan(X), axis=0) > 0
        X = X[np.ix_(valid_samples, valid_features)]

        n_samples, n_features = X.shape
        comp_feature          = np.sum(~np.isnan(X), axis=0) / n_samples
        comp_sample           = np.sum(~np.isnan(X), axis=1) / n_features

        retention_curves[lab] = np.array([
            np.sum(comp_feature >= t) / n_features * 100
            for t in thresholds
        ])
        sample_curves[lab] = np.array([
            np.sum(comp_sample >= t) / n_samples * 100
            for t in thresholds
        ])

    return thresholds * 100, retention_curves, sample_curves


def compute_combined_curve(
    adata,
    conditions: list[str] | None  = None,
    layer:      str               = '',
    thresholds: np.ndarray | None = None,
    axis:       int               = 0,
) -> np.ndarray:
    """Compute completeness curve across all labs combined."""
    if thresholds is None:
        thresholds = np.arange(0, 101, 1)

    obs = adata.obs.copy()
    if conditions is not None:
        obs = obs[obs['source name'].isin(conditions)]

    X = adata[obs.index, :].layers[layer] if layer else adata[obs.index, :].X
    X = X[np.sum(~np.isnan(X), axis=1) > 0, :]

    n_samples, n_features = X.shape
    completeness = (
        np.sum(~np.isnan(X), axis=axis) /
        (n_samples if axis == 0 else n_features)
    )
    total = n_features if axis == 0 else n_samples
    return np.array([
        np.sum(completeness >= t / 100) / total * 100
        for t in thresholds
    ])


def store_completeness(
    adata,
    conditions:        list[str] | None = None,
    layer:             str              = '',
    conditions_suffix: str              = '',
) -> None:
    """
    Store per-sample and per-feature completeness in ``adata.obs`` and ``adata.var``.

    Parameters
    ----------
    adata:
        AnnData object with ``plate_idx`` in ``obs``.
    conditions:
        Source names to include.  ``None`` includes all.
    layer:
        Layer to use.  Empty string uses ``X``.
    conditions_suffix:
        Suffix appended to column names, e.g. ``'SC'`` or ``'QC'``.
    """
    suffix  = f'_{conditions_suffix}' if conditions_suffix else ''
    obs_col = f'completeness_sample{suffix}'
    adata.obs[obs_col] = np.nan

    obs = adata.obs.copy()
    if conditions is not None:
        obs = obs[obs['source name'].isin(conditions)]

    for lab in LABS:
        idx = obs[obs['plate_idx'] == lab].index
        if len(idx) == 0:
            continue

        X          = adata[idx, :].layers[layer] if layer else adata[idx, :].X
        n_samples  = X.shape[0]
        n_features = X.shape[1]

        adata.obs.loc[idx, obs_col]                      = np.sum(~np.isnan(X), axis=1) / n_features
        adata.var[f'completeness_feature_{lab}{suffix}'] = np.sum(~np.isnan(X), axis=0) / n_samples