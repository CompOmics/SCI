from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm import tqdm


def build_design_matrix(
    obs:                    pd.DataFrame,
    numeric_covariates:     list[str] = (),
    categorical_covariates: list[str] = (),
) -> tuple[np.ndarray, list[str]]:
    """
    Build an intercept + numeric + reference-coded-categorical design matrix from
    columns of `obs`. Numeric covariates are mean-centered (does not change OLS
    residuals, improves conditioning).

    Returns
    -------
    (D, design_columns):
        D               : (n_samples x n_params) array
        design_columns  : names of D's columns, in order
    """
    n_samples = len(obs)
    blocks    = [np.ones((n_samples, 1))]
    columns   = ["intercept"]

    if numeric_covariates:
        numeric = obs[list(numeric_covariates)].astype(float)
        numeric = numeric - numeric.mean()
        if numeric.isna().any().any():
            raise ValueError(
                f"Numeric covariates contain NaN: "
                f"{numeric.columns[numeric.isna().any()].tolist()}"
            )
        blocks.append(numeric.to_numpy())
        columns.extend(numeric.columns)

    if categorical_covariates:
        dummies = pd.get_dummies(
            obs[list(categorical_covariates)].astype(str), drop_first=True
        ).astype(float)
        blocks.append(dummies.to_numpy())
        columns.extend(dummies.columns)

    D = np.column_stack(blocks)
    if not np.isfinite(D).all():
        raise ValueError("Design matrix contains non-finite values.")
    return D, columns


def regress_out_covariates(
    matrix:                 np.ndarray,
    obs:                    pd.DataFrame,
    numeric_covariates:     list[str] = (),
    categorical_covariates: list[str] = (),
    min_obs_per_param:      int       = 2,
    show_progress:          bool      = True,
) -> dict:
    """
    Per-feature OLS regression of `matrix` on the covariates in `obs`, keeping the
    per-feature mean (i.e. intercept) so the global feature-level order within a
    sample is preserved. Equivalent to scp's
    ``scpRemoveBatchEffect(..., intercept = TRUE)`` removing all listed effects.

    For each feature (column), the model is fit only on samples where that feature
    is observed (NaN is treated as missing, not zero). The corrected value is
    ``per-feature mean + residual``. Features observed in fewer than
    ``min_obs_per_param * n_params`` samples are left uncorrected (raw values kept)
    since there are too few points to separate the covariate effects from signal.
    Missing values remain missing.

    Parameters
    ----------
    matrix:
        (n_samples x n_features) array, e.g. log2 intensities. NaN marks a missing
        value.
    obs:
        Sample metadata, aligned row-for-row with `matrix` (same order, same
        length). Must contain every column named in `numeric_covariates` and
        `categorical_covariates`.
    numeric_covariates:
        Column names in `obs` to include as centered numeric terms.
    categorical_covariates:
        Column names in `obs` to include as reference-coded dummy terms.
    min_obs_per_param:
        Minimum observations per design-matrix parameter required to fit a feature.
    show_progress:
        Whether to display a tqdm progress bar over features.

    Returns
    -------
    Dict with:
        'corrected'      : (n_samples x n_features) corrected matrix
        'design_matrix'  : (n_samples x n_params) array
        'design_columns' : list of design matrix column names
        'min_obs'        : the min-observations threshold used
        'n_regressed'    : number of features fit
        'n_kept_raw'     : number of features left uncorrected (too few observations)
        'n_empty'        : number of all-missing features
    """
    if len(obs) != matrix.shape[0]:
        raise ValueError(
            f"obs has {len(obs)} rows but matrix has {matrix.shape[0]} samples."
        )

    D, design_columns = build_design_matrix(obs, numeric_covariates, categorical_covariates)
    n_params = D.shape[1]
    min_obs  = min_obs_per_param * n_params

    n_features = matrix.shape[1]
    corrected  = np.full_like(matrix, np.nan, dtype=float)
    n_regressed = n_kept_raw = n_empty = 0

    iterator = tqdm(range(n_features), desc="Regressing features") if show_progress else range(n_features)
    for j in iterator:
        y    = matrix[:, j]
        mask = ~np.isnan(y)
        n_obs = int(mask.sum())

        if n_obs == 0:
            n_empty += 1
            continue
        if n_obs < min_obs:
            corrected[mask, j] = y[mask]
            n_kept_raw += 1
            continue

        Dj = D[mask]
        yj = y[mask]
        beta, *_ = np.linalg.lstsq(Dj, yj, rcond=None)
        resid = yj - Dj @ beta
        corrected[mask, j] = resid + yj.mean()
        n_regressed += 1

    return {
        'corrected':      corrected,
        'design_matrix':  D,
        'design_columns': design_columns,
        'min_obs':        min_obs,
        'n_regressed':    n_regressed,
        'n_kept_raw':     n_kept_raw,
        'n_empty':        n_empty,
    }
