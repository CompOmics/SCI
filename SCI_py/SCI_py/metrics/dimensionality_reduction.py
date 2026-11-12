from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA


def select_pca_features(matrix: np.ndarray, min_completeness: float = 0.8) -> np.ndarray:
    """
    Boolean mask of features (columns) observed in at least `min_completeness`
    fraction of samples (rows). NaN marks a missing value.
    """
    return (~np.isnan(matrix)).mean(axis=0) >= min_completeness


def _impute_and_center(matrix: np.ndarray, feature_mask: np.ndarray) -> np.ndarray:
    """Per-feature mean-impute NaNs, then center. Returns (n_samples x n_selected_features)."""
    sub      = matrix[:, feature_mask]
    col_mean = np.nanmean(sub, axis=0)
    filled   = np.where(np.isnan(sub), col_mean, sub)
    return filled - col_mean


def compute_pca_embeddings(
    matrices:                 dict[str, np.ndarray],
    min_feature_completeness: float       = 0.8,
    n_components:             int         = 10,
    reference:                str | None  = None,
    random_state:             int         = 0,
) -> dict[str, dict]:
    """
    Fit PCA independently on each of several aligned log2-intensity matrices, e.g. a
    raw matrix and a batch-corrected version of it, so the resulting embeddings are
    directly comparable.

    All matrices must share the same samples (rows) and features (columns). The
    feature filter is computed once from `reference` and applied identically to every
    matrix, so every embedding is built on the same feature set.

    Parameters
    ----------
    matrices:
        Dict mapping a name (e.g. 'raw', 'corrected') to an (n_samples x n_features)
        array. NaN marks a missing value.
    min_feature_completeness:
        Minimum fraction of samples a feature must be observed in (in `reference`) to
        be included in the PCA.
    n_components:
        Number of principal components to compute.
    reference:
        Key of `matrices` whose missingness defines the shared feature mask. Defaults
        to the first key.
    random_state:
        Passed to sklearn's PCA for reproducibility.

    Returns
    -------
    Dict keyed like `matrices`, each value a dict with:
        'pcs'                      : (n_samples x n_components) array of scores
        'explained_variance_ratio' : (n_components,) array
        'model'                    : fitted sklearn PCA object
        'feature_mask'             : boolean mask of features used (shared across matrices)
    """
    if not matrices:
        raise ValueError("`matrices` must contain at least one entry.")

    keys      = list(matrices.keys())
    reference = reference or keys[0]
    if reference not in matrices:
        raise ValueError(f"reference '{reference}' not found in matrices: {keys}")

    feature_mask = select_pca_features(matrices[reference], min_feature_completeness)

    results = {}
    for name, matrix in matrices.items():
        Z     = _impute_and_center(matrix, feature_mask)
        model = PCA(n_components=n_components, random_state=random_state).fit(Z)
        results[name] = {
            'pcs':                      model.transform(Z),
            'explained_variance_ratio': model.explained_variance_ratio_,
            'model':                    model,
            'feature_mask':             feature_mask,
        }

    # every entry must share the exact same feature mask, since it is derived
    # once from `reference` and never recomputed per matrix.
    assert all(np.array_equal(feature_mask, r['feature_mask']) for r in results.values()), \
        "Internal error: embeddings do not share the same feature mask."

    return results


def compute_embeddings(
    matrices:                 dict[str, np.ndarray],
    methods:                  list[str]   = ('pca', 'umap', 'tsne'),
    min_feature_completeness: float       = 0.8,
    n_pca_components:         int         = 10,
    reference:                str | None  = None,
    random_state:             int         = 0,
) -> dict[str, dict]:
    """
    Compute one or more 2D-plottable embeddings for several aligned matrices, so
    every matrix/method combination is directly comparable.

    All methods are built on the same basis: `compute_pca_embeddings` is run first
    (same feature mask across matrices, as there), then 'umap' and/or 'tsne' — if
    requested — are fit on those PCA scores rather than on the raw feature space.
    This matches standard scRNA-seq/scp practice (PCA first, then UMAP/t-SNE on the
    PCA-reduced representation) and avoids a second, independent imputation step.

    Parameters
    ----------
    matrices:
        Dict mapping a name (e.g. 'Raw', 'Corrected') to an (n_samples x n_features)
        array. NaN marks a missing value.
    methods:
        Subset of ('pca', 'umap', 'tsne').
    min_feature_completeness / reference:
        See `compute_pca_embeddings`.
    n_pca_components:
        Number of PCA components to compute, and the number fed into UMAP/t-SNE.
    random_state:
        Passed to PCA/UMAP/t-SNE for reproducibility.

    Returns
    -------
    Dict keyed like `matrices`; each value is a dict keyed by method name, each
    holding at least:
        'coords'       : (n_samples x 2) array, ready to scatter-plot
        'feature_mask' : boolean mask of features used (shared across matrices/methods)
    'pca' additionally holds 'explained_variance_ratio' and 'model'.
    """
    pca_results = compute_pca_embeddings(
        matrices,
        min_feature_completeness=min_feature_completeness,
        n_components=n_pca_components,
        reference=reference,
        random_state=random_state,
    )

    results = {}
    for name, pca_res in pca_results.items():
        entry = {}
        if 'pca' in methods:
            entry['pca'] = {
                'coords':                   pca_res['pcs'][:, :2],
                'explained_variance_ratio': pca_res['explained_variance_ratio'],
                'model':                    pca_res['model'],
                'feature_mask':             pca_res['feature_mask'],
            }
        if 'umap' in methods:
            import umap
            reducer = umap.UMAP(n_components=2, random_state=random_state)
            entry['umap'] = {
                'coords':       reducer.fit_transform(pca_res['pcs']),
                'model':        reducer,
                'feature_mask': pca_res['feature_mask'],
            }
        if 'tsne' in methods:
            from sklearn.manifold import TSNE
            tsne = TSNE(n_components=2, init='pca', random_state=random_state)
            entry['tsne'] = {
                'coords':       tsne.fit_transform(pca_res['pcs']),
                'model':        tsne,
                'feature_mask': pca_res['feature_mask'],
            }
        results[name] = entry

    return results
