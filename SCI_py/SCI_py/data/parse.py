import pandas as pd
import anndata
from typing import cast
import numpy as np
from scipy.sparse import issparse

def df_to_anndata(
    df:                pd.DataFrame,
    feature_col:       str,
    value_col:         str,
    sample_col:        str,
    obs_cols:          list[str],
    var_cols:          list[str],
) -> tuple["anndata.AnnData", pd.DataFrame]:
    """
    Convert a long-format dataframe to AnnData.

    Parameters
    ----------
    df          : long-format dataframe
    feature_col : column to use as features (var axis), e.g. 'FG.PrecursorIon'
    value_col   : column containing intensity values
    sample_col  : column to use as observations (obs axis), e.g. 'sample_id'
    obs_cols    : columns to include in obs
    var_cols    : columns to include in var
    """
    df = df.dropna(subset=[feature_col]).drop_duplicates(
        subset=[sample_col, feature_col, value_col]
    )

    X = df.pivot(index=sample_col, columns=feature_col, values=value_col)

    obs_cols = [c for c in obs_cols if c!=sample_col]
    obs = (df.drop_duplicates(sample_col)
             .set_index(sample_col)
             [obs_cols]
             .reindex(X.index))

    var_cols = [c for c in var_cols if c != feature_col]
    var = (df.drop_duplicates(feature_col)
             .set_index(feature_col)
             [[c for c in var_cols if c in df.columns]]
             .reindex(X.columns))

    adata = anndata.AnnData(X=X, obs=obs, var=var)
    return adata, X


def align_obs_across_anndata(
    adata_1: "anndata.AnnData",
    X_1:     pd.DataFrame,
    adata_2: "anndata.AnnData",
    X_2:     pd.DataFrame,
) -> tuple["anndata.AnnData", "anndata.AnnData"]:

    all_samples = X_1.index.union(X_2.index)

    X_1_aligned = X_1.reindex(all_samples)
    X_2_aligned = X_2.reindex(all_samples)

    obs_1 = cast(pd.DataFrame, adata_1.obs)
    obs_2 = cast(pd.DataFrame, adata_2.obs)
    var_1 = cast(pd.DataFrame, adata_1.var)
    var_2 = cast(pd.DataFrame, adata_2.var)

    obs_all = (obs_1.reindex(all_samples)
                    .combine_first(obs_2.reindex(all_samples)))

    def _align_layers(adata, X_aligned):
        return {
            key: pd.DataFrame(
                layer.toarray() if hasattr(layer, 'toarray') else layer,
                index   = adata.obs_names,
                columns = adata.var_names,
            ).reindex(X_aligned.index)
            for key, layer in adata.layers.items()
        }

    adata_1 = anndata.AnnData(
        X      = X_1_aligned,
        obs    = obs_all.loc[X_1_aligned.index],
        var    = var_1,
        layers = _align_layers(adata_1, X_1_aligned),
    )
    adata_2 = anndata.AnnData(
        X      = X_2_aligned,
        obs    = obs_all.loc[X_2_aligned.index],
        var    = var_2,
        layers = _align_layers(adata_2, X_2_aligned),
    )
    return adata_1, adata_2


def replace_zeros_with_nan(adata, layer=None, dtype=np.float32, verbose=True):
    """Replace zeros with NaN in an AnnData matrix.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to modify in place.
    layer : str, optional
        Name of the layer to modify. If None, operates on `adata.X`.
    dtype : numpy dtype, default np.float32
        Float dtype to cast the result to. Must support NaN.
    verbose : bool, default True
        Whether to print the number of zeros replaced.

    Returns
    -------
    AnnData
        The same AnnData object, modified in place.
    """
    X = adata.layers[layer] if layer is not None else adata.X

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X).copy()

    n_zeros = int((X == 0).sum())
    n_total = X.size

    if verbose:
        target = f"layers['{layer}']" if layer is not None else "X"
        print(
            f"Replacing {n_zeros} zeros out of {n_total} values in {target} "
            f"({100 * n_zeros / n_total:.2f}%)"
        )

    X = X.astype(dtype, copy=False)
    X[X == 0] = np.nan

    if layer is not None:
        adata.layers[layer] = X
    else:
        adata.X = X

    return adata