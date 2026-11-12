"""Build, align and inspect AnnData objects of search-engine results.

- :func:`df_to_anndata` -- pivot a long report table into an AnnData object
- :func:`align_obs_across_anndata` -- make ``obs`` identical across datasets
- :func:`replace_zeros_with_nan` -- treat zero intensities as missing
- :func:`inspect_anndata_size` -- memory footprint per AnnData component
- :func:`get_proteins_by_condition` -- proteins detected per lab for given conditions
- :func:`simplify_protein_group_ids` -- strip ``db|accession|entry`` ids to accessions
"""

from __future__ import annotations

import sys
from typing import cast

import anndata
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import issparse

from ..constants.labs import LABS


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
    """
    Give two AnnData objects the same samples (rows) and a shared ``obs``.

    Both intensity tables are reindexed to the union of their sample indices
    (samples missing from one dataset get all-NaN rows). ``obs`` is merged,
    preferring values from ``adata_1``. Features (``var``) are not aligned.

    Parameters
    ----------
    adata_1, adata_2:
        AnnData objects to align, e.g. as returned by :func:`df_to_anndata`.
    X_1, X_2:
        Their sample x feature intensity DataFrames (the second return value of
        :func:`df_to_anndata`), indexed by sample.

    Returns
    -------
    tuple[AnnData, AnnData]
        New AnnData objects with ``X``, ``obs`` and all ``layers`` reindexed to
        the union of samples, in the same order.
    """
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
        """Reindex every layer of `adata` to the sample order of `X_aligned`."""
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


def inspect_anndata_size(adata, name='adata'):
    """
    Estimate the memory used by each component of an AnnData object.

    Prints the AnnData shape and total estimated size.

    Parameters
    ----------
    adata:
        AnnData object to inspect.
    name:
        Label used in the printed summary.

    Returns
    -------
    pd.DataFrame
        One row per component (``X``, each layer, ``obs``, ``var``, each
        ``obsm``/``varm``/``uns`` entry) with columns ``component``, ``detail``
        and ``size_mb``, sorted by size (largest first).
    """
    rows = []

    # X matrix
    X = adata.X
    if scipy.sparse.issparse(X):
        size = X.data.nbytes + X.indices.nbytes + X.indptr.nbytes
    else:
        size = X.nbytes if X is not None else 0 
    rows.append({'component': 'X', 'detail': str(type(X).__name__), 'size_mb': size / 1e6})

    # layers
    for key, layer in adata.layers.items():
        if scipy.sparse.issparse(layer):
            size = layer.data.nbytes + layer.indices.nbytes + layer.indptr.nbytes
        else:
            size = layer.nbytes
        rows.append({'component': 'layers', 'detail': key, 'size_mb': size / 1e6})

    # obs / var
    for attr in ('obs', 'var'):
        df   = getattr(adata, attr)
        size = sum(col.memory_usage(deep=True) for _, col in df.items())
        rows.append({'component': attr, 'detail': f'{df.shape[0]} × {df.shape[1]}', 'size_mb': size / 1e6})

    # obsm / varm
    for attr in ('obsm', 'varm'):
        mapping = getattr(adata, attr)
        for key, arr in mapping.items():
            size = arr.nbytes if hasattr(arr, 'nbytes') else sys.getsizeof(arr)
            rows.append({'component': attr, 'detail': key, 'size_mb': size / 1e6})

    # uns
    for key, val in adata.uns.items():
        if hasattr(val, 'nbytes'):
            size = val.nbytes
        elif isinstance(val, pd.DataFrame):
            size = sum(col.memory_usage(deep=True) for _, col in val.items())
        else:
            size = sys.getsizeof(val)
        rows.append({'component': 'uns', 'detail': key, 'size_mb': size / 1e6})

    df = pd.DataFrame(rows).sort_values('size_mb', ascending=False).reset_index(drop=True)
    df['size_mb'] = df['size_mb'].round(2)

    total = df['size_mb'].sum()
    print(f'{name}: {adata.n_obs} obs × {adata.n_vars} vars')
    print(f'Total estimated size: {total:.1f} MB\n')

    return df


def get_proteins_by_condition(
    adata,
    source_name: list,
):
    """
    List the proteins detected in each lab for the given conditions.

    Contaminants (``var['Contaminant']``) are excluded. A protein counts as
    detected in a lab if it is non-NaN in ``X`` in at least one of that lab's
    samples. Prints the number of proteins per lab.

    Parameters
    ----------
    adata:
        Protein-level AnnData with ``plate_idx`` and ``source name`` in ``obs``
        and ``Contaminant`` and ``PG.ProteinAccessions`` in ``var``.
    source_name:
        ``source name`` values to include (e.g. blank conditions).

    Returns
    -------
    pd.DataFrame
        One column per lab holding its detected ``PG.ProteinAccessions``;
        shorter columns are padded with NaN.
    """
    adata_selection = adata[
        adata.obs['source name'].isin(source_name),
        ~adata.var['Contaminant']
    ]

    # For each lab, get the protein accessions detected in at least one blank cell
    protein_list = {}
    for lab in LABS:
        lab_idx = adata_selection.obs[adata_selection.obs['plate_idx'] == lab].index
        if len(lab_idx) == 0:
            continue

        X_lab        = adata_selection[lab_idx, :].X
        detected     = np.sum(~np.isnan(X_lab), axis=0) > 0
        protein_ids  = adata_selection.var.loc[detected, 'PG.ProteinAccessions'].tolist()
        protein_list[lab] = protein_ids
        print(f"{lab}: {len(protein_ids)} proteins detected in blanks")

    # Show as a DataFrame for easy inspection
    max_len = max(len(v) for v in protein_list.values())
    protein_df = pd.DataFrame({
        lab: pd.Series(proteins)
        for lab, proteins in protein_list.items()
    })
    return protein_df


def simplify_protein_group_ids(protein_ids) -> list[str]:
    """
    Simplify UniProt-style protein(-group) identifiers of the form
    ``db|accession|entry_name`` (``;``-joined for protein groups) down to just the
    accession(s), e.g.::

        'sp|A0A075B6H7|KV37_HUMAN;sp|P01624|KV315_HUMAN' -> 'A0A075B6H7;P01624'

    A member that doesn't match the ``db|accession|entry_name`` pattern (exactly two
    ``|``) is left unchanged.

    Parameters
    ----------
    protein_ids:
        Iterable of protein(-group) identifier strings.

    Returns
    -------
    List of simplified identifier strings, same length and order as `protein_ids`.
    """
    simplified = []
    for group in protein_ids:
        accessions = []
        for member in group.split(';'):
            parts = member.split('|')
            accessions.append(parts[1] if len(parts) == 3 else member)
        simplified.append(';'.join(accessions))
    return simplified
