import scipy
import sys
import pandas as pd

def inspect_anndata_size(adata, name='adata'):
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