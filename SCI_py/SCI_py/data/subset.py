import numpy as np
import pandas as pd

from ..constants.maps import LABS

def get_proteins_by_condition(
    adata,
    source_name: list,
):
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