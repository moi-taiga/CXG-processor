import scanpy as sc

def check_scanpy_metadata(adata):
    """
    Check which metadata columns are present in a Scanpy AnnData object.
    Special focus on the 'celltype' column.
    
    Parameters
    ----------
    adata : AnnData
        The Scanpy object to check.
    """
    print("Metadata columns available in .obs:")
    print(list(adata.obs.columns))
    print()

    if "celltype" in adata.obs.columns:
        print("'celltype' column found ✅")
        print("Available cell types:")
        print(sorted(adata.obs["celltype"].unique()))
    else:
        print("'celltype' column NOT found ❌")

