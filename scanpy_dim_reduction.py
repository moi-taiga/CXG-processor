import scanpy as sc

def run_dim_reduction_and_plots(adata, pca_components=50, umap_neighbors=15):
    """
    Check for existing PCA/UMAP in a Scanpy object, run them if missing,
    and generate UMAP + PCA loadings plots.
 Steps performed:
    ----------------
    1. Check `.obsm` for existing 'X_pca' and 'X_umap' entries.
    2. If PCA is missing:
       - Normalize counts (total = 1e4) and log-transform.
       - Scale the data.
       - Run PCA with the specified number of components.
    3. If UMAP is missing:
       - Compute nearest neighbors using the PCA space.
       - Run UMAP with the specified number of neighbors.
    4. Plot UMAP (colored by 'celltype' if present in `.obs`).
    5. Plot PCA loadings for the first 3 principal components.

    Parameters
    ----------
    adata : AnnData
        The Scanpy object to process.
    pca_components : int, optional (default=50)
        Number of principal components to compute in PCA.
    umap_neighbors : int, optional (default=15)
        Number of neighbors to use in UMAP computation.

    Notes
    -----
    - This function modifies the AnnData object in-place.
    - If your data is already normalized/logged, the normalization steps will still run unless PCA exists.
    - For best results, run `sc.pp.highly_variable_genes()` before using this function on raw count data.
    - Plots are shown immediately; to save them, use Scanpy's `save` argument in `sc.pl` functions.
    """
    print("\n🔍 Checking dimensionality reductions in `adata.obsm`...")
    if adata.obsm:
        for key in adata.obsm.keys():
            print(f"   - Found: {key}")
    else:
        print("   - None found.")

    # PCA
    if "X_pca" not in adata.obsm:
        print("\n⚙️ Running PCA...")
        if 'n_top_genes' not in adata.var.keys():
            print("   → Normalizing and log-transforming data...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        print("   → Scaling data...")
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=pca_components, svd_solver='arpack')
    else:
        print("\n✅ PCA already present.")

    # UMAP
    if "X_umap" not in adata.obsm:
        print("\n⚙️ Running neighbors and UMAP...")
        sc.pp.neighbors(adata, n_neighbors=umap_neighbors, n_pcs=min(pca_components, adata.obsm['X_pca'].shape[1]))
        sc.tl.umap(adata)
    else:
        print("\n✅ UMAP already present.")

    # Plots
    print("\n📊 Plotting UMAP...")
    if "celltype" in adata.obs.columns:
        sc.pl.umap(adata, color="celltype", legend_loc='on data', frameon=False)
    else:
        sc.pl.umap(adata, frameon=False)

    print("\n📊 Plotting PCA loadings...")
    sc.pl.pca_loadings(adata, components=[1, 2, 3], show=True)
    
    print("\n✅ Dimensionality check and plots complete.")

	
