import sys
import numpy as np
import scanpy as sc

def hartigans_diptest(x):
    try:
        from diptest import diptest
        dip, pval = diptest(x)
        return dip, pval
    except ImportError:
        return None, None

def print_h5ad_metrics(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
    print(f"File: {h5ad_path}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")

    # Count matrix
    X = adata.X
    if not isinstance(X, np.ndarray):
        X = X.toarray()
    total_entries = X.size
    zero_entries = np.sum(X == 0)
    print(f"Percent zero entries: {zero_entries / total_entries * 100:.2f}%")

    cell_sums = X.sum(axis=1)
    gene_sums = X.sum(axis=0)
    print(f"Mean counts per cell: {np.mean(cell_sums):.2f}")
    print(f"Median counts per cell: {np.median(cell_sums):.2f}")
    print(f"Mean counts per gene: {np.mean(gene_sums):.2f}")
    print(f"Median counts per gene: {np.median(gene_sums):.2f}")

    # Bimodality test on log(cell sums)
    log_cell_sums = np.log1p(cell_sums)
    dip, pval = hartigans_diptest(log_cell_sums)
    if dip is not None:
        print(f"Hartigan's Dip Test statistic (cell sums): {dip:.4f}, p-value: {pval:.4g}")
        if pval < 0.05:
            print("Gene expression per cell appears bimodal (p < 0.05).")
        else:
            print("Gene expression per cell does not appear bimodal (p >= 0.05).")
    else:
        print("Install 'diptest' package for bimodality test: pip install diptest")

    # Metadata checks
    print(f"Cell metadata columns: {list(adata.obs.columns)}")
    print(f"Gene metadata columns: {list(adata.var.columns)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python h5ad_metrics.py <path_to_h5ad>")
        sys.exit(1)
    print_h5ad_metrics(sys.argv[1])
