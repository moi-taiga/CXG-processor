#!/usr/bin/env python3
"""
Cellxgene Census Demo

This script demonstrates how to use the cellxgene-census package directly
to query and download data from the CZ CELLxGENE Discover Census.
"""

import sys
import os
from pathlib import Path

try:
    import cellxgene_census
    from cellxgene_census import get_anndata
    import anndata
except ImportError as e:
    print(f"Error: cellxgene-census package not installed.")
    print(f"Please install it with: pip install cellxgene-census")
    print(f"Missing package: {e}")
    sys.exit(1)


def demo_basic_query():
    """Demonstrate a basic census query."""
    print("=" * 60)
    print("BASIC CENSUS QUERY DEMO")
    print("=" * 60)
    
    try:
        print("Querying human PBMC data from census...")
        
        # Use the correct API pattern with context manager
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            adata = cellxgene_census.get_anndata(
                census,
                "Homo sapiens",
                obs_value_filter='tissue_general=="blood"'
            )
        
        print(f"✓ Successfully queried data:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Save the data
        output_path = Path("demo_pbmc_sample.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error in basic query: {e}")
        return None


def demo_mouse_brain_query():
    """Demonstrate a mouse brain census query."""
    print("\n" + "=" * 60)
    print("MOUSE BRAIN CENSUS QUERY DEMO")
    print("=" * 60)
    
    try:
        print("Querying mouse brain data from census...")
        
        # Use the correct API pattern with context manager
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            adata = cellxgene_census.get_anndata(
                census,
                "Mus musculus",
                obs_value_filter='tissue_general=="brain"'
            )
        
        print(f"✓ Successfully queried data:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Save the data
        output_path = Path("demo_mouse_brain_sample.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error in mouse brain query: {e}")
        return None


def demo_data_exploration(adata, name):
    """Demonstrate basic data exploration."""
    print(f"\n{'='*60}")
    print(f"DATA EXPLORATION: {name}")
    print(f"{'='*60}")
    
    if adata is None:
        print("No data to explore")
        return
    
    print(f"Dataset: {name}")
    print(f"Shape: {adata.shape}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    
    # Show available observations (cell metadata)
    if adata.obs.columns.tolist():
        print(f"Cell annotations: {list(adata.obs.columns)}")
        
        # Show unique values for some key columns
        for col in ['tissue_general', 'cell_type', 'assay']:
            if col in adata.obs.columns:
                unique_vals = adata.obs[col].unique()
                print(f"  {col}: {list(unique_vals)}")
    
    # Show available variables (gene metadata)
    if adata.var.columns.tolist():
        print(f"Gene annotations: {list(adata.var.columns)}")
    
    # Show available layers
    if adata.layers:
        print(f"Available layers: {list(adata.layers.keys())}")
    
    # Show unstructured data
    if adata.uns:
        print(f"Unstructured data keys: {list(adata.uns.keys())}")


def demo_census_info():
    """Show information about available census data."""
    print("=" * 60)
    print("CENSUS INFORMATION")
    print("=" * 60)
    
    try:
        # Get census info using context manager
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            print("Available organisms:")
            for organism in census.census_info["soma"]["builds"]:
                print(f"  - {organism}")
            
            print(f"\nCensus version: {census.census_info['census_version']}")
            print(f"Release date: {census.census_info['release_date']}")
        
    except Exception as e:
        print(f"Error getting census info: {e}")


def main():
    """Main demonstration function."""
    print("Cellxgene Census Demo")
    print("=" * 50)
    print("This demo shows how to use cellxgene-census to query data")
    print("from the CZ CELLxGENE Discover Census.")
    print()
    
    # Show census information
    demo_census_info()
    
    # Basic query demo
    pbmc_data = demo_basic_query()
    demo_data_exploration(pbmc_data, "Human PBMC Sample")
    
    # Mouse brain query demo
    mouse_data = demo_mouse_brain_query()
    demo_data_exploration(mouse_data, "Mouse Brain Sample")
    
    print("\n" + "=" * 60)
    print("DEMO COMPLETE")
    print("=" * 60)
    print("Generated files:")
    
    for file_path in ["demo_pbmc_sample.h5ad", "demo_mouse_brain_sample.h5ad"]:
        if Path(file_path).exists():
            print(f"  ✓ {file_path}")
        else:
            print(f"  ✗ {file_path} (not created)")
    
    print("\nYou can now use these h5ad files with:")
    print("1. Regular cellxgene for visualization")
    print("2. Scanpy for analysis")
    print("3. Other single-cell analysis tools")


if __name__ == "__main__":
    main()
