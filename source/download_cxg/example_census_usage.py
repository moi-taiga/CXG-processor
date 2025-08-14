#!/usr/bin/env python3
"""
Example: Cellxgene Census Usage

This script demonstrates the exact API usage pattern for cellxgene-census
as shown in the user's example.
"""

import sys
from pathlib import Path

try:
    import cellxgene_census
    import anndata
except ImportError as e:
    print(f"Error: cellxgene-census package not installed.")
    print(f"Please install it with: pip install cellxgene-census")
    print(f"Missing package: {e}")
    sys.exit(1)


def example_exact_pattern():
    """Example using the exact pattern from the user's code."""
    print("=" * 60)
    print("EXACT PATTERN EXAMPLE")
    print("=" * 60)
    
    try:
        # This matches the user's example exactly
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            adata = cellxgene_census.get_anndata(
                census,
                "Homo sapiens",
                obs_value_filter='dataset_id=="8e47ed12-c658-4252-b126-381df8d52a3d"',
            )
        
        print(f"✓ Successfully queried data:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Get unique cell types (as in the user's example)
        uq_cells = sorted(adata.obs["cell_type"].unique().tolist())
        print(f"  - Unique cell types: {uq_cells}")
        
        # Save the data
        output_path = Path("example_dataset.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error in exact pattern example: {e}")
        return None


def example_custom_query():
    """Example with custom query parameters."""
    print("\n" + "=" * 60)
    print("CUSTOM QUERY EXAMPLE")
    print("=" * 60)
    
    try:
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            # Custom query for blood tissue
            adata = cellxgene_census.get_anndata(
                census,
                "Homo sapiens",
                obs_value_filter='tissue_general=="blood"'
            )
        
        print(f"✓ Successfully queried blood tissue data:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Show available cell types
        if "cell_type" in adata.obs.columns:
            uq_cells = sorted(adata.obs["cell_type"].unique().tolist())
            print(f"  - Available cell types: {len(uq_cells)}")
            print(f"  - Sample cell types: {uq_cells[:5]}")  # Show first 5
        
        # Save the data
        output_path = Path("custom_blood_sample.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error in custom query example: {e}")
        return None


def example_mouse_data():
    """Example with mouse data."""
    print("\n" + "=" * 60)
    print("MOUSE DATA EXAMPLE")
    print("=" * 60)
    
    try:
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            # Query mouse brain data
            adata = cellxgene_census.get_anndata(
                census,
                "Mus musculus",
                obs_value_filter='tissue_general=="brain"'
            )
        
        print(f"✓ Successfully queried mouse brain data:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Show available cell types
        if "cell_type" in adata.obs.columns:
            uq_cells = sorted(adata.obs["cell_type"].unique().tolist())
            print(f"  - Available cell types: {len(uq_cells)}")
            print(f"  - Sample cell types: {uq_cells[:5]}")
        
        # Save the data
        output_path = Path("mouse_brain_sample.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error in mouse data example: {e}")
        return None


def main():
    """Main function to run all examples."""
    print("Cellxgene Census Examples")
    print("=" * 50)
    print("This script demonstrates the correct API usage pattern")
    print("for cellxgene-census, matching the user's example.")
    print()
    
    # Run the exact pattern example
    exact_data = example_exact_pattern()
    
    # Run custom query example
    custom_data = example_custom_query()
    
    # Run mouse data example
    mouse_data = example_mouse_data()
    
    print("\n" + "=" * 60)
    print("EXAMPLES COMPLETE")
    print("=" * 60)
    print("Generated files:")
    
    for file_path in ["example_dataset.h5ad", "custom_blood_sample.h5ad", "mouse_brain_sample.h5ad"]:
        if Path(file_path).exists():
            print(f"  ✓ {file_path}")
        else:
            print(f"  ✗ {file_path} (not created)")
    
    print("\nKey points about the API:")
    print("1. Use context manager: with cellxgene_census.open_soma() as census:")
    print("2. Pass census object to get_anndata()")
    print("3. Use organism names like 'Homo sapiens' or 'Mus musculus'")
    print("4. Use double quotes in filters: 'tissue_general==\"blood\"'")
    print("5. Always close the census connection (handled by context manager)")


if __name__ == "__main__":
    main()
