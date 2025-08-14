#!/usr/bin/env python3
"""
Download Specific Datasets by Dataset ID

This script demonstrates how to download specific datasets from cellxgene-census
using their dataset IDs.
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


def download_dataset_by_id(dataset_id: str, organism: str = "Homo sapiens", output_name: str = None):
    """
    Download a specific dataset by its dataset_id.
    
    Args:
        dataset_id (str): The specific dataset ID to download
        organism (str): Organism name (default: "Homo sapiens")
        output_name (str): Name for the output file (default: dataset_id)
    
    Returns:
        bool: True if successful, False otherwise
    """
    if output_name is None:
        output_name = f"dataset_{dataset_id[:8]}"  # Use first 8 chars of dataset_id
    
    print(f"Downloading dataset: {dataset_id}")
    print(f"Organism: {organism}")
    print(f"Output name: {output_name}")
    
    try:
        # Use the correct API pattern with context manager
        CENSUS_VERSION = "2023-12-15"
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            # Query the specific dataset
            adata = cellxgene_census.get_anndata(
                census,
                organism,
                obs_value_filter=f'dataset_id=="{dataset_id}"'
            )
        
        print(f"✓ Successfully downloaded dataset:")
        print(f"  - Shape: {adata.shape}")
        print(f"  - Cells: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        # Show available cell types if present
        if "cell_type" in adata.obs.columns:
            uq_cells = sorted(adata.obs["cell_type"].unique().tolist())
            print(f"  - Available cell types: {len(uq_cells)}")
            print(f"  - Sample cell types: {uq_cells[:5]}")
        
        # Save the data
        output_path = Path(f"{output_name}.h5ad")
        adata.write_h5ad(output_path)
        print(f"  - Saved to: {output_path.absolute()}")
        
        return True, adata
        
    except Exception as e:
        print(f"✗ Error downloading dataset {dataset_id}: {e}")
        return False, None


def main():
    """Main function to demonstrate specific dataset downloads."""
    print("Download Specific Datasets by Dataset ID")
    print("=" * 50)
    
    # Example dataset IDs (you can replace these with your own)
    example_datasets = [
        {
            "dataset_id": "8e47ed12-c658-4252-b126-381df8d52a3d",
            "organism": "Homo sapiens",
            "output_name": "example_dataset_1"
        },
        # Add more datasets here
        # {
        #     "dataset_id": "your-dataset-id-here",
        #     "organism": "Homo sapiens",
        #     "max_cells": 500,
        #     "output_name": "your_dataset_name"
        # }
    ]
    
    print(f"Found {len(example_datasets)} datasets to download")
    print()
    
    successful_downloads = []
    
    for i, dataset_info in enumerate(example_datasets, 1):
        print(f"Processing dataset {i}/{len(example_datasets)}:")
        print("-" * 40)
        
        success, adata = download_dataset_by_id(
            dataset_id=dataset_info["dataset_id"],
            organism=dataset_info["organism"],
            output_name=dataset_info["output_name"]
        )
        
        if success:
            successful_downloads.append(dataset_info["output_name"])
        
        print()
    
    print("=" * 50)
    print("DOWNLOAD SUMMARY")
    print("=" * 50)
    print(f"Successfully downloaded: {len(successful_downloads)}/{len(example_datasets)} datasets")
    
    if successful_downloads:
        print("Downloaded files:")
        for name in successful_downloads:
            file_path = Path(f"{name}.h5ad")
            if file_path.exists():
                print(f"  ✓ {file_path}")
            else:
                print(f"  ✗ {file_path} (not found)")
    
    print("\nTo add your own datasets:")
    print("1. Find the dataset_id from the cellxgene-census data")
    print("2. Add it to the example_datasets list in this script")
    print("3. Or use the main downloader with config.yaml")


if __name__ == "__main__":
    main()
