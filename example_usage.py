#!/usr/bin/env python3
"""
Example usage of downloaded h5ad files with cellxgene package.

This script demonstrates how to:
1. Load downloaded h5ad files
2. Launch cellxgene for visualization
3. Basic data exploration
"""

import os
import sys
from pathlib import Path
import yaml

# Add the current directory to Python path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    import cellxgene_census
    import anndata
except ImportError as e:
    print(f"Error: Required packages not installed. Please run: pip install cellxgene-census anndata")
    print(f"Missing package: {e}")
    sys.exit(1)


def load_config(config_path: str = "config.yaml"):
    """Load configuration file."""
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Configuration file {config_path} not found")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML configuration: {e}")
        return None


def list_downloaded_files(output_dir: str):
    """List all downloaded h5ad files."""
    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"Output directory {output_dir} does not exist")
        return []
    
    h5ad_files = list(output_path.glob("*.h5ad"))
    return h5ad_files


def explore_dataset(file_path: Path):
    """Basic exploration of an h5ad dataset."""
    print(f"\n{'='*60}")
    print(f"EXPLORING DATASET: {file_path.name}")
    print(f"{'='*60}")
    
    try:
        # Load the dataset
        print("Loading dataset...")
        adata = anndata.read_h5ad(file_path)
        
        # Basic information
        print(f"Dataset shape: {adata.shape}")
        print(f"Number of cells: {adata.n_obs}")
        print(f"Number of genes: {adata.n_vars}")
        
        # Available observations (cell metadata)
        if adata.obs.columns.tolist():
            print(f"Cell annotations: {list(adata.obs.columns)}")
        
        # Available variables (gene metadata)
        if adata.var.columns.tolist():
            print(f"Gene annotations: {list(adata.var.columns)}")
        
        # Available layers
        if adata.layers:
            print(f"Available layers: {list(adata.layers.keys())}")
        
        # Available uns (unstructured data)
        if adata.uns:
            print(f"Unstructured data keys: {list(adata.uns.keys())}")
        
        return adata
        
    except Exception as e:
        print(f"Error loading dataset {file_path}: {e}")
        return None


def launch_cellxgene(file_path: Path, port: int = 5005):
    """Launch cellxgene for the specified dataset."""
    print(f"\nLaunching cellxgene for {file_path.name}...")
    print(f"Access the visualization at: http://localhost:{port}")
    print("Press Ctrl+C to stop the server")
    
    try:
        # Note: cellxgene-census doesn't have a launch function like cellxgene
        # You would need to use the regular cellxgene package for visualization
        print("Note: cellxgene-census is for data access, not visualization.")
        print("To visualize the data, you can:")
        print("1. Use the regular cellxgene package: cellxgene launch <file>")
        print("2. Use scanpy for analysis and plotting")
        print("3. Use the downloaded h5ad file with other visualization tools")
        
        # For now, just show the file path
        print(f"Dataset saved at: {file_path.absolute()}")
        
    except Exception as e:
        print(f"Error with dataset: {e}")


def main():
    """Main function demonstrating usage."""
    print("Cellxgene H5AD File Usage Example")
    print("=" * 50)
    
    # Load configuration
    config = load_config()
    if not config:
        print("Using default output directory: ./downloads")
        output_dir = "./downloads"
    else:
        output_dir = config.get('output_directory', './downloads')
    
    # List downloaded files
    downloaded_files = list_downloaded_files(output_dir)
    
    if not downloaded_files:
        print(f"No h5ad files found in {output_dir}")
        print("Please run the downloader first:")
        print("python cellxgene_downloader.py")
        return
    
    print(f"Found {len(downloaded_files)} downloaded h5ad files:")
    for i, file_path in enumerate(downloaded_files, 1):
        print(f"{i}. {file_path.name}")
    
    # Explore each dataset
    datasets = {}
    for file_path in downloaded_files:
        adata = explore_dataset(file_path)
        if adata is not None:
            datasets[file_path.name] = adata
    
    if not datasets:
        print("No valid datasets found")
        return
    
    # Interactive selection for cellxgene launch
    print(f"\n{'='*60}")
    print("CELLXGENE LAUNCH")
    print(f"{'='*60}")
    
    print("Available datasets for cellxgene launch:")
    for i, file_path in enumerate(downloaded_files, 1):
        print(f"{i}. {file_path.name}")
    
    try:
        choice = input(f"\nSelect dataset to launch with cellxgene (1-{len(downloaded_files)}) or 'q' to quit: ")
        
        if choice.lower() == 'q':
            print("Exiting...")
            return
        
        choice_idx = int(choice) - 1
        if 0 <= choice_idx < len(downloaded_files):
            selected_file = downloaded_files[choice_idx]
            
            # Get port number
            port_input = input("Enter port number (default: 5005): ").strip()
            port = int(port_input) if port_input else 5005
            
            # Launch cellxgene
            launch_cellxgene(selected_file, port)
        else:
            print("Invalid selection")
            
    except ValueError:
        print("Invalid input. Please enter a number.")
    except KeyboardInterrupt:
        print("\nExiting...")


if __name__ == "__main__":
    main()
