#!/usr/bin/env python3
"""
H5AD Data Type Checker

A Python script that uses anndata and scanpy to check if h5ad files
contain raw or normalized single-cell data.

Author: AI Assistant
Date: 2024
"""

import sys
import os
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

try:
    import anndata
    import scanpy as sc
    import scipy.sparse
except ImportError as e:
    print(f"Error: Required packages not installed.")
    print(f"Please install them with: pip install anndata scanpy scipy")
    print(f"Missing package: {e}")
    sys.exit(1)


class H5ADDataChecker:
    """Class to check h5ad files for raw vs normalized data."""
    
    def __init__(self):
        """Initialize the checker."""
        self.results = {}
    
    def check_data_type(self, file_path: str) -> dict:
        """
        Check if an h5ad file contains raw or normalized data.
        
        Args:
            file_path (str): Path to the h5ad file
            
        Returns:
            dict: Dictionary containing analysis results
        """
        print(f"\n{'='*60}")
        print(f"ANALYZING: {Path(file_path).name}")
        print(f"{'='*60}")
        
        try:
            # Load the h5ad file
            print("Loading h5ad file...")
            adata = anndata.read_h5ad(file_path)
            
            result = {
                'file_path': file_path,
                'file_name': Path(file_path).name,
                'shape': adata.shape,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'has_raw': False,
                'has_normalized': False,
                'data_layers': [],
                'raw_data_info': {},
                'normalized_data_info': {},
                'data_type_indicators': {},
                'recommendations': []
            }
            
            print(f"Dataset shape: {adata.shape}")
            print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
            
            # Check for raw data
            self._check_raw_data(adata, result)
            
            # Check for normalized data
            self._check_normalized_data(adata, result)
            
            # Check data layers
            self._check_data_layers(adata, result)
            
            # Analyze data characteristics
            self._analyze_data_characteristics(adata, result)
            
            # Generate recommendations
            self._generate_recommendations(result)
            
            # Print summary
            self._print_summary(result)
            
            return result
            
        except Exception as e:
            print(f"Error analyzing {file_path}: {e}")
            return {
                'file_path': file_path,
                'file_name': Path(file_path).name,
                'error': str(e)
            }
    
    def _check_raw_data(self, adata, result):
        """Check for raw data in the AnnData object."""
        print("\n--- RAW DATA ANALYSIS ---")
        
        # Check if raw data exists
        if adata.raw is not None:
            result['has_raw'] = True
            result['raw_data_info'] = {
                'shape': adata.raw.shape,
                'type': type(adata.raw.X).__name__,
                'sparse': scipy.sparse.issparse(adata.raw.X) if hasattr(adata.raw, 'X') else False
            }
            print(f"✓ Raw data found: {adata.raw.shape}")
            
            # Analyze raw data characteristics
            if hasattr(adata.raw, 'X') and adata.raw.X is not None:
                raw_data = adata.raw.X
                if scipy.sparse.issparse(raw_data):
                    raw_data = raw_data.toarray()
                
                result['raw_data_info'].update({
                    'min_value': float(np.min(raw_data)),
                    'max_value': float(np.max(raw_data)),
                    'mean_value': float(np.mean(raw_data)),
                    'median_value': float(np.median(raw_data)),
                    'std_value': float(np.std(raw_data)),
                    'zero_fraction': float(np.sum(raw_data == 0) / raw_data.size)
                })
                
                print(f"  Raw data range: {result['raw_data_info']['min_value']:.2f} - {result['raw_data_info']['max_value']:.2f}")
                print(f"  Raw data mean: {result['raw_data_info']['mean_value']:.2f}")
                print(f"  Zero fraction: {result['raw_data_info']['zero_fraction']:.2%}")
        else:
            print("✗ No raw data found")
    
    def _check_normalized_data(self, adata, result):
        """Check for normalized data in the AnnData object."""
        print("\n--- NORMALIZED DATA ANALYSIS ---")
        
        # Check main X matrix
        if adata.X is not None:
            print(f"Main data matrix (X): {adata.X.shape}")
            
            # Analyze X matrix characteristics
            x_data = adata.X
            if scipy.sparse.issparse(x_data):
                x_data = x_data.toarray()
            
            x_stats = {
                'min_value': float(np.min(x_data)),
                'max_value': float(np.max(x_data)),
                'mean_value': float(np.mean(x_data)),
                'median_value': float(np.median(x_data)),
                'std_value': float(np.std(x_data)),
                'zero_fraction': float(np.sum(x_data == 0) / x_data.size)
            }
            
            print(f"  X data range: {x_stats['min_value']:.2f} - {x_stats['max_value']:.2f}")
            print(f"  X data mean: {x_stats['mean_value']:.2f}")
            print(f"  Zero fraction: {x_stats['zero_fraction']:.2%}")
            
            # Determine if X is likely normalized
            is_normalized = self._is_likely_normalized(x_stats)
            result['normalized_data_info']['X'] = {
                'stats': x_stats,
                'is_normalized': is_normalized
            }
            
            if is_normalized:
                result['has_normalized'] = True
                print("  ✓ X matrix appears to be normalized")
            else:
                print("  ⚠ X matrix appears to be raw data")
    
    def _check_data_layers(self, adata, result):
        """Check for additional data layers."""
        print("\n--- DATA LAYERS ANALYSIS ---")
        
        if adata.layers:
            print(f"Found {len(adata.layers)} data layers:")
            for layer_name, layer_data in adata.layers.items():
                result['data_layers'].append(layer_name)
                
                if layer_data is not None:
                    layer_stats = {
                        'shape': layer_data.shape,
                        'type': type(layer_data).__name__,
                        'sparse': scipy.sparse.issparse(layer_data),
                        'min_value': float(np.min(layer_data)),
                        'max_value': float(np.max(layer_data)),
                        'mean_value': float(np.mean(layer_data)),
                        'zero_fraction': float(np.sum(layer_data == 0) / layer_data.size)
                    }
                    
                    print(f"  {layer_name}: {layer_data.shape}")
                    print(f"    Range: {layer_stats['min_value']:.2f} - {layer_stats['max_value']:.2f}")
                    print(f"    Mean: {layer_stats['mean_value']:.2f}")
                    
                    # Check if layer is normalized
                    is_normalized = self._is_likely_normalized(layer_stats)
                    if is_normalized:
                        result['has_normalized'] = True
                        print(f"    ✓ {layer_name} appears to be normalized")
                    else:
                        print(f"    ⚠ {layer_name} appears to be raw data")
                else:
                    print(f"  {layer_name}: None/empty")
        else:
            print("No additional data layers found")
    
    def _is_likely_normalized(self, stats):
        """
        Determine if data is likely normalized based on statistics.
        
        Args:
            stats (dict): Dictionary containing data statistics
            
        Returns:
            bool: True if data appears to be normalized
        """
        # Criteria for normalized data:
        # 1. Mean close to 0 or 1
        # 2. Standard deviation close to 1 (for z-score) or small (for log)
        # 3. Non-integer values
        # 4. Reasonable range (not extremely large)
        
        mean = stats['mean_value']
        std = stats['std_value']
        max_val = stats['max_value']
        
        # Check for log-normalized data (log1p)
        if 0 <= mean <= 2 and 0.5 <= std <= 2:
            return True
        
        # Check for z-score normalized data
        if abs(mean) < 1 and 0.5 <= std <= 2:
            return True
        
        # Check for min-max normalized data
        if 0 <= mean <= 1 and max_val <= 1:
            return True
        
        # Check for CPM/TPM normalized data (typically 0-1000 range)
        if 0 <= mean <= 1000 and max_val <= 10000:
            return True
        
        return False
    
    def _analyze_data_characteristics(self, adata, result):
        """Analyze overall data characteristics."""
        print("\n--- DATA CHARACTERISTICS ---")
        
        # Check for common normalization indicators
        indicators = {}
        
        # Check uns metadata for normalization info
        if 'log1p' in adata.uns:
            indicators['log1p_normalized'] = True
            print("✓ Log1p normalization detected in metadata")
        
        if 'normalization' in adata.uns:
            indicators['normalization_info'] = adata.uns['normalization']
            print(f"✓ Normalization info found: {adata.uns['normalization']}")
        
        # Check var metadata for gene info
        if 'gene_ids' in adata.var.columns:
            indicators['has_gene_ids'] = True
            print("✓ Gene IDs found in metadata")
        
        # Check obs metadata for cell info
        if 'cell_ids' in adata.obs.columns:
            indicators['has_cell_ids'] = True
            print("✓ Cell IDs found in metadata")
        
        result['data_type_indicators'] = indicators
    
    def _generate_recommendations(self, result):
        """Generate recommendations based on analysis."""
        print("\n--- RECOMMENDATIONS ---")
        
        recommendations = []
        
        if result['has_raw'] and result['has_normalized']:
            recommendations.append("✓ Both raw and normalized data available - ideal for analysis")
            recommendations.append("  Use adata.raw.X for raw counts, adata.X for normalized data")
        
        elif result['has_raw'] and not result['has_normalized']:
            recommendations.append("⚠ Only raw data available - consider normalizing before analysis")
            recommendations.append("  Use scanpy.pp.normalize_total() and scanpy.pp.log1p()")
        
        elif not result['has_raw'] and result['has_normalized']:
            recommendations.append("⚠ Only normalized data available - raw data may be needed for some analyses")
            recommendations.append("  Consider if you need raw counts for differential expression")
        
        else:
            recommendations.append("❌ No clear raw or normalized data detected")
            recommendations.append("  Check if the h5ad file contains the expected data")
        
        # Additional recommendations based on data layers
        if result['data_layers']:
            recommendations.append(f"✓ {len(result['data_layers'])} additional data layers available")
            recommendations.append(f"  Layers: {', '.join(result['data_layers'])}")
        
        for rec in recommendations:
            print(rec)
        
        result['recommendations'] = recommendations
    
    def _print_summary(self, result):
        """Print a summary of the analysis."""
        print(f"\n{'='*60}")
        print("SUMMARY")
        print(f"{'='*60}")
        
        print(f"File: {result['file_name']}")
        print(f"Shape: {result['shape']}")
        print(f"Raw data: {'✓ Yes' if result['has_raw'] else '✗ No'}")
        print(f"Normalized data: {'✓ Yes' if result['has_normalized'] else '✗ No'}")
        print(f"Data layers: {len(result['data_layers'])}")
        
        if result['data_layers']:
            print(f"  Available layers: {', '.join(result['data_layers'])}")
    
    def check_multiple_files(self, file_paths: list) -> dict:
        """
        Check multiple h5ad files.
        
        Args:
            file_paths (list): List of file paths to check
            
        Returns:
            dict: Dictionary containing results for all files
        """
        all_results = {}
        
        for file_path in file_paths:
            result = self.check_data_type(file_path)
            all_results[file_path] = result
        
        return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Check h5ad files for raw vs normalized data"
    )
    parser.add_argument(
        'files',
        nargs='+',
        help='H5AD file(s) to analyze'
    )
    parser.add_argument(
        '--output',
        '-o',
        help='Output file for results (JSON format)'
    )
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Check if files exist
    valid_files = []
    for file_path in args.files:
        if os.path.exists(file_path):
            valid_files.append(file_path)
        else:
            print(f"Warning: File {file_path} not found, skipping...")
    
    if not valid_files:
        print("No valid files to analyze!")
        sys.exit(1)
    
    print(f"Analyzing {len(valid_files)} h5ad file(s)...")
    
    # Create checker and analyze files
    checker = H5ADDataChecker()
    results = checker.check_multiple_files(valid_files)
    
    # Save results if output file specified
    if args.output:
        import json
        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        # Convert results to JSON-serializable format
        json_results = {}
        for file_path, result in results.items():
            json_results[file_path] = {}
            for key, value in result.items():
                if isinstance(value, dict):
                    json_results[file_path][key] = {}
                    for k, v in value.items():
                        json_results[file_path][key][k] = convert_numpy(v)
                else:
                    json_results[file_path][key] = convert_numpy(value)
        
        with open(args.output, 'w') as f:
            json.dump(json_results, f, indent=2)
        print(f"\nResults saved to {args.output}")


if __name__ == "__main__":
    main()
