#!/usr/bin/env python3
"""
Example: Using the H5AD Data Checker

This script demonstrates how to use the check_h5ad_data.py script
to analyze h5ad files for raw vs normalized data.
"""

import sys
from pathlib import Path

# Import the checker
from check_h5ad_data import H5ADDataChecker


def example_single_file():
    """Example of checking a single h5ad file."""
    print("Example: Checking a single h5ad file")
    print("=" * 50)
    
    # Example file path (replace with your actual file)
    file_path = "example.h5ad"
    
    if not Path(file_path).exists():
        print(f"File {file_path} not found. Please provide a valid h5ad file path.")
        return
    
    # Create checker and analyze
    checker = H5ADDataChecker()
    result = checker.check_data_type(file_path)
    
    print(f"\nAnalysis complete for {file_path}")
    return result


def example_multiple_files():
    """Example of checking multiple h5ad files."""
    print("Example: Checking multiple h5ad files")
    print("=" * 50)
    
    # Example file paths (replace with your actual files)
    file_paths = [
        "file1.h5ad",
        "file2.h5ad",
        "file3.h5ad"
    ]
    
    # Filter to existing files
    existing_files = [f for f in file_paths if Path(f).exists()]
    
    if not existing_files:
        print("No valid h5ad files found. Please provide valid file paths.")
        return
    
    # Create checker and analyze
    checker = H5ADDataChecker()
    results = checker.check_multiple_files(existing_files)
    
    print(f"\nAnalysis complete for {len(existing_files)} files")
    return results


def example_with_custom_analysis():
    """Example with custom analysis of results."""
    print("Example: Custom analysis of results")
    print("=" * 50)
    
    # Example file path
    file_path = "example.h5ad"
    
    if not Path(file_path).exists():
        print(f"File {file_path} not found. Please provide a valid h5ad file path.")
        return
    
    # Create checker and analyze
    checker = H5ADDataChecker()
    result = checker.check_data_type(file_path)
    
    # Custom analysis
    print("\n" + "="*60)
    print("CUSTOM ANALYSIS")
    print("="*60)
    
    if 'error' in result:
        print(f"Error analyzing file: {result['error']}")
        return
    
    # Check data types
    print(f"File: {result['file_name']}")
    print(f"Dataset shape: {result['shape']}")
    
    # Raw data analysis
    if result['has_raw']:
        print("✓ Raw data available")
        raw_info = result['raw_data_info']
        print(f"  Raw data shape: {raw_info['shape']}")
        print(f"  Raw data type: {raw_info['type']}")
        if 'mean_value' in raw_info:
            print(f"  Raw data mean: {raw_info['mean_value']:.2f}")
    else:
        print("✗ No raw data available")
    
    # Normalized data analysis
    if result['has_normalized']:
        print("✓ Normalized data available")
        if 'X' in result['normalized_data_info']:
            x_info = result['normalized_data_info']['X']
            if x_info['is_normalized']:
                print("  ✓ X matrix is normalized")
            else:
                print("  ⚠ X matrix appears to be raw data")
    else:
        print("✗ No normalized data available")
    
    # Data layers analysis
    if result['data_layers']:
        print(f"✓ {len(result['data_layers'])} data layers available:")
        for layer in result['data_layers']:
            print(f"  - {layer}")
    
    # Recommendations
    print("\nRecommendations:")
    for rec in result['recommendations']:
        print(f"  {rec}")


def main():
    """Main function to run examples."""
    print("H5AD Data Checker Examples")
    print("=" * 50)
    
    # Check if any h5ad files exist in current directory
    h5ad_files = list(Path('.').glob('*.h5ad'))
    
    if h5ad_files:
        print(f"Found {len(h5ad_files)} h5ad files in current directory:")
        for f in h5ad_files:
            print(f"  - {f}")
        print()
        
        # Run example with first file
        example_file = str(h5ad_files[0])
        print(f"Running example with: {example_file}")
        
        checker = H5ADDataChecker()
        result = checker.check_data_type(example_file)
        
    else:
        print("No h5ad files found in current directory.")
        print("Please place some h5ad files in this directory or modify the script to use your file paths.")
        print("\nExample usage:")
        print("  python check_h5ad_data.py your_file.h5ad")
        print("  python check_h5ad_data.py file1.h5ad file2.h5ad --output results.json")


if __name__ == "__main__":
    main()
