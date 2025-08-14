#!/usr/bin/env python3
"""
Cellxgene Census H5AD Downloader

A Python script to download h5ad files using the cellxgene-census package
with configuration from a YAML file.

Author: AI Assistant
Date: 2024
"""

import os
import sys
import yaml
import logging
import requests
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlparse
import time
from tqdm import tqdm


class CellxgeneCensusDownloader:
    """Main class for downloading h5ad files using cellxgene-census package."""
    
    def __init__(self, config_path: str = "config.yaml"):
        """
        Initialize the downloader with configuration.
        
        Args:
            config_path (str): Path to the YAML configuration file
        """
        self.config_path = config_path
        self.config = self._load_config()
        self._setup_logging()
        self.session = requests.Session()
        
        # Log that configuration was loaded successfully
        self.logger.info(f"Configuration loaded from {self.config_path}")
        
        # Set custom headers if specified in config
        if 'headers' in self.config:
            self.session.headers.update(self.config['headers'])
    
    def _load_config(self) -> Dict:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r') as file:
                config = yaml.safe_load(file)
            return config
        except FileNotFoundError:
            print(f"Configuration file {self.config_path} not found")
            sys.exit(1)
        except yaml.YAMLError as e:
            print(f"Error parsing YAML configuration: {e}")
            sys.exit(1)
    
    def _setup_logging(self):
        """Setup logging configuration."""
        log_config = self.config.get('logging', {})
        log_level = getattr(logging, log_config.get('level', 'INFO'))
        
        # Create logger
        self.logger = logging.getLogger('cellxgene_census_downloader')
        self.logger.setLevel(log_level)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # File handler
        if 'file' in log_config:
            file_handler = logging.FileHandler(log_config['file'])
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
        
        # Console handler
        if log_config.get('console_output', True):
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)
    
    def _create_output_directory(self) -> Path:
        """Create output directory if it doesn't exist."""
        output_dir = Path(self.config.get('output_directory', './downloads'))
        output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output directory: {output_dir.absolute()}")
        return output_dir
    
    def _get_filename_from_url(self, url: str) -> str:
        """Extract filename from URL."""
        parsed_url = urlparse(url)
        filename = os.path.basename(parsed_url.path)
        if not filename.endswith('.h5ad'):
            filename += '.h5ad'
        return filename
    
    def _download_file(self, url: str, output_path: Path, dataset_name: str) -> bool:
        """
        Download a single file with progress bar and retry logic.
        
        Args:
            url (str): URL to download from
            output_path (Path): Local path to save the file
            dataset_name (str): Name of the dataset for logging
            
        Returns:
            bool: True if download successful, False otherwise
        """
        download_settings = self.config.get('download_settings', {})
        timeout = download_settings.get('timeout', 300)
        retry_attempts = download_settings.get('retry_attempts', 3)
        chunk_size = download_settings.get('chunk_size', 8192)
        verify_ssl = download_settings.get('verify_ssl', True)
        
        for attempt in range(retry_attempts):
            try:
                self.logger.info(f"Downloading {dataset_name} (attempt {attempt + 1}/{retry_attempts})")
                
                # Get file size for progress bar
                response = self.session.head(url, timeout=timeout, verify=verify_ssl)
                response.raise_for_status()
                total_size = int(response.headers.get('content-length', 0))
                
                # Download with progress bar
                response = self.session.get(url, stream=True, timeout=timeout, verify=verify_ssl)
                response.raise_for_status()
                
                with open(output_path, 'wb') as file:
                    with tqdm(
                        total=total_size,
                        unit='B',
                        unit_scale=True,
                        desc=f"Downloading {dataset_name}",
                        disable=not self.config.get('logging', {}).get('console_output', True)
                    ) as pbar:
                        for chunk in response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                file.write(chunk)
                                pbar.update(len(chunk))
                
                self.logger.info(f"Successfully downloaded {dataset_name} to {output_path}")
                return True
                
            except requests.exceptions.RequestException as e:
                self.logger.warning(f"Download attempt {attempt + 1} failed for {dataset_name}: {e}")
                if attempt < retry_attempts - 1:
                    wait_time = 2 ** attempt  # Exponential backoff
                    self.logger.info(f"Waiting {wait_time} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    self.logger.error(f"Failed to download {dataset_name} after {retry_attempts} attempts")
                    return False
        
        return False
    
    def download_datasets(self) -> Dict[str, bool]:
        """
        Download all datasets specified in the configuration.
        
        Returns:
            Dict[str, bool]: Dictionary mapping dataset names to download success status
        """
        output_dir = self._create_output_directory()
        datasets = self.config.get('datasets', [])
        results = {}
        
        self.logger.info(f"Starting download of {len(datasets)} datasets")
        
        for dataset in datasets:
            name = dataset['name']
            url = dataset['url']
            description = dataset.get('description', 'No description available')
            
            self.logger.info(f"Processing dataset: {name} - {description}")
            
            # Create filename
            filename = self._get_filename_from_url(url)
            output_path = output_dir / filename
            
            # Check if file already exists
            if output_path.exists():
                self.logger.info(f"File {output_path} already exists, skipping download")
                results[name] = True
                continue
            
            # Download the file
            success = self._download_file(url, output_path, name)
            results[name] = success
        
        return results
    
    def query_census_datasets(self) -> Dict[str, bool]:
        """
        Query datasets from cellxgene-census and save as h5ad files.
        
        Returns:
            Dict[str, bool]: Dictionary mapping dataset names to query success status
        """
        census_queries = self.config.get('census_queries', [])
        specific_datasets = self.config.get('specific_datasets', [])
        results = {}
        
        total_queries = len(census_queries) + len(specific_datasets)
        if total_queries == 0:
            self.logger.info("No census queries configured")
            return results
        
        self.logger.info(f"Starting census queries for {total_queries} datasets")
        
        # Process regular census queries
        for query in census_queries:
            name = query['name']
            description = query.get('description', 'No description available')
            query_params = query.get('query_params', {})
            
            self.logger.info(f"Processing census query: {name} - {description}")
            
            # Check if file already exists
            output_dir = self._create_output_directory()
            output_path = output_dir / f"{name}.h5ad"
            
            if output_path.exists():
                self.logger.info(f"File {output_path} already exists, skipping query")
                results[name] = True
                continue
            
            # Query the census data
            success = self.query_census_data(name, query_params)
            results[name] = success
        
        # Process specific dataset downloads
        for dataset in specific_datasets:
            name = dataset['name']
            description = dataset.get('description', 'No description available')
            dataset_id = dataset.get('dataset_id')
            organism = dataset.get('organism', 'Homo sapiens')
            max_cells = dataset.get('max_cells', 1000)
            
            self.logger.info(f"Processing specific dataset: {name} - {description}")
            self.logger.info(f"Dataset ID: {dataset_id}")
            
            # Check if file already exists
            output_dir = self._create_output_directory()
            output_path = output_dir / f"{name}.h5ad"
            
            if output_path.exists():
                self.logger.info(f"File {output_path} already exists, skipping download")
                results[name] = True
                continue
            
            # Download the specific dataset
            success = self.download_specific_dataset(name, dataset_id, organism)
            results[name] = success
        
        return results
    
    def query_census_data(self, dataset_name: str, query_params: Dict = None) -> bool:
        """
        Query data from cellxgene-census and save as h5ad file.
        
        Args:
            dataset_name (str): Name of the dataset to query
            query_params (Dict): Query parameters for census data
            
        Returns:
            bool: True if query and save successful, False otherwise
        """
        try:
            import cellxgene_census
            
            self.logger.info(f"Querying census data for dataset: {dataset_name}")
            
            # Default query parameters if none provided
            if query_params is None:
                query_params = {
                    'organism': 'Homo sapiens',
                    'obs_value_filter': None,
                    'max_cells': 1000
                }
            
            # Extract parameters
            organism = query_params.get('organism', 'Homo sapiens')
            obs_value_filter = query_params.get('obs_value_filter')
            
            # Use the correct API pattern with context manager
            CENSUS_VERSION = "2023-12-15"
            with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
                # Build the query parameters
                query_kwargs = {
                    'census': census,
                    'organism': organism
                }
                
                # Add optional filters
                if obs_value_filter:
                    query_kwargs['obs_value_filter'] = obs_value_filter
                
                # Query the census data
                adata = cellxgene_census.get_anndata(**query_kwargs)
            
            # Save as h5ad file
            output_dir = self._create_output_directory()
            output_path = output_dir / f"{dataset_name}.h5ad"
            
            adata.write_h5ad(output_path)
            self.logger.info(f"Successfully saved census data to {output_path}")
            
            return True
            
        except ImportError:
            self.logger.error("cellxgene-census package not installed. Please install it first.")
            return False
        except Exception as e:
            self.logger.error(f"Error querying census data: {e}")
            return False
    
    def download_specific_dataset(self, dataset_name: str, dataset_id: str, organism: str = "Homo sapiens") -> bool:
        """
        Download a specific dataset by dataset_id from cellxgene-census.
        
        Args:
            dataset_name (str): Name for the dataset
            dataset_id (str): Specific dataset ID to download
            organism (str): Organism name (default: "Homo sapiens")
            
        Returns:
            bool: True if download successful, False otherwise
        """
        try:
            import cellxgene_census
            
            self.logger.info(f"Downloading specific dataset: {dataset_name}")
            self.logger.info(f"Dataset ID: {dataset_id}")
            self.logger.info(f"Organism: {organism}")
            
            # Use the correct API pattern with context manager
            CENSUS_VERSION = "2023-12-15"
            with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
                # Query the specific dataset
                adata = cellxgene_census.get_anndata(
                    census,
                    organism,
                    obs_value_filter=f'dataset_id=="{dataset_id}"'
                )
            
            # Save as h5ad file
            output_dir = self._create_output_directory()
            output_path = output_dir / f"{dataset_name}.h5ad"
            
            adata.write_h5ad(output_path)
            self.logger.info(f"Successfully downloaded dataset {dataset_id} to {output_path}")
            self.logger.info(f"Dataset shape: {adata.shape}")
            
            return True
            
        except ImportError:
            self.logger.error("cellxgene-census package not installed. Please install it first.")
            return False
        except Exception as e:
            self.logger.error(f"Error downloading specific dataset {dataset_id}: {e}")
            return False
    
    def validate_datasets(self) -> List[str]:
        """
        Validate that all dataset URLs are accessible.
        
        Returns:
            List[str]: List of dataset names that failed validation
        """
        failed_datasets = []
        datasets = self.config.get('datasets', [])
        timeout = self.config.get('download_settings', {}).get('timeout', 300)
        verify_ssl = self.config.get('download_settings', {}).get('verify_ssl', True)
        
        self.logger.info("Validating dataset URLs...")
        
        for dataset in datasets:
            name = dataset['name']
            url = dataset['url']
            
            try:
                response = self.session.head(url, timeout=timeout, verify=verify_ssl)
                response.raise_for_status()
                self.logger.info(f"✓ {name}: URL accessible")
            except requests.exceptions.RequestException as e:
                self.logger.error(f"✗ {name}: URL not accessible - {e}")
                failed_datasets.append(name)
        
        return failed_datasets
    
    def list_available_datasets(self):
        """List all datasets configured for download."""
        datasets = self.config.get('datasets', [])
        census_queries = self.config.get('census_queries', [])
        specific_datasets = self.config.get('specific_datasets', [])
        
        print(f"\n{'='*60}")
        print("CONFIGURED DATASETS")
        print(f"{'='*60}")
        
        # Direct URL downloads
        if datasets:
            print("DIRECT URL DOWNLOADS:")
            for i, dataset in enumerate(datasets, 1):
                print(f"{i}. {dataset['name']}")
                print(f"   URL: {dataset['url']}")
                print(f"   Description: {dataset.get('description', 'No description')}")
                print()
        
        # Census queries
        if census_queries:
            print("CENSUS QUERIES:")
            for i, query in enumerate(census_queries, len(datasets) + 1):
                print(f"{i}. {query['name']}")
                print(f"   Type: Census Query")
                print(f"   Description: {query.get('description', 'No description')}")
                print(f"   Parameters: {query.get('query_params', {})}")
                print()
        
        # Specific datasets
        if specific_datasets:
            print("SPECIFIC DATASETS (by dataset_id):")
            for i, dataset in enumerate(specific_datasets, len(datasets) + len(census_queries) + 1):
                print(f"{i}. {dataset['name']}")
                print(f"   Dataset ID: {dataset.get('dataset_id', 'Not specified')}")
                print(f"   Organism: {dataset.get('organism', 'Homo sapiens')}")
                print(f"   Description: {dataset.get('description', 'No description')}")
                print()
        
        total_datasets = len(datasets) + len(census_queries) + len(specific_datasets)
        print(f"Total datasets: {total_datasets}")
        print(f"Output directory: {self.config.get('output_directory', './downloads')}")


def main():
    """Main function to run the downloader."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Download h5ad files using cellxgene-census package with YAML configuration"
    )
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Path to configuration file (default: config.yaml)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only validate dataset URLs without downloading'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List configured datasets and exit'
    )
    parser.add_argument(
        '--census-only',
        action='store_true',
        help='Only query census data, skip direct downloads'
    )
    parser.add_argument(
        '--download-only',
        action='store_true',
        help='Only download direct URLs, skip census queries'
    )
    
    args = parser.parse_args()
    
    # Initialize downloader
    downloader = CellxgeneCensusDownloader(args.config)
    
    if args.list:
        downloader.list_available_datasets()
        return
    
    if args.validate_only:
        failed_datasets = downloader.validate_datasets()
        if failed_datasets:
            print(f"\nValidation failed for {len(failed_datasets)} datasets:")
            for dataset in failed_datasets:
                print(f"  - {dataset}")
            sys.exit(1)
        else:
            print("\nAll dataset URLs are accessible!")
        return
    
    # Handle different operation modes
    if args.census_only:
        # Only query census data
        print("Querying census data only...")
        results = downloader.query_census_datasets()
    elif args.download_only:
        # Only download direct URLs
        print("Downloading direct URLs only...")
        # Validate datasets before downloading
        failed_datasets = downloader.validate_datasets()
        if failed_datasets:
            print(f"\nWarning: {len(failed_datasets)} datasets failed validation:")
            for dataset in failed_datasets:
                print(f"  - {dataset}")
            
            response = input("\nContinue with download anyway? (y/N): ")
            if response.lower() != 'y':
                print("Download cancelled.")
                sys.exit(1)
        
        results = downloader.download_datasets()
    else:
        # Both census queries and direct downloads
        print("Processing both census queries and direct downloads...")
        
        # First, query census data
        census_results = downloader.query_census_datasets()
        
        # Then, download direct URLs
        failed_datasets = downloader.validate_datasets()
        if failed_datasets:
            print(f"\nWarning: {len(failed_datasets)} datasets failed validation:")
            for dataset in failed_datasets:
                print(f"  - {dataset}")
            
            response = input("\nContinue with download anyway? (y/N): ")
            if response.lower() != 'y':
                print("Download cancelled.")
                sys.exit(1)
        
        download_results = downloader.download_datasets()
        
        # Combine results
        results = {**census_results, **download_results}
    
    # Print summary
    print(f"\n{'='*60}")
    print("OPERATION SUMMARY")
    print(f"{'='*60}")
    
    successful = sum(1 for success in results.values() if success)
    total = len(results)
    
    for name, success in results.items():
        status = "✓ SUCCESS" if success else "✗ FAILED"
        print(f"{name}: {status}")
    
    print(f"\nTotal: {successful}/{total} datasets processed successfully")
    
    if successful < total:
        sys.exit(1)


if __name__ == "__main__":
    main()
