#!/usr/bin/env python3
"""
Test script for the cellxgene downloader.

This script tests the functionality without downloading large files.
"""

import os
import sys
import tempfile
import yaml
from pathlib import Path
import unittest
from unittest.mock import patch, MagicMock

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from source.download_cxg.cellxgene_downloader import CellxgeneCensusDownloader


class TestCellxgeneCensusDownloader(unittest.TestCase):
    """Test cases for CellxgeneCensusDownloader class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.config_data = {
            'output_directory': self.temp_dir,
            'datasets': [
                {
                    'name': 'test_dataset',
                    'url': 'https://example.com/test.h5ad',
                    'description': 'Test dataset'
                }
            ],
            'download_settings': {
                'timeout': 30,
                'retry_attempts': 2,
                'chunk_size': 1024,
                'verify_ssl': True
            },
            'logging': {
                'level': 'INFO',
                'file': 'test.log',
                'console_output': False
            }
        }
        
        # Create temporary config file
        self.config_path = os.path.join(self.temp_dir, 'test_config.yaml')
        with open(self.config_path, 'w') as f:
            yaml.dump(self.config_data, f)
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_load_config(self):
        """Test configuration loading."""
        downloader = CellxgeneCensusDownloader(self.config_path)
        self.assertEqual(downloader.config, self.config_data)
    
    def test_load_nonexistent_config(self):
        """Test loading non-existent configuration file."""
        with self.assertRaises(SystemExit):
            CellxgeneCensusDownloader('nonexistent.yaml')
    
    def test_create_output_directory(self):
        """Test output directory creation."""
        downloader = CellxgeneCensusDownloader(self.config_path)
        output_dir = downloader._create_output_directory()
        self.assertTrue(output_dir.exists())
        self.assertEqual(str(output_dir), self.temp_dir)
    
    def test_get_filename_from_url(self):
        """Test filename extraction from URL."""
        downloader = CellxgeneCensusDownloader(self.config_path)
        
        # Test with .h5ad extension
        url = "https://example.com/dataset.h5ad"
        filename = downloader._get_filename_from_url(url)
        self.assertEqual(filename, "dataset.h5ad")
        
        # Test without .h5ad extension
        url = "https://example.com/dataset"
        filename = downloader._get_filename_from_url(url)
        self.assertEqual(filename, "dataset.h5ad")
    
    @patch('requests.Session')
    def test_validate_datasets_success(self, mock_session):
        """Test dataset validation with successful responses."""
        # Mock successful HEAD response
        mock_response = MagicMock()
        mock_response.raise_for_status.return_value = None
        mock_session.return_value.head.return_value = mock_response
        
        downloader = CellxgeneCensusDownloader(self.config_path)
        failed_datasets = downloader.validate_datasets()
        
        self.assertEqual(failed_datasets, [])
    
    @patch('requests.Session')
    def test_validate_datasets_failure(self, mock_session):
        """Test dataset validation with failed responses."""
        # Mock failed HEAD response
        mock_session.return_value.head.side_effect = Exception("Connection error")
        
        downloader = CellxgeneCensusDownloader(self.config_path)
        failed_datasets = downloader.validate_datasets()
        
        self.assertEqual(failed_datasets, ['test_dataset'])
    
    def test_list_available_datasets(self):
        """Test listing available datasets."""
        downloader = CellxgeneCensusDownloader(self.config_path)
        
        # Capture stdout to test output
        from io import StringIO
        import sys
        
        captured_output = StringIO()
        sys.stdout = captured_output
        
        try:
            downloader.list_available_datasets()
            output = captured_output.getvalue()
            
            # Check if output contains expected information
            self.assertIn("CONFIGURED DATASETS", output)
            self.assertIn("test_dataset", output)
            self.assertIn("https://example.com/test.h5ad", output)
            self.assertIn("Test dataset", output)
        finally:
            sys.stdout = sys.__stdout__


def run_quick_test():
    """Run a quick functionality test."""
    print("Running quick functionality test...")
    
    # Create a minimal test configuration
    test_config = {
        'output_directory': './test_downloads',
        'datasets': [
            {
                'name': 'test_pbmc3k',
                'url': 'https://cellxgene-example-data.czi.technology/pbmc3k.h5ad',
                'description': 'Test PBMC dataset'
            }
        ],
        'download_settings': {
            'timeout': 10,
            'retry_attempts': 1,
            'chunk_size': 1024,
            'verify_ssl': True
        },
        'logging': {
            'level': 'INFO',
            'console_output': True
        }
    }
    
    # Write test config
    with open('test_config.yaml', 'w') as f:
        yaml.dump(test_config, f)
    
    try:
        # Test configuration loading
        downloader = CellxgeneCensusDownloader('test_config.yaml')
        print("✓ Configuration loading: PASSED")
        
        # Test output directory creation
        output_dir = downloader._create_output_directory()
        print(f"✓ Output directory creation: PASSED ({output_dir})")
        
        # Test filename extraction
        filename = downloader._get_filename_from_url("https://example.com/test.h5ad")
        assert filename == "test.h5ad"
        print("✓ Filename extraction: PASSED")
        
        # Test dataset listing
        downloader.list_available_datasets()
        print("✓ Dataset listing: PASSED")
        
        print("\nAll basic functionality tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists('test_config.yaml'):
            os.remove('test_config.yaml')
        if os.path.exists('./test_downloads'):
            import shutil
            shutil.rmtree('./test_downloads')
    
    return True


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "quick":
        success = run_quick_test()
        sys.exit(0 if success else 1)
    else:
        # Run unit tests
        unittest.main(argv=[''], exit=False, verbosity=2)
