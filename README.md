# CXG Processor Pipeline

A pipeline for pulling data from cellXgene and undertaking some basic preprocessing.

## User Story

As a single-cell Bioinformatician, I want to automate my single-cell analysis, so that my workflow is faster.

## MVP

Automated pipeline to take CXG data and apply scanpy preprocessing.

## Stretch Goal

Continue the pipeline to GeneSwitches.

## Output Structure

Each run generates a directory containing:
- the data,
- a logfile,
- a config file,
- any graphical outputs (e.g. UMAPs, QC metrics),
- and the processed object saved as a `.h5ad` file.

-----------------------

# Cellxgene Census H5AD Downloader & Data Checker

A Python application for downloading and querying h5ad files using the [cellxgene-census](https://github.com/chanzuckerberg/cellxgene-census) package with YAML configuration support, plus tools to analyze h5ad files for raw vs normalized data. This tool is designed to efficiently access single-cell transcriptomics datasets from the CZ CELLxGENE Discover Census and download datasets from various sources.

## Features

### Download & Query Features
- **YAML Configuration**: Easy-to-use YAML configuration file for managing datasets
- **Progress Tracking**: Real-time download progress with progress bars
- **Retry Logic**: Automatic retry with exponential backoff for failed downloads
- **Validation**: URL validation before downloading
- **Logging**: Comprehensive logging to both file and console
- **Resume Support**: Skip already downloaded files
- **Flexible Output**: Configurable output directory and file naming

### Data Analysis Features
- **Raw vs Normalized Detection**: Automatically detect if h5ad files contain raw counts or normalized data
- **Statistical Analysis**: Comprehensive statistical analysis of data characteristics
- **Data Layer Analysis**: Examine additional data layers in h5ad files
- **Smart Recommendations**: Provide actionable recommendations based on data type
- **Batch Processing**: Analyze multiple h5ad files at once
- **JSON Output**: Save analysis results in structured JSON format

## Installation

### Prerequisites

- Python 3.10 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install cellxgene-census>=1.17.0 PyYAML>=6.0 requests>=2.28.0 pathlib2>=2.3.7 tqdm>=4.64.0
```

## Configuration

The application uses a YAML configuration file (`config.yaml`) to define datasets and settings:

### Basic Configuration Structure

```yaml
# Output directory for downloaded files
output_directory: "./downloads"

# List of datasets to download (direct URLs)
datasets:
  - name: "pbmc3k"
    url: "https://cellxgene-example-data.czi.technology/pbmc3k.h5ad"
    description: "Peripheral blood mononuclear cells (PBMC) dataset with 3k cells"
    
  - name: "pbmc10k"
    url: "https://cellxgene-example-data.czi.technology/pbmc10k.h5ad"
    description: "Peripheral blood mononuclear cells (PBMC) dataset with 10k cells"

# Census query configurations
census_queries:
  - name: "human_pbmc_sample"
    description: "Sample human PBMC data from census"
    query_params:
      organism: "Homo sapiens"
      obs_value_filter: 'tissue_general=="blood"'
      max_cells: 1000  # Limit for testing

# Specific dataset IDs to download
specific_datasets:
  - name: "dataset_8e47ed12"
    description: "Specific dataset with ID 8e47ed12-c658-4252-b126-381df8d52a3d"
    dataset_id: "8e47ed12-c658-4252-b126-381df8d52a3d"
    organism: "Homo sapiens"
    max_cells: 1000

# Download settings
download_settings:
  timeout: 300  # seconds
  retry_attempts: 3
  chunk_size: 8192  # bytes
  verify_ssl: true
  
# Logging configuration
logging:
  level: "INFO"  # DEBUG, INFO, WARNING, ERROR
  file: "download.log"
  console_output: true

# Optional: Custom headers for requests
headers:
  User-Agent: "cellxgene-census-downloader/1.0"
```

### Configuration Options

#### Direct URL Downloads
- `name`: Unique identifier for the dataset
- `url`: Direct URL to the h5ad file
- `description`: Optional description of the dataset

#### Census Queries
- `name`: Unique identifier for the query
- `description`: Optional description of the query
- `query_params`: Dictionary of query parameters
  - `organism`: Organism name (e.g., "Homo sapiens", "Mus musculus")
  - `obs_value_filter`: Filter string for observations

#### Specific Datasets (by dataset_id)
- `name`: Unique identifier for the dataset
- `description`: Optional description of the dataset
- `dataset_id`: Specific dataset ID from the census
- `organism`: Organism name (e.g., "Homo sapiens", "Mus musculus")

#### Download Settings
- `timeout`: Request timeout in seconds (default: 300)
- `retry_attempts`: Number of retry attempts for failed downloads (default: 3)
- `chunk_size`: Chunk size for streaming downloads in bytes (default: 8192)
- `verify_ssl`: Whether to verify SSL certificates (default: true)

#### Logging
- `level`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `file`: Log file path (optional)
- `console_output`: Whether to output logs to console (default: true)

## Usage

### Download & Query Usage

Download all datasets defined in the configuration:

```bash
python cellxgene_downloader.py
```

### Data Analysis Usage

Check if h5ad files contain raw or normalized data:

```bash
# Check a single file
python check_h5ad_data.py your_file.h5ad

# Check multiple files
python check_h5ad_data.py file1.h5ad file2.h5ad file3.h5ad

# Save results to JSON
python check_h5ad_data.py your_file.h5ad --output results.json
```

### Command Line Options

```bash
python cellxgene_downloader.py [OPTIONS]
```

#### Available Options

- `-c, --config`: Path to configuration file (default: config.yaml)
- `--validate-only`: Only validate dataset URLs without downloading
- `--list`: List configured datasets and exit
- `--census-only`: Only query census data, skip direct downloads
- `--download-only`: Only download direct URLs, skip census queries

### Examples

#### List Configured Datasets

```bash
python cellxgene_downloader.py --list
```

Output:
```
============================================================
CONFIGURED DATASETS
============================================================
DIRECT URL DOWNLOADS:
1. pbmc3k
   URL: https://cellxgene-example-data.czi.technology/pbmc3k.h5ad
   Description: Peripheral blood mononuclear cells (PBMC) dataset with 3k cells

2. pbmc10k
   URL: https://cellxgene-example-data.czi.technology/pbmc10k.h5ad
   Description: Peripheral blood mononuclear cells (PBMC) dataset with 10k cells

CENSUS QUERIES:
3. human_pbmc_sample
   Type: Census Query
   Description: Sample human PBMC data from census
   Parameters: {'organism': 'Homo sapiens', 'obs_value_filter': 'tissue_general=="blood"'}

SPECIFIC DATASETS (by dataset_id):
4. dataset_8e47ed12
   Dataset ID: 8e47ed12-c658-4252-b126-381df8d52a3d
   Organism: Homo sapiens
   Description: Specific dataset with ID 8e47ed12-c658-4252-b126-381df8d52a3d


Total datasets: 4
Output directory: ./downloads
```

#### Validate URLs Only

```bash
python cellxgene_downloader.py --validate-only
```

#### Query Census Data Only

```bash
python cellxgene_downloader.py --census-only
```

#### Download Direct URLs Only

```bash
python cellxgene_downloader.py --download-only
```

#### Use Custom Configuration File

```bash
python cellxgene_downloader.py -c my_config.yaml
```

#### Data Analysis Examples

```bash
# Analyze downloaded files
python check_h5ad_data.py downloads/*.h5ad

# Get detailed analysis with recommendations
python check_h5ad_data.py your_dataset.h5ad --verbose

# Save analysis results for later use
python check_h5ad_data.py *.h5ad --output analysis_results.json
```

## Output

### Downloaded Files

Files are downloaded to the specified output directory with their original filenames:

```
downloads/
├── pbmc3k.h5ad
├── pbmc10k.h5ad
└── mouse_brain.h5ad
```

### Analysis Results

The data checker provides detailed analysis output:

```
============================================================
ANALYZING: example.h5ad
============================================================
Loading h5ad file...
Dataset shape: (1000, 2000)
Cells: 1000, Genes: 2000

--- RAW DATA ANALYSIS ---
✓ Raw data found: (1000, 2000)
  Raw data range: 0.00 - 150.00
  Raw data mean: 2.45
  Zero fraction: 85.23%

--- NORMALIZED DATA ANALYSIS ---
Main data matrix (X): (1000, 2000)
  X data range: 0.00 - 5.01
  X data mean: 0.85
  Zero fraction: 85.23%
  ✓ X matrix appears to be normalized

--- RECOMMENDATIONS ---
✓ Both raw and normalized data available - ideal for analysis
  Use adata.raw.X for raw counts, adata.X for normalized data
```

### Log Files

If logging to file is enabled, logs are written to the specified log file:

```
2024-01-15 10:30:15 - cellxgene_downloader - INFO - Configuration loaded from config.yaml
2024-01-15 10:30:15 - cellxgene_downloader - INFO - Output directory: /path/to/downloads
2024-01-15 10:30:15 - cellxgene_downloader - INFO - Starting download of 3 datasets
2024-01-15 10:30:15 - cellxgene_downloader - INFO - Processing dataset: pbmc3k - Peripheral blood mononuclear cells (PBMC) dataset with 3k cells
2024-01-15 10:30:16 - cellxgene_downloader - INFO - Downloading pbmc3k (attempt 1/3)
2024-01-15 10:30:45 - cellxgene_downloader - INFO - Successfully downloaded pbmc3k to /path/to/downloads/pbmc3k.h5ad
```

## Integration with Cellxgene Census

The downloaded h5ad files can be used with various tools:

### Using cellxgene-census for data access:

```python
import cellxgene_census

# Query census data directly
CENSUS_VERSION = "2023-12-15"
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        "Homo sapiens",
        obs_value_filter='tissue_general=="blood"'
    )
```

### Using regular cellxgene for visualization:

```python
import cellxgene

# Launch cellxgene with downloaded file
cellxgene.launch("downloads/pbmc3k.h5ad")
```

Or from command line:

```bash
cellxgene launch downloads/pbmc3k.h5ad
```

## Data Analysis Integration

### Using the data checker in your analysis pipeline:

```python
from check_h5ad_data import H5ADDataChecker

# Check data type before analysis
checker = H5ADDataChecker()
result = checker.check_data_type("your_file.h5ad")

if result['has_raw'] and result['has_normalized']:
    print("Perfect! Both raw and normalized data available")
    # Use adata.raw.X for differential expression
    # Use adata.X for clustering and visualization
elif result['has_raw']:
    print("Only raw data - consider normalizing")
    import scanpy as sc
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
else:
    print("Check if data is properly loaded")
```

## Error Handling

The application includes comprehensive error handling:

### Download & Query Errors
- **Network Errors**: Automatic retry with exponential backoff
- **Invalid URLs**: Validation before download attempts
- **File System Errors**: Graceful handling of permission and space issues
- **Configuration Errors**: Clear error messages for malformed YAML

### Data Analysis Errors
- **File Format Errors**: Clear messages for invalid h5ad files
- **Memory Errors**: Graceful handling of large datasets
- **Missing Dependencies**: Helpful installation instructions
- **Data Type Errors**: Clear identification of data format issues

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Cellxgene Census](https://github.com/chanzuckerberg/cellxgene-census) - The main package for accessing CZ CELLxGENE Discover Census data
- [Chan Zuckerberg Initiative](https://chanzuckerberg.com/) - Developers of cellxgene-census
- [PyYAML](https://pyyaml.org/) - YAML parser for Python
- [Requests](https://requests.readthedocs.io/) - HTTP library for Python
- [tqdm](https://tqdm.github.io/) - Progress bar library

## Support

For issues and questions:

### Download & Query Issues
1. Check the [cellxgene-census documentation](https://chanzuckerberg.github.io/cellxgene-census/)
2. Review the log files for detailed error information
3. Ensure your configuration file is properly formatted
4. Verify that all URLs are accessible from your network

### Data Analysis Issues
1. Check the [anndata documentation](https://anndata.readthedocs.io/)
2. Verify h5ad file format and integrity
3. Ensure sufficient memory for large datasets
4. Check scanpy version compatibility

## Changelog

### Version 1.1.0
- Added h5ad data type checker functionality
- Raw vs normalized data detection
- Statistical analysis of data characteristics
- Data layer analysis
- Smart recommendations system
- JSON output support
- Batch processing capabilities

### Version 1.0.0
- Initial release
- YAML configuration support
- Progress tracking with tqdm
- Retry logic with exponential backoff
- URL validation
- Comprehensive logging
- Command-line interface




