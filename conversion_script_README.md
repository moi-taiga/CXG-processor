Gene ID Conversion Script for Cellxgene H5AD Files
Purpose
Converts gene identifiers in Cellxgene .h5ad files to Ensembl gene IDs with automatic detection and robust fallback handling.
Logic Flow
1. Priority-Based Detection

Input H5AD → Check for existing Ensembl IDs → Use directly (100% success)
                ↓ (if not found)
            Map from gene symbols using TSV → Report success rate
                ↓ (if incomplete)
            Optional online lookup via mygene

2. Column Search Order
Ensembl IDs: id, _index, ensembl_gene_id, Accession, Gene stable ID
Gene Symbols: gene_short_name, feature_name, gene_symbol, gene_name

3. Automatic Format Handling
Categorical data (groups with categories/codes)
Direct datasets
Mixed storage formats

4. Usage
Arguments
Argument	Description	Default
--input-h5ad	Input H5AD file path	Required
--output-h5ad	Output H5AD file path	Required
--species	Species for online lookup	human
--set-var-names	Set var_names to	ensembl
--mapping-tsv	Symbol→Ensembl mapping file	Optional
--sanity-check	Run validation checks	False

5. Output
H5AD file with ensembl_gene_id column added
Optional var_names update to Ensembl IDs
Comprehensive conversion statistics
Input/output validation checks

6. Edge Cases Handled
Missing Python modules → H5py fallback
Different data storage formats → Automatic detection
Column name variations → Flexible matching
Species mismatches → Clear error reporting
Partial mappings → Detailed success statistics

7. Success Scenarios
Existing Ensembl IDs: 100% success, no mapping needed
Symbol mapping: Success rate depends on mapping file completeness
Online lookup: Fallback when mapping files are incomplete