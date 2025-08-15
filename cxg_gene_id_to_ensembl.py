#!/usr/bin/env python3
"""
Convert gene identifiers in a Cellxgene .h5ad file to Ensembl gene IDs.
Designed to work with Python 3.12.6 without conda dependencies.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np
import pandas as pd
import h5py

# Optional import guarded for online lookup
try:
    from mygene import MyGeneInfo  # type: ignore
    _HAS_MYG = True
except Exception:
    _HAS_MYG = False

PREFERRED_ENSEMBL_COLS = [
    "Gene stable ID", "gene_stable_id", "Accession", "ensembl_gene_id", 
    "ensembl_id", "feature_id", "gene_id", "Ensembl_ID", "id", "_index"
]

PREFERRED_SYMBOL_COLS = [
    "Gene name", "HGNC symbol", "feature_name", "gene_symbol", "gene", 
    "symbol", "hgnc_symbol", "gene_short_name", "gene_name"
]

def normalize_ensembl(ensembl_value: object) -> Optional[str]:
    """Return normalized Ensembl gene id (without version) or None."""
    if ensembl_value is None or (isinstance(ensembl_value, float) and np.isnan(ensembl_value)):
        return None
    text = str(ensembl_value)
    if not text or text == "nan":
        return None
    # Remove version suffix like ".12"
    return text.split(".")[0]

def build_symbol_to_ensembl_map(mapping_tsv: str) -> Dict[str, str]:
    """Build symbol->ensembl map from a TSV/CSV file."""
    print(f"📁 Loading mapping file: {mapping_tsv}")
    
    if not os.path.exists(mapping_tsv):
        raise FileNotFoundError(f"Mapping file not found: {mapping_tsv}")
    
    df = pd.read_csv(mapping_tsv, sep="\t" if mapping_tsv.endswith(".tsv") else ",")
    print(f"📊 Available columns: {list(df.columns)}")
    
    # Find symbol column - prioritize your file format
    symbol_col = None
    for cand in PREFERRED_SYMBOL_COLS:
        if cand in df.columns:
            symbol_col = cand
            break
    
    if symbol_col is None:
        available = list(df.columns)
        raise ValueError(f"Mapping file missing a symbol column. Available: {available}. Expected one of: {PREFERRED_SYMBOL_COLS}")
    
    # Find ensembl column - prioritize your file format
    ensembl_col = None
    for cand in PREFERRED_ENSEMBL_COLS:
        if cand in df.columns:
            ensembl_col = cand
            break
    
    if ensembl_col is None:
        available = list(df.columns)
        raise ValueError(f"Mapping file missing an Ensembl ID column. Available: {available}. Expected one of: {PREFERRED_ENSEMBL_COLS}")
    
    print(f"✅ Using symbol column: '{symbol_col}', ensembl column: '{ensembl_col}'")
    
    # Build mapping dictionary
    mapping: Dict[str, str] = {}
    mapped_count = 0
    
    for _, row in df.iterrows():
        ensg = normalize_ensembl(row.get(ensembl_col))
        if not ensg:
            continue
            
        # Get primary symbol
        primary = str(row.get(symbol_col, "")).strip()
        if primary and primary != "nan" and ensg:
            if primary not in mapping:  # First occurrence wins
                mapping[primary] = ensg
                mapped_count += 1
    
    print(f"📋 Created mapping with {mapped_count} symbol->ensembl pairs")
    return mapping

def _read_var_column_h5(var_group: h5py.Group, column_name: str) -> Optional[pd.Series]:
    """Read a column from var using h5py, handling categorical storage."""
    if column_name not in var_group:
        return None
    
    obj = var_group[column_name]
    
    # Handle categorical encoding
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        cats = obj["categories"][()]
        cats = [c.decode() if isinstance(c, (bytes, bytearray)) else str(c) for c in cats]
        codes = obj["codes"][()]
        values = [cats[int(i)] if 0 <= int(i) < len(cats) else "" for i in codes]
        return pd.Series(values, dtype=str)
    
    # Handle direct dataset
    if isinstance(obj, h5py.Dataset):
        arr = obj[()]
        values = [a.decode() if isinstance(a, (bytes, bytearray)) else str(a) for a in arr]
        return pd.Series(values, dtype=str)
    
    return None

def _write_string_column_h5(var_group: h5py.Group, column_name: str, values: pd.Series) -> None:
    """Write/overwrite a string dataset under var."""
    vals = values.fillna("").astype(str).tolist()
    maxlen = max(1, max(len(v) for v in vals))
    dt = f'S{maxlen}'
    data = np.array([v.encode('utf-8') for v in vals], dtype=dt)
    
    if column_name in var_group:
        del var_group[column_name]
    var_group.create_dataset(column_name, data=data, dtype=data.dtype)

def sanity_check_h5py(input_file: str, output_file: str = None):
    """Print first 10 gene IDs using h5py for sanity checking."""
    print(f"\n{'='*60}")
    print("INPUT FILE SANITY CHECK")
    print(f"{'='*60}")
    
    with h5py.File(input_file, 'r') as f:
        if 'var' not in f:
            print("❌ No 'var' group found in input file")
            return
        
        var = f['var']
        print(f"📊 Available var columns: {list(var.keys())}")
        
        # Check for existing Ensembl IDs first
        ensembl_found = False
        ensembl_column = None
        for col_name in PREFERRED_ENSEMBL_COLS:
            if col_name in var:
                ensembl_data = _read_var_column_h5(var, col_name)
                if ensembl_data is not None:
                    # Check if this looks like Ensembl IDs
                    sample_values = ensembl_data.head(5).astype(str)
                    if any('ENS' in str(val) for val in sample_values):
                        ensembl_found = True
                        ensembl_column = col_name
                        print(f"✅ Found existing Ensembl IDs in column: {col_name}")
                        print(f"🧬 Total genes: {len(ensembl_data)}")
                        print(f"🔍 First 10 Ensembl IDs:")
                        for i, gene in enumerate(ensembl_data[:10]):
                            print(f"  {i+1:2d}. {gene}")
                        break
        
        # If no Ensembl IDs found, look for gene names/symbols
        if not ensembl_found:
            gene_names = None
            for col_name in ['_index', 'index', 'gene_names', 'gene_short_name', 'feature_name']:
                if col_name in var:
                    gene_names = _read_var_column_h5(var, col_name)
                    if gene_names is not None:
                        print(f"📝 Found gene names in column: {col_name}")
                        print(f"🧬 Total genes: {len(gene_names)}")
                        print(f"🔍 First 10 gene names:")
                        for i, gene in enumerate(gene_names[:10]):
                            print(f"  {i+1:2d}. {gene}")
                        break
            
            if gene_names is None:
                print("❌ No Ensembl IDs or gene names found in input file")
            else:
                print("ℹ️  No existing Ensembl IDs found, will attempt mapping from gene names")
    
    # Check output file if it exists
    if output_file and os.path.exists(output_file):
        print(f"\n{'='*60}")
        print("OUTPUT FILE SANITY CHECK")
        print(f"{'='*60}")
        
        with h5py.File(output_file, 'r') as f:
            if 'var' in f:
                var = f['var']
                
                # Get output Ensembl IDs
                ensembl_ids = _read_var_column_h5(var, 'ensembl_gene_id')
                
                if ensembl_ids is not None:
                    print(f"🧬 Total genes: {len(ensembl_ids)}")
                    print(f"🔍 First 10 Ensembl IDs:")
                    
                    for i in range(min(10, len(ensembl_ids))):
                        ensg = ensembl_ids.iloc[i] if i < len(ensembl_ids) else "N/A"
                        status = "✅" if ensg != "N/A" and ensg != "" and ensg != "nan" else "❌"
                        print(f"  {i+1:2d}. {ensg} {status}")
                    
                    # Calculate mapping success
                    mapped = ensembl_ids.notna().sum()
                    total = len(ensembl_ids)
                    print(f"\n📈 Mapping success: {mapped}/{total} ({mapped/total*100:.1f}%)")
                    
                    if mapped == total:
                        print("🎉 All genes successfully converted to Ensembl IDs!")
                    elif mapped > 0:
                        print(f"⚠️  {total-mapped} genes failed to convert")
                    else:
                        print("❌ No genes were successfully converted")
                else:
                    print("❌ Could not read Ensembl IDs from output file")
    
    print(f"{'='*60}\n")

def convert_with_h5py_complete(input_h5ad: str, output_h5ad: str, mapping_tsv: Optional[str], 
                              set_var_names: str, species: str, allow_online: bool, 
                              sanity_check: bool = False) -> Tuple[int, int]:
    """Complete h5py-based conversion that avoids anndata entirely."""
    
    print("🔄 Using h5py-only conversion path...")
    
    if sanity_check:
        sanity_check_h5py(input_h5ad)
    
    # Copy entire file first
    with h5py.File(input_h5ad, 'r') as fin:
        with h5py.File(output_h5ad, 'w') as fout:
            # Copy all groups and datasets
            for key in fin.keys():
                fin.copy(fin[key], fout, name=key)
            # Copy root attributes
            for akey, aval in fin.attrs.items():
                fout.attrs[akey] = aval

    # Re-open output for modification
    with h5py.File(output_h5ad, 'a') as f:
        if 'var' not in f:
            raise RuntimeError("No 'var' group in h5ad file")
        
        var = f['var']
        
        # 1) Check if ensembl-like column exists
        ensembl_series = None
        for cand in PREFERRED_ENSEMBL_COLS:
            ensembl_series = _read_var_column_h5(var, cand)
            if ensembl_series is not None:
                print(f"✅ Found existing Ensembl column: {cand}")
                ens = ensembl_series.map(normalize_ensembl)
                _write_string_column_h5(var, 'ensembl_gene_id', ens)
                mapped = int(ens.notna().sum())
                total = int(len(ens))
                
                if set_var_names == 'ensembl' and '_index' in var:
                    _write_string_column_h5(var, '_index', ens.fillna("UNMAPPED"))
                
                if sanity_check:
                    sanity_check_h5py(input_h5ad, output_h5ad)
                
                return mapped, total
        
        # 2) Map from symbols
        sym_series = None
        used_column = None
        for cand in PREFERRED_SYMBOL_COLS:
            sym_series = _read_var_column_h5(var, cand)
            if sym_series is not None:
                used_column = cand
                print(f"✅ Using symbol column: {cand}")
                break
        
        # Fallback to _index if no symbol column found
        if sym_series is None:
            sym_series = _read_var_column_h5(var, '_index')
            if sym_series is not None:
                used_column = '_index'
                print("✅ Using _index as symbol column")
        
        if sym_series is None:
            raise RuntimeError("Could not find any symbol column in var to map from")
        
        # Build mapping
        mapping: Dict[str, str] = {}
        if mapping_tsv and os.path.exists(mapping_tsv):
            try:
                mapping.update(build_symbol_to_ensembl_map(mapping_tsv))
            except Exception as e:
                print(f"❌ Error loading mapping file: {e}")
        
        # Online resolution if enabled
        if allow_online and len(mapping) < sym_series.nunique() and _HAS_MYG:
            try:
                print("🌐 Attempting online gene ID resolution...")
                mg = MyGeneInfo()
                unmapped_symbols = sym_series[~sym_series.isin(mapping.keys())].unique()
                res = mg.querymany(unmapped_symbols.tolist(), species=species, 
                                 scopes=["symbol"], fields=["ensembl.gene"], 
                                 as_dataframe=True, returnall=False, verbose=False)
                
                if isinstance(res, pd.DataFrame):
                    for sym, row in res.iterrows():
                        val = row.get("ensembl.gene")
                        ensg = None
                        if isinstance(val, list) and val:
                            ensg = normalize_ensembl(val[0])
                        elif isinstance(val, str):
                            ensg = normalize_ensembl(val)
                        if ensg:
                            mapping[str(sym)] = ensg
                    print(f"🌐 Online resolution found {len(res)} additional mappings")
            except Exception as e:
                print(f"⚠️  Online resolution failed: {e}")
        
        # Apply mapping
        ens = sym_series.map(mapping).map(normalize_ensembl)
        _write_string_column_h5(var, 'ensembl_gene_id', ens)
        mapped = int(ens.notna().sum())
        total = int(len(ens))
        
        print(f"📊 Mapping results: {mapped}/{total} ({mapped/total*100:.1f}%)")
        
        # Update var_names if requested
        if set_var_names == 'ensembl' and '_index' in var:
            _write_string_column_h5(var, '_index', ens.fillna("UNMAPPED"))
            print("✅ Set var_names to Ensembl IDs")
        elif set_var_names == 'symbol' and used_column != '_index' and '_index' in var:
            _write_string_column_h5(var, '_index', sym_series)
            print(f"✅ Set var_names to symbols from {used_column}")
        
        if sanity_check:
            sanity_check_h5py(input_h5ad, output_h5ad)
        
        return mapped, total

def main():
    ap = argparse.ArgumentParser(description="Convert gene identifiers in Cellxgene .h5ad to Ensembl gene IDs (Python 3.12.6 compatible)")
    ap.add_argument("--input-h5ad", required=True, help="Input .h5ad from Cellxgene")
    ap.add_argument("--output-h5ad", required=True, help="Output .h5ad path with Ensembl IDs added")
    ap.add_argument("--mapping-tsv", help="TSV/CSV with symbol->ensembl mapping")
    ap.add_argument("--species", default="human", choices=["human", "mouse", "zebrafish"], 
                   help="Species for online lookup (if enabled)")
    ap.add_argument("--allow-online", action="store_true", help="Allow online lookup via mygene")
    ap.add_argument("--set-var-names", choices=["none", "ensembl", "symbol"], default="ensembl",
                   help="What to set .var_names to in output (default: ensembl)")
    ap.add_argument("--report-csv", help="Optional path to write a mapping report CSV")
    ap.add_argument("--sanity-check", action="store_true", help="Print first 10 gene IDs before and after conversion")
    
    args = ap.parse_args()

    if not os.path.exists(args.input_h5ad):
        print(f"❌ Error: input not found: {args.input_h5ad}")
        sys.exit(1)

    print("🐍 Running on Python 3.12.6 without conda dependencies")
    print("📦 Using h5py-only conversion path to avoid anndata/bz2 issues")
    
    try:
        # Use h5py-only conversion
        mapped, total = convert_with_h5py_complete(
            input_h5ad=args.input_h5ad,
            output_h5ad=args.output_h5ad,
            mapping_tsv=args.mapping_tsv,
            set_var_names=args.set_var_names,
            species=args.species,
            allow_online=args.allow_online,
            sanity_check=args.sanity_check
        )
        
        # Write report if requested
        if args.report_csv:
            report = pd.DataFrame({
                "total_genes": [int(total)],
                "mapped_genes": [int(mapped)],
                "percent_mapped": [float(mapped/total*100.0)],
                "failed_genes": [int(total-mapped)]
            })
            report.to_csv(args.report_csv, index=False)
            print(f"📄 Report written: {args.report_csv}")

        # Final summary
        print(f"\n{'='*60}")
        print("🎉 CONVERSION COMPLETED SUCCESSFULLY")
        print(f"{'='*60}")
        print(f"📁 Input file: {args.input_h5ad}")
        print(f"📁 Output file: {args.output_h5ad}")
        print(f"🧬 Total genes: {total}")
        print(f"✅ Successfully mapped: {mapped} ({mapped/total*100:.1f}%)")
        print(f"❌ Failed mappings: {total-mapped} ({(total-mapped)/total*100:.1f}%)")
        print(f"🔬 Species: {args.species}")
        print(f"🐍 Python version: 3.12.6 (no conda)")
        print(f"{'='*60}")
        
    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        print(f"📋 Make sure your mapping file has the correct columns:")
        print(f"   - Symbol columns: {PREFERRED_SYMBOL_COLS}")
        print(f"   - Ensembl columns: {PREFERRED_ENSEMBL_COLS}")
        sys.exit(1)

if __name__ == "__main__":
    main()
