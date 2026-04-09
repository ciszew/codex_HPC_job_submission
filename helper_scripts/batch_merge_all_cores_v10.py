#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================================================
Batch 4-Compartment Quant Merge (v10 - Final)
================================================================================
Purpose:
  This script uses a robust, index-based join method to merge the four
  quantification files. This version fixes an UnboundLocalError from v9.

Usage:
  python batch_merge_all_cores_v10.py --project-dir /path/to/your/project
"""

import pandas as pd
import argparse
from skimage import io as skio
from pathlib import Path
import sys
import re

def process_single_core(core_id, quant_dir, seg_dir, output_dir):
    """Loads, joins, and saves the 4 quant files for a single core."""
    wc_quant_path = quant_dir / f"{core_id}--mesmer_cell.csv"
    nuc_quant_path = quant_dir / f"{core_id}--mesmer_nuclear.csv"
    cyto_quant_path = quant_dir / f"{core_id}--mesmer_cytoplasm.csv"
    mem_quant_path = quant_dir / f"{core_id}--mesmer_membrane.csv"
    wc_mask_path = seg_dir / f"mesmer-{core_id}/cell.tif"

    required_files = [wc_quant_path, nuc_quant_path, wc_mask_path]
    if not all(f.exists() for f in required_files):
        print(f"  [!] Skipping core {core_id}: Missing one or more required base files (cell, nuclear).")
        return False

    df_wc = pd.read_csv(wc_quant_path)
    df_nuc = pd.read_csv(nuc_quant_path)
    df_cyto = pd.read_csv(cyto_quant_path) if cyto_quant_path.exists() else pd.DataFrame({'CellID': []})
    df_mem = pd.read_csv(mem_quant_path) if mem_quant_path.exists() else pd.DataFrame({'CellID': []})
    
    wc_mask = skio.imread(wc_mask_path)
    
    nuc_centroids = df_nuc[['CellID', 'X_centroid', 'Y_centroid']].copy()
    nuc_centroids['x_int'] = nuc_centroids['X_centroid'].round().astype(int)
    nuc_centroids['y_int'] = nuc_centroids['Y_centroid'].round().astype(int)
    height, width = wc_mask.shape
    nuc_centroids = nuc_centroids[
        (nuc_centroids['y_int'] >= 0) & (nuc_centroids['y_int'] < height) &
        (nuc_centroids['x_int'] >= 0) & (nuc_centroids['x_int'] < width)
    ]
    wc_ids = wc_mask[nuc_centroids['y_int'].values, nuc_centroids['x_int'].values]
    nuc_centroids['wc_CellID'] = wc_ids
    id_map = nuc_centroids[nuc_centroids['wc_CellID'] > 0][['CellID', 'wc_CellID']].rename(columns={'CellID': 'nuc_CellID'})
    
    id_map = id_map.drop_duplicates(subset=['wc_CellID'], keep=False)
    id_map = id_map.drop_duplicates(subset=['nuc_CellID'], keep=False)
    
    df_wc = df_wc.set_index('CellID')
    df_cyto = df_cyto.set_index('CellID')
    df_mem = df_mem.set_index('CellID')
    
    df_nuc = df_nuc.rename(columns={'CellID': 'nuc_CellID'})
    df_nuc = pd.merge(df_nuc, id_map, on='nuc_CellID').set_index('wc_CellID')
    
    df_wc = df_wc.add_suffix('_wc')
    df_nuc = df_nuc.add_suffix('_nuc')
    df_cyto = df_cyto.add_suffix('_cyto')
    df_mem = df_mem.add_suffix('_mem')

    master_df = df_wc.join([df_nuc, df_cyto, df_mem], how='left')
    master_df = master_df.reset_index().rename(columns={'index': 'CellID'})

    # --- FIX: Define output_path BEFORE using it ---
    output_path = output_dir / f"{core_id}_master_quant.csv"
    master_df.to_csv(output_path, index=False)
    print(f"  -> Master file saved to: {output_path} ({len(master_df)} cells)")
    
    return True

def main():
    """Main function to parse arguments and loop through all cores."""
    parser = argparse.ArgumentParser(description='Batch merge 4-compartment quantification files for all cores.')
    parser.add_argument('--project-dir', required=True, help='Path to the main mcmicro project directory.')
    args = parser.parse_args()

    project_dir = Path(args.project_dir)
    quant_dir = project_dir / 'quantification'
    seg_dir = project_dir / 'segmentation'
    output_dir = project_dir / 'quantification_merged'

    if not quant_dir.is_dir() or not seg_dir.is_dir():
        print(f"Error: Could not find 'quantification' and/or 'segmentation' subdirectories in '{project_dir}'", file=sys.stderr)
        sys.exit(1)
        
    core_ids = set()
    pattern = re.compile(r"(\d+)--mesmer_cell\.csv")
    for f in quant_dir.glob('*--mesmer_cell.csv'):
        match = pattern.match(f.name)
        if match:
            core_ids.add(match.group(1))

    if not core_ids:
        print(f"Error: No quantification files found in '{quant_dir}'.", file=sys.stderr)
        sys.exit(1)

    print(f"Detected {len(core_ids)} cores to process in '{project_dir}'.")
    output_dir.mkdir(exist_ok=True)
    
    processed_count = 0
    failed_cores = []
    for core_id in sorted(list(core_ids), key=int):
        print(f"--- Processing Core {core_id} ---")
        try:
            if process_single_core(core_id, quant_dir, seg_dir, output_dir):
                processed_count += 1
            else:
                failed_cores.append(core_id)
        except Exception as e:
            print(f"  [!] An unexpected error occurred while processing core {core_id}: {e}")
            failed_cores.append(core_id)

    print("-" * 40)
    print(f"Batch processing complete. ✨")
    print(f"Successfully processed {processed_count}/{len(core_ids)} cores.")
    if failed_cores:
        print(f"The following cores failed or were skipped: {', '.join(failed_cores)}")

if __name__ == "__main__":
    main()