#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Standalone Ashlar Submission Wrapper for Slurm HPC
================================================================================

Purpose:
  This script automates the staging, renaming, and submission of Ashlar image 
  stitching jobs on a Slurm cluster. It is engineered to be highly flexible, 
  supporting everything from massive multi-cycle multiplexed datasets down to 
  simple single-cycle or single-channel images.

Key Features:
  - Smart Directory Scanning: Automatically detects flat single-cycle directories, 
    nested single-cycle folders, or complex multi-cycle datasets.
  - Chronological Sorting: Sorts multi-cycle folders safely based on exact 
    timestamp regex (YYYY-MM-DD HH.MM.SS) embedded in the folder names.
  - Automated File Renaming: Strips arbitrary prefixes to create clean Ashlar-
    compliant file names (e.g., F1.ims, F2.ims) in-place.
  - Recursive XML Search: Hunts down the microscope metadata (.xml) regardless 
    of where it is nested inside the data folder tree.
  - Dynamic Resource Allocation: Estimates Slurm time limits based on total data size.

Usage:
  python3 run_ashlar_standalone2.py /path/to/raw_data/ [--align-channel INDEX] --align-channel 0 would be 1st channel in image stack, --align-channel -1 would be last

Positional Arguments:
  data_dir          Path to the target directory containing your image files 
                    or cycle subdirectories.

Optional Arguments:
  --align-channel   The 0-based index of the channel to use for image alignment. 
                    - Default: -1 (Automatically aligns using the LAST channel in the stack).
                    - Use 0 for single-channel datasets (e.g., pure DAPI).
                    - Pass a specific integer (e.g., --align-channel 2) to force 
                      alignment on a specific physical channel.

Expected Data Structures:
  1. Flat Directory: A single folder containing all .ims files and the .xml file.
  2. Nested/Multi-Cycle: A parent folder containing subdirectories (e.g., "Cycle_1", 
     "Cycle_2"), each containing .ims files, with the .xml file located anywhere inside.

Outputs:
  The script generates a `submit_ashlar.slurm` file in the target directory, 
  pauses to allow you to manually review/edit Slurm parameters (RAM, queues, MAX_SHIFT), 
  and then submits the job. The final stitched image will be saved to:
  /path/to/raw_data/registration/stitched.ome.tiff

================================================================================
"""
"""
Standalone Ashlar Submission Wrapper
Execute this script from within the directory containing your raw data folders.
for example: python run_ashlar_standalone.py /path/to/single_cycle_data --align-channel 0
"""

import os
import re
import sys
import subprocess
from datetime import datetime
import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser(description="Standalone Ashlar Orchestrator")
parser.add_argument('data_dir', type=str, help='Path to the raw data directory.')
parser.add_argument('--align-channel', type=int, default=-1, 
                    help='0-based index for the alignment channel (default: -1 for last channel. Use 0 for single-channel data).')

args = parser.parse_args()

# Extract command line arguments
user_align_channel = args.align_channel
source_data_dir = os.path.abspath(args.data_dir)

# --- Constants based on standard HDID environment ---
CONDA_ENV = "/gpfs/data/hdid-share/conda/ashlar_group"

def get_datetime_key(dir_name):
    """Parses directory name for date and time to create a sort key."""
    date_match = re.search(r"(\d{4}-\d{2}-\d{2})", dir_name)
    time_match = re.search(r"(\d{2}\.\d{2}\.\d{2})", dir_name)

    if date_match and time_match:
        try:
            datetime_str = f"{date_match.group(1)} {time_match.group(1)}"
            return datetime.strptime(datetime_str, "%Y-%m-%d %H.%M.%S")
        except ValueError:
            return datetime.max
    elif date_match:
        try:
            return datetime.strptime(date_match.group(1), "%Y-%m-%d")
        except ValueError:
            return datetime.max
    return datetime.max

def rename_and_collect_data(base_path):
    print("\n--- Step 1: Scanning & Renaming Source Data ---")
    
    # 1. FLAT DIRECTORY CHECK (Files directly in base_path)
    files_in_root = [f for f in os.listdir(base_path) if f.endswith((".ims", ".ome.tif"))]
    if files_in_root:
        print("  Detected flat directory structure (files in root).")
        total_size_bytes = 0
        file_extension = None
        for f in files_in_root:
            match = re.search(r'_F(\d+)(\.ims|\.ome\.tif)', f, re.IGNORECASE)
            if match:
                new_f = f"F{match.group(1)}{match.group(2)}"
                if f != new_f:
                    os.rename(os.path.join(base_path, f), os.path.join(base_path, new_f))
                    f = new_f
            total_size_bytes += os.path.getsize(os.path.join(base_path, f))
            if file_extension is None: 
                file_extension = ".ims" if f.endswith(".ims") else ".ome.tif"
                
        total_gb = total_size_bytes / (1024**3)
        print(f"  Identified 1 cycle.")
        print(f"  Total Data Size: {total_gb:.2f} GB")
        return 1, total_gb, ["."], file_extension

    # 2. SUBDIRECTORY CHECK (Multi-cycle or nested single cycle)
    all_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
    
    # Only process directories that have _Field in the name or actually contain image files
    valid_dirs = []
    for d in all_dirs:
        dir_path = os.path.join(base_path, d)
        if '_Field' in d or any(f.endswith((".ims", ".ome.tif")) for f in os.listdir(dir_path)):
            valid_dirs.append(d)
    
    dirs_with_date = [d for d in valid_dirs if re.search(r'\d{4}-\d{2}-\d{2}', d)]
    sorted_dirs = sorted(dirs_with_date, key=get_datetime_key)
    
    # Generate mapping for multi-cycle renaming
    dirs_to_rename_map = {}
    for i, d in enumerate(sorted_dirs):
        parts = d.split("_")
        base_name = parts[0]
        field_part = next((p for p in parts if "Field" in p), f"Field_{i+1}")
        new_name = f"{base_name}_{i + 1}_{field_part}"
        dirs_to_rename_map[d] = new_name

    # Rename directories
    for old, new in dirs_to_rename_map.items():
        if old != new:
            os.rename(os.path.join(base_path, old), os.path.join(base_path, new))
            print(f"  Renamed Dir: {old} -> {new}")

    def safe_sort_key(dir_name):
        match = re.search(r'_(\d+)_Field', dir_name)
        return int(match.group(1)) if match else -1
    
    # Collect newly named prospective directories
    current_all_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
    prospective_dirs_unsorted = [d for d in current_all_dirs if re.search(r'_\d+_Field', d) or any(f.endswith((".ims", ".ome.tif")) for f in os.listdir(os.path.join(base_path, d)))]
    prospective_dirs = sorted(prospective_dirs_unsorted, key=safe_sort_key)
    
    # Rename files inside directories
    total_size_bytes = 0
    file_extension = None
    valid_prospective_dirs = []
    
    for dir_name in prospective_dirs:
        dir_path = os.path.join(base_path, dir_name)
        has_files = False
        for f in os.listdir(dir_path):
            if f.endswith((".ims", ".ome.tif")):
                has_files = True
                match = re.search(r'_F(\d+)(\.ims|\.ome\.tif)', f, re.IGNORECASE)
                if match:
                    new_f = f"F{match.group(1)}{match.group(2)}"
                    if f != new_f:
                        os.rename(os.path.join(dir_path, f), os.path.join(dir_path, new_f))
                        f = new_f
                
                total_size_bytes += os.path.getsize(os.path.join(dir_path, f))
                if file_extension is None: 
                    file_extension = ".ims" if f.endswith(".ims") else ".ome.tif"
                    
        if has_files:
            valid_prospective_dirs.append(dir_name)

    total_gb = total_size_bytes / (1024**3)
    print(f"  Identified {len(valid_prospective_dirs)} cycles.")
    print(f"  Total Data Size: {total_gb:.2f} GB")
    
    return len(valid_prospective_dirs), total_gb, valid_prospective_dirs, file_extension

def parse_xml(base_path):
    print("\n--- Step 2: Parsing Metadata ---")
    
    # Use os.walk to search the root and all subdirectories for the XML file
    xml_path = next((os.path.join(root, f) for root, _, files in os.walk(base_path) for f in files if f.endswith('.xml')), None)
    
    if not xml_path:
        print("Error: No XML metadata file found anywhere in the data directory or its subdirectories.")
        sys.exit(1)
        
    tree = ET.parse(xml_path)
    root = tree.getroot()
    width = root.find('dimensions').attrib['stack_columns']
    height = root.find('dimensions').attrib['stack_rows']
    pixel_size = root.find('voxel_dims').attrib['H']
    
    print(f"  Found XML: {os.path.basename(xml_path)}")
    print(f"  Width: {width}, Height: {height}, Pixel Size: {pixel_size}")
    
    return width, height, pixel_size

def calculate_time(total_gb):
    ref_gb, ref_ashlar_hours = 226.0, 10.25
    ashlar_hours_per_gb = ref_ashlar_hours / ref_gb
    est_ashlar_hours = (total_gb * ashlar_hours_per_gb) * 1.44 + 2
    
    days, rem_hours = divmod(int(est_ashlar_hours), 24)
    if days > 0:
        return f"{days}-{rem_hours:02d}:00:00"
    else:
        hrs = max(2, int(est_ashlar_hours))
        return f"{hrs:02d}:00:00"

def generate_slurm_script(base_path, job_name, est_time, width, height, pixel_size, prospective_dirs, file_ext):
    slurm_file = os.path.join(base_path, "submit_ashlar.slurm")
    dir_strings = " ".join(prospective_dirs)
    
    script_content = f"""#!/bin/bash -l

# ==============================================================================
# USER CONFIGURATION - EDIT THESE VALUES BEFORE PROCEEDING
# ==============================================================================

# --- SLURM RESOURCES ---
#SBATCH --job-name=Ashlar_{job_name}
#SBATCH --account=hdid-share
#SBATCH --partition=tier1q
#SBATCH --time={est_time}
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=Ashlar_{job_name}.%J.out
#SBATCH --error=Ashlar_{job_name}.%J.err

# --- ASHLAR PARAMETERS ---
ALIGN_CHANNEL={user_align_channel}  # Passed securely via command line argument
OVERLAP=0.1
LAYOUT=snake
DIRECTION=vertical
MAX_SHIFT=30
FILTER_SIGMA=2

# ==============================================================================
# DO NOT EDIT BELOW (Unless adjusting core module paths)
# ==============================================================================

# Hardware specifics extracted from XML
WIDTH={width}
HEIGHT={height}
PIXEL_SIZE={pixel_size}

export PYTHONNOUSERSITE=1
cd {base_path}

# Load required modules
module load go/1.20.1 miniconda3 openjdk/17.0.2

# Activate environment
source activate {CONDA_ENV}

# Build the Ashlar command array
c=()
for dir in {dir_strings}; do
    c+=( "fileseries|${{dir}}|pattern=F{{series}}{file_ext}|width=${{WIDTH}}|height=${{HEIGHT}}|pixel_size=${{PIXEL_SIZE}}|overlap=${{OVERLAP}}|layout=${{LAYOUT}}|direction=${{DIRECTION}}" )
done

# Ensure output directory exists
mkdir -p registration

echo "Starting Ashlar processing..."
echo "Aligning on Channel: $ALIGN_CHANNEL"

ashlar "${{c[@]}}" \\
    --flip-y \\
    --pyramid \\
    --maximum-shift $MAX_SHIFT \\
    --filter-sigma $FILTER_SIGMA \\
    --align-channel $ALIGN_CHANNEL \\
    -o "registration/stitched.ome.tiff"

echo "Ashlar processing complete."
"""
    
    with open(slurm_file, "w") as f:
        f.write(script_content)
        
    return slurm_file

def main():
    base_path = source_data_dir  # Fixed: Now accurately uses the command line argument!
    job_name = os.path.basename(base_path.rstrip('/'))
    
    # 1. Rename & Collect
    num_cycles, total_gb, prospective_dirs, file_ext = rename_and_collect_data(base_path)
    if num_cycles == 0:
        print("Error: No valid cycle directories found.")
        sys.exit(1)
        
    # 2. Parse XML
    width, height, pixel_size = parse_xml(base_path)
    
    # 3. Calc Resources
    est_time = calculate_time(total_gb)
    print(f"\n  Calculated SLURM Time Limit: {est_time}")
    
    # 4. Generate Script
    slurm_file = generate_slurm_script(
        base_path, job_name, est_time, width, height, 
        pixel_size, prospective_dirs, file_ext
    )
    
    # 5. Pause & Prompt
    print("\n" + "="*70)
    print(" READY FOR REVIEW")
    print("="*70)
    print(f"  I have generated the submission script: {os.path.basename(slurm_file)}")
    print(f"  It is located in: {base_path}")
    print("  ")
    print("  Please open another terminal or use an editor to review/edit ")
    print("  the SLURM resources or Ashlar parameters in the USER CONFIGURATION")
    print("  block at the top of the file.")
    print("  ")
    print("  Press [ENTER] when you are ready to submit to the cluster.")
    print("  Press [Ctrl+C] to abort.")
    print("="*70)
    
    try:
        input()
    except KeyboardInterrupt:
        print("\nAborted. The SLURM script remains in the directory if you wish to submit it manually later.")
        sys.exit(0)
        
    # 6. Submit
    print("Submitting job...")
    os.chdir(base_path)
    try:
        result = subprocess.run(["sbatch", "submit_ashlar.slurm"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
        print(f"Success! Job Submitted: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"\nERROR SUBMITTING JOB\nStderr: {e.stderr.strip()}")

if __name__ == "__main__":
    main()