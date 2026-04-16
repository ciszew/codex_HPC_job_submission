import os
import re
import sys
import h5py
import numpy as np
import tifffile
from datetime import datetime

# =============================================================================
# --- CONFIGURATION ---
# =============================================================================

# 1. The QA/QC Report file
QA_QC_REPORT_PATH = "/gpfs/data/hdid-share/Codex/OlopadeLab/04022026/QA_QC_Report.txt"

# 2. The RAW data directory (Contains the 13 chronologically named folders)
RAW_DATA_DIR = "/gpfs/data/hdid-share/Codex/OlopadeLab/04022026/"

# 3. The CONVERTED data directory (Where the patched files currently live for Ashlar)
ASHLAR_DATA_DIR = "/ess/scratch/scratch1/angreene/codex_runs/04022026/"

# Crop & Channel settings
EXPECTED_CHANNELS = 5
TRUE_HEIGHT = 1817
TRUE_WIDTH = 1979

# =============================================================================

def parse_report_for_corrupted_files(report_path):
    """Scans the QA/QC report and extracts the relative paths of corrupted files."""
    corrupted_paths = []
    with open(report_path, 'r') as f:
        for line in f:
            if "WARNING: Corrupted file detected:" in line:
                # Extract the path string: "Sample2_31ab_4_Field_1/F1373.ims"
                match = re.search(r'detected:\s*(\S+\.ims)', line)
                if match:
                    corrupted_paths.append(match.group(1))
    return corrupted_paths

def map_to_raw_file(relative_path, raw_dir):
    """Translates 'Sample2_31ab_4_Field_1/F1373.ims' to the actual raw file path."""
    
    match = re.search(r'_(\d+)_Field.*?/(F\d+)\.ims', relative_path)
    if not match:
        raise ValueError(f"Could not parse cycle number or tile ID from: {relative_path}")
        
    cycle_num = int(match.group(1))
    tile_id = match.group(2)
    
    raw_subdirs = [d for d in os.listdir(raw_dir) if os.path.isdir(os.path.join(raw_dir, d))]
    raw_subdirs.sort()
    
    if cycle_num > len(raw_subdirs) or cycle_num < 1:
        raise IndexError(f"Report asked for Cycle {cycle_num}, but only {len(raw_subdirs)} raw folders exist!")
        
    # Python 0-indexing: Cycle 1 = Index 0, Cycle 4 = Index 3
    target_raw_folder = raw_subdirs[cycle_num - 1]
    full_target_raw_folder = os.path.join(raw_dir, target_raw_folder)
    
    for file in os.listdir(full_target_raw_folder):
        if file.endswith(f"_{tile_id}.ims"):
            return os.path.join(full_target_raw_folder, file)
            
    raise FileNotFoundError(f"Could not find raw file ending in _{tile_id}.ims in {target_raw_folder}")

def rescue_and_replace(raw_file_path, converted_target_path, audit_log):
    """Extracts data from raw IMS, pads channels, and overwrites target as a fake .ims"""
    
    base_name = os.path.basename(converted_target_path)
    missing_channels = []
    dead_channels = []
    
    try:
        # Open source file in strict read-only mode
        with h5py.File(raw_file_path, 'r') as f:
            time_point = f['DataSet']['ResolutionLevel 0']['TimePoint 0']
            found_channels = list(time_point.keys())
            
            image_stack = []
            
            for i in range(EXPECTED_CHANNELS):
                expected_ch = f"Channel {i}"
                
                if expected_ch in found_channels:
                    data = time_point[expected_ch]['Data'][:]
                    cropped_data = np.squeeze(data)[0:TRUE_HEIGHT, 0:TRUE_WIDTH]
                    
                    if np.max(cropped_data) == 0:
                        dead_channels.append(i)
                        image_stack.append(cropped_data) 
                    else:
                        image_stack.append(cropped_data) 
                else:
                    missing_channels.append(i)
                    blank_layer = np.zeros((TRUE_HEIGHT, TRUE_WIDTH), dtype=np.uint16)
                    image_stack.append(blank_layer)
            
            final_image = np.stack(image_stack)
            
            # --- THE TROJAN HORSE OVERWRITE ---
            # We use tifffile to save the array, but we overwrite the exact .ims file path.
            # This completely replaces the patched file with our new TIFF-structured data.
            tifffile.imwrite(converted_target_path, final_image, imagej=True)
            
            # --- LOGGING ---
            status_msg = f"{datetime.now().strftime('%H:%M:%S')} | OVERWROTE: {base_name} "
            if missing_channels:
                status_msg += f"| Missing (Injected Blanks): {missing_channels} "
            if dead_channels:
                status_msg += f"| Dead (All Zeros): {dead_channels}"
            if not missing_channels and not dead_channels:
                status_msg += "| 100% Intact"
                
            print(status_msg)
            audit_log.write(status_msg + "\n")
            
    except Exception as e:
        err_msg = f"{datetime.now().strftime('%H:%M:%S')} | FAILED to rescue {base_name}: {e}"
        print(err_msg)
        audit_log.write(err_msg + "\n")

def main():
    print("--- Starting Retroactive Ashlar Data Rescue (Uniform Extension Mode) ---")
    
    corrupted_files = parse_report_for_corrupted_files(QA_QC_REPORT_PATH)
    print(f"Found {len(corrupted_files)} corrupted files in the QA/QC report.\n")
    
    audit_log_path = os.path.join(ASHLAR_DATA_DIR, f"Rescue_Audit_Log_{datetime.now().strftime('%Y%m%d_%H%M')}.txt")
    
    with open(audit_log_path, 'w') as audit_log:
        audit_log.write("--- RETROACTIVE DATA RESCUE AUDIT LOG ---\n")
        audit_log.write(f"Executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for rel_path in corrupted_files:
            try:
                raw_path = map_to_raw_file(rel_path, RAW_DATA_DIR)
                target_converted_path = os.path.join(ASHLAR_DATA_DIR, rel_path)
                
                # Execute the overwrite
                rescue_and_replace(raw_path, target_converted_path, audit_log)
                
            except Exception as e:
                err_msg = f"ERROR mapping or processing {rel_path}: {e}"
                print(err_msg)
                audit_log.write(err_msg + "\n")

    print(f"\n--- Rescue Complete. Audit log saved to: {audit_log_path} ---")

if __name__ == "__main__":
    main()