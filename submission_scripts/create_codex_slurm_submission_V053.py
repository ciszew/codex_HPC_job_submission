#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
MCMICRO Pipeline Submission Script for Slurm HPC
================================================================================
Version: 53 (Overlay Restoration & Archive Error Handling)

Description:
  This script automates the MCMICRO spatial multiplexed image analysis pipeline 
  on a Slurm-based High-Performance Computing (HPC) cluster. It structures the 
  workflow into a chain of 7 sequentially dependent Slurm jobs, staging data 
  to fast local scratch storage, registering/stitching tile images, running 
  preprocessing, segmenting cells and nuclei, generating QC report html/overlays, 
  and archiving final outputs back to persistent storage.

  A core highlight of the script is simulated index mapping. Users specify nuclei 
  and membrane markers by their exact string names. The script simulates the 
  omission of skipped (Ashlar skip) and removed (Backsub remove) channels from 
  markers_ashlar.csv to map exact 0-based indices for downstream steps, ensuring 
  perfect tool configuration.

The 7 Jobs in the Chain:
  1. 1_staging: 
     Creates a workspace in scratch storage and copies the raw cycle directories 
     and metadata (XML, CSV) from the source directory using rsync.
  2. 2_renaming: 
     Standardizes cycle directories and tile filenames (renaming cycle folders 
     to chronological order and files to F1.ims, F2.ims, etc.). Runs a size-based 
     corruption check, replacing files < 85% median size with healthy adjacent tiles.
  3. 3_ashlar: 
     Runs Ashlar tile registration and stitching to generate stitched.ome.tiff. 
     Uses an inline monkey-patched reader to skip/remove channels dynamically.
  4. 4_preprocess: 
     Runs Nextflow-based MCMICRO preprocessing, performing background subtraction. 
     For TMA runs, it also dearrays cores using Coreograph.
  5. 5_segment: 
     Runs Cell and Nuclear segmentation sequentially via Mesmer (Singularity containers). 
     Extracts outline/label masks and calls run_overlay_qc.sh to generate 
     outline overlays in the deepcell environment. Merges nuclei and cell masks into 
     standard cell.tif/nuclear.tif masks.
  6. 6_quant_merge: 
     Calculates expression metrics for each segmented cell (mcquant). Splits whole-cell 
     and nuclear quant outputs, evaluates segmentation quality via mcmicro_qc_analyzer.py 
     generating HTML reports, and aggregates results (batch merge for TMA, copying for WSI).
  7. 7_archive: 
     Copies registration, background subtraction, segmentation, quantification, QC 
     reports, and diagnostic logs back to the source folder using rsync, bypassing 
     missing folders cleanly.

Data Structure:
  - Input Source Directory:
      - Raw cycle directories containing fields/tiles (e.g. contain "_Field_" in folder name)
      - An XML metadata file (e.g., scan settings containing columns/rows dimensions)
      - markers_ashlar.csv: CSV mapping channels, cycle index, markers, and ashlar/remove flags
  - Scratch Run Directory:
      - standardized cycle folders (e.g. <Core_Run>_<Cycle>_Field)
      - registration/: Contains stitched.ome.tiff
      - background/: Background subtracted stitched TIFF
      - dearray/: (TMA only) Cropped core TIFFs
      - segmentation/: Final cell.tif and nuclear.tif mask files
      - quantification/: Cell and nuclear expression tables (CSV)
      - quantification_merged/: Consolidated quantification CSV tables
      - qc_overlays_wc/ / qc_overlays_nuc/: Multi-channel outline overlays
      - qc_report_wc.html / qc_report_nuc.html: Interactive HTML reports

Execution Instructions:
  Run the script directly on the login node passing the raw data directory as the target:
  
    python create_codex_slurm_submission_V053.py /path/to/raw_data_dir
    
  To generate the slurm job scripts and parameters without submitting them to the scheduler:
  
    python create_codex_slurm_submission_V053.py /path/to/raw_data_dir --dry-run

Troubleshooting:
  1. Job Failures: 
     Identify the failing job in the chain. The script prints submission logs indicating 
     job IDs. Check standard Slurm logs named `<Job_Key>.<Job_ID>.out` and `<Job_Key>.<Job_ID>.err` 
     placed directly in the source data directory.
  2. Diagnostic Logs: 
     Review detail-rich step logs written to the scratch run directory:
       `diag_<Job_Key>.<Job_ID>.log`
     These logs contain system hostname, scratch disk usage, path checks, and payload stdout.
  3. Resource Monitoring: 
     Check `<Job_Key>_resource_log.<Job_ID>.txt` in scratch to evaluate memory (RSS) and CPU 
     utilization curves. Useful to check if a job failed due to OOM or CPU throttling.
  4. Stale Mounts / Missing Tools: 
     Check `diag_1_staging.*.log` or job error files. Ensure nodes in EXCLUDE_NODES list 
     exclude problematic compute nodes.

Changes from V052:
  - BUG FIX: Restored dynamic '--is-tma' flag to run_overlay_qc.sh in Job 5,
             fixing the regression where overlays were skipped in WSI mode.
  - IMPROVEMENT: Removed silent error suppression from Job 7 archiving rsync commands,
             and added checks to skip missing directories/files (e.g. dearray on WSI)
             to prevent false failures while ensuring actual rsync failures propagate
             to Slurm via PAYLOAD_EXIT.
"""

import os
import re
import csv
from datetime import datetime
import xml.etree.ElementTree as ET
import subprocess
import sys
import argparse

# ==============================================================================
# --- USER-CONFIGURABLE VARIABLES ---
# ==============================================================================
# Login node mounts shared scratch at /scratch/; compute nodes mount the SAME
# filesystem at /ess/scratch/. Use the LOGIN NODE path here so Python makedirs
# works. #SBATCH --output/--error use $HOME (universally writable, see below).
SCRATCH_BASE_DIR = "/scratch/cciszews/nextflow_runs/06042026/"
CONDA_BASE_PATH = "/gpfs/data/hdid-share/conda"
QC_SCRIPTS_DIR = "/gpfs/data/hdid-share/Codex/HDID/scripts/current_working_scripts/"
EXCLUDE_NODES = "cri22cn140"  # Comma-separated list of nodes to exclude (e.g. stale mount nodes)

IS_TMA_WORKFLOW = False
IMAGE_MPP = 0.16
ARCSINH_COFACTOR = 5

# Set Segmentation Channels by EXACT MARKER NAME (as written in markers_ashlar.csv)
NUC_MARKER_NAME = "UV_high"
MEMBRANE_MARKER_NAMES = "CD45_Atto550 PanCK_AF488"

ASHLAR_OVERLAP = 0.1
ASHLAR_LAYOUT = 'snake'
ASHLAR_DIRECTION = 'vertical'
ASHLAR_MAX_SHIFT = 30
ASHLAR_FILTER_SIGMA = 2

# ==============================================================================
# --- DO NOT EDIT BELOW THIS LINE ---
# ==============================================================================

QC_ENV_DEEPCELL = os.path.join(CONDA_BASE_PATH, 'deepcell')
QC_ENV_REPORTS = os.path.join(CONDA_BASE_PATH, 'qc')

QC_SCRIPT_PATHS = {
    'overlay_sh': os.path.join(QC_SCRIPTS_DIR, 'run_overlay_qc.sh'),
    'overlay_py': os.path.join(QC_SCRIPTS_DIR, 'create_overlay_final.py'),
    'analyzer_py': os.path.join(QC_SCRIPTS_DIR, 'mcmicro_qc_analyzer.py'),
    'merger_py': os.path.join(QC_SCRIPTS_DIR, 'batch_merge_all_cores_v10.py')
}

def get_final_indices(markers_path, nuc_name, memb_names_str):
    """Calculates the final 0-based image indices by simulating dropped channels."""
    memb_names = memb_names_str.split()
    surviving_markers = []

    with open(markers_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ashlar_val = row.get('ashlar', '').strip().lower()
            remove_val = row.get('remove', '').strip().lower()

            # Channel survives if it isn't skipped by Ashlar AND isn't removed by Backsub
            if ashlar_val != 'skip' and remove_val not in ['true', 't', '1', 'yes']:
                surviving_markers.append(row.get('marker_name', '').strip())

    # Find the exact indices in the surviving list
    nuc_idx = surviving_markers.index(nuc_name) if nuc_name in surviving_markers else -1

    memb_indices = []
    for name in memb_names:
        if name in surviving_markers:
            memb_indices.append(surviving_markers.index(name))

    return nuc_idx, memb_indices

def get_data_info(path):
    all_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and '_Field_' in d]
    total_size_bytes, file_extension, dirs_to_rename_map = 0, None, {}

    def get_datetime_key(dir_name):
        date_match = re.search(r"(\d{4}-\d{2}-\d{2})", dir_name)
        time_match = re.search(r"(\d{2}\.\d{2}\.\d{2})", dir_name)
        if date_match and time_match:
            try: return datetime.strptime(f"{date_match.group(1)} {time_match.group(1)}", "%Y-%m-%d %H.%M.%S")
            except ValueError: pass
        return datetime.max

    dirs_with_date = [d for d in all_dirs if re.search(r'\d{4}-\d{2}-\d{2}', d)]
    sorted_dirs = sorted(dirs_with_date, key=get_datetime_key)

    for i, d in enumerate(sorted_dirs):
        name_match = re.search(r'^(.*?)_?(?=\d{4}-\d{2}-\d{2})', d)
        if name_match and name_match.group(1):
            raw_base = name_match.group(1).strip()
            base_name = re.sub(r'[\s\-]+', '_', raw_base)
        else: base_name = f"Core_Run"
        field_match = re.search(r'(Field_?\d*)', d, re.IGNORECASE)
        field_part = field_match.group(1) if field_match else "Field"
        dirs_to_rename_map[d] = f"{base_name}_{i + 1}_{field_part}"

    def safe_sort_key(dir_name):
        match = re.search(r'_(\d+)_Field', dir_name)
        return int(match.group(1)) if match else -1

    unsorted_dirs = list(set(all_dirs) - set(dirs_with_date))
    prospective_dirs = sorted(list(dirs_to_rename_map.values()) + unsorted_dirs, key=safe_sort_key)

    for d_name in all_dirs:
        dir_path = os.path.join(path, d_name)
        for f in os.listdir(dir_path):
            if f.endswith((".ims", ".ome.tif")):
                total_size_bytes += os.path.getsize(os.path.join(dir_path, f))
                if file_extension is None: file_extension = ".ims" if f.endswith(".ims") else ".ome.tif"

    return len(all_dirs), total_size_bytes / (1024**3), prospective_dirs, file_extension, dirs_to_rename_map

def format_slurm_time(hours):
    import math
    ceil_hours = math.ceil(hours)
    days, rem_hours = divmod(ceil_hours, 24)
    return f"{days}-{rem_hours:02d}:00:00" if days > 0 else f"{rem_hours:02d}:00:00"

def estimate_resources(total_gb, num_cores, is_tma):
    # Sub-linear time scaling
    est_ashlar_hours = 1.5 + 0.3 * (total_gb ** 0.6)
    
    if is_tma:
        # TMAs are cropped into small pieces, memory footprint is low
        est_segmentation_hours = 1.5 + 0.4 * (num_cores ** 0.6)
        est_quant_hours = 1.5 + 0.3 * (num_cores ** 0.6)
        seg_mem = "48G"
        ashlar_mem = "32G"
        quant_mem = "24G"
        prep_mem = "24G"
    else:
        # 40x Whole Slide Imaging loads massive arrays into RAM, scaled sub-linearly
        est_segmentation_hours = 2.0 + 0.35 * (total_gb ** 0.6)
        est_quant_hours = 1.5 + 0.15 * (total_gb ** 0.6)
        
        # Sub-linear memory scaling (GB)
        seg_mem = f"{int(24 + 4 * (total_gb ** 0.5))}G"
        ashlar_mem = f"{int(24 + 3 * (total_gb ** 0.5))}G"
        quant_mem = f"{int(16 + 2 * (total_gb ** 0.5))}G"
        prep_mem = f"{int(16 + 2 * (total_gb ** 0.5))}G"
        
    staging_hours = 1.0 + 0.1 * (total_gb ** 0.5)
    archive_hours = 1.0 + 0.1 * (total_gb ** 0.5)

    return {
        "1_staging":     {"time": format_slurm_time(staging_hours), "mem": "4G",   "cpu": "1"},
        "2_renaming":    {"time": "01:00:00", "mem": "4G",   "cpu": "1"},
        "3_ashlar":      {"time": format_slurm_time(est_ashlar_hours),      "mem": ashlar_mem,  "cpu": "4"},
        "4_preprocess":  {"time": "03:00:00", "mem": prep_mem,  "cpu": "4"},
        "5_segment":     {"time": format_slurm_time(est_segmentation_hours), "mem": seg_mem, "cpu": "4"},
        "6_quant_merge": {"time": format_slurm_time(est_quant_hours),        "mem": quant_mem,  "cpu": "4"},
        "7_archive":     {"time": format_slurm_time(archive_hours), "mem": "4G",   "cpu": "1"}
    }

def get_monitoring_script_block(log_prefix):
    return f"""
LOG_FILE="{log_prefix}_resource_log.${{SLURM_JOB_ID}}.txt"
echo "TimeElapsed|AveCPU|AveRSS|MaxRSS|AveDiskRead|MaxDiskRead|AveDiskWrite|MaxDiskWrite" > "${{LOG_FILE}}"
(while true; do
    ELAPSED=$(squeue -h -j $SLURM_JOB_ID -o %M);
    STATS=$(sstat --format=AveCPU,AveRSS,MaxRSS,AveDiskRead,MaxDiskRead,AveDiskWrite,MaxDiskWrite -P -n -j "${{SLURM_JOB_ID}}.batch" | tail -n1);
    echo "${{ELAPSED}}|${{STATS}}" >> $LOG_FILE;
    sleep 900;
done) &
MONITOR_PID=$!
trap "echo '>>> Cleaning up monitor process PID $MONITOR_PID...'; kill $MONITOR_PID" EXIT
"""

def get_diag_start_block(job_key, scratch_run_dir, check_paths=None):
    """
    Returns a bash prologue that initializes diagnostic logging for a job.
    Defines: DIAG_LOG path, diag_log() function, and PAYLOAD_EXIT tracker.
    check_paths: optional list of filesystem paths to verify at job start.
    """
    path_checks = ""
    if check_paths:
        for p in check_paths:
            path_checks += (
                f'\ndiag_log "  Checking path: {p}"'
                f'\nls -lhd "{p}" >> "${{DIAG_LOG}}" 2>&1'
                f' || diag_log "  WARNING: path not found or inaccessible: {p}"'
            )

    return f"""
# =============================================================================
# DIAGNOSTIC LOGGING PROLOGUE  [V049]
# Log file: {scratch_run_dir}/diag_{job_key}.<JOBID>.log
# =============================================================================
DIAG_LOG="{scratch_run_dir}/diag_{job_key}.${{SLURM_JOB_ID}}.log"
diag_log() {{ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${{DIAG_LOG}}"; }}

# PAYLOAD_EXIT tracks the highest-priority failure exit code across the job.
# It is set after each critical command and used as the final job exit code.
PAYLOAD_EXIT=0

diag_log "========================================"
diag_log "=== JOB START: {job_key} ==="
diag_log "========================================"
diag_log "SLURM_JOB_ID       : ${{SLURM_JOB_ID}}"
diag_log "SLURM_JOB_NAME     : ${{SLURM_JOB_NAME}}"
diag_log "SLURM_NODELIST     : ${{SLURM_NODELIST}}"
diag_log "SLURM_CPUS_ON_NODE : ${{SLURM_CPUS_ON_NODE}}"
diag_log "SLURM_MEM_PER_NODE : ${{SLURM_MEM_PER_NODE}} MB"
diag_log "Hostname           : $(hostname)"
diag_log "Start time         : $(date)"
diag_log "--- Scratch disk space ---"
df -h "{scratch_run_dir}" >> "${{DIAG_LOG}}" 2>&1
diag_log "--- Key path availability checks ---"{path_checks}
diag_log "--- Startup diagnostics complete. Beginning job payload. ---"
# =============================================================================
"""

def get_diag_end_block(job_key):
    """
    Returns a bash epilogue that finalizes the diagnostic log and exits
    with PAYLOAD_EXIT (set throughout the job payload after critical commands).
    Must be placed as the very last block in the job script.
    """
    return f"""
# =============================================================================
# DIAGNOSTIC LOGGING EPILOGUE  [V049]
# =============================================================================
diag_log "End time           : $(date)"
diag_log "========================================"
diag_log "=== JOB END: {job_key} === PAYLOAD_EXIT: $PAYLOAD_EXIT"
diag_log "========================================"
exit $PAYLOAD_EXIT
"""

def generate_rename_script(path, dirs_to_rename_map, xml_path):
    return f'''
import os, re, statistics, shutil

path = "{path}"
dirs_map = {dirs_to_rename_map}

for old, new in dirs_map.items():
    try: os.rename(os.path.join(path, old), os.path.join(path, new))
    except OSError: pass

all_cycles = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and re.search(r'_\\d+_Field', d)]
all_cycles.sort(key=lambda x: int(re.search(r'_(\\d+)_Field', x).group(1)))

log_file = open(os.path.join(path, "QA_QC_Report.txt"), "w")
log_file.write("--- Data Integrity Check ---\\n")

for dir_name in all_cycles:
    dir_path = os.path.join(path, dir_name)

    for f in os.listdir(dir_path):
        match = re.search(r'_F(\\d+)(\\.ims|\\.ome\\.tif)', f, re.IGNORECASE)
        if match:
            new_name = f"F{{match.group(1)}}{{match.group(2)}}"
            try: os.rename(os.path.join(dir_path, f), os.path.join(dir_path, new_name))
            except OSError: pass

    files_in_dir = [f for f in os.listdir(dir_path) if f.startswith('F') and f.endswith('.ims')]
    if not files_in_dir: continue

    file_sizes = [os.path.getsize(os.path.join(dir_path, f)) for f in files_in_dir]
    median_size = statistics.median(file_sizes)
    threshold = median_size * 0.85

    files_in_dir.sort(key=lambda x: int(re.search(r'F(\\d+)', x).group(1)))

    for i, f_name in enumerate(files_in_dir):
        file_path = os.path.join(dir_path, f_name)
        size = os.path.getsize(file_path)

        if size < threshold:
            log_file.write(f"WARNING: Corrupted file detected: {{dir_name}}/{{f_name}} ({{size}} bytes. Median: {{median_size}})\\n")
            replacement = None
            if i > 0 and os.path.getsize(os.path.join(dir_path, files_in_dir[i-1])) > threshold: replacement = files_in_dir[i-1]
            elif i < len(files_in_dir) - 1 and os.path.getsize(os.path.join(dir_path, files_in_dir[i+1])) > threshold: replacement = files_in_dir[i+1]

            if replacement:
                shutil.copy2(os.path.join(dir_path, replacement), file_path)
                log_file.write(f"  -> Successfully patched by duplicating adjacent tile: {{replacement}}\\n")
            else: log_file.write(f"  -> ERROR: Could not find healthy adjacent tile in this cycle.\\n")

log_file.close()
'''

def create_ashlar_script(dirs, ext, xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    width = root.find('dimensions').attrib['stack_columns']
    height = root.find('dimensions').attrib['stack_rows']
    pixel_size = root.find('voxel_dims').attrib['H']
    align_channel = -1

    return f"""
WIDTH={width}
HEIGHT={height}
PIXEL_SIZE={pixel_size}
OVERLAP={ASHLAR_OVERLAP}
LAYOUT={ASHLAR_LAYOUT}
DIRECTION={ASHLAR_DIRECTION}
MAX_SHIFT={ASHLAR_MAX_SHIFT}
SIGMA={ASHLAR_FILTER_SIGMA}
ALIGN_CHANNEL={align_channel}

cat << 'EOF' > ashlar_wrapper.py
import sys, csv, ashlar.fileseries
from ashlar.scripts.ashlar import main as ashlar_main

cycle_to_kept = {{}}
clean_rows = []

with open('markers_ashlar.csv', 'r') as f:
    reader = csv.reader(f)
    headers = next(reader)
    cycle_idx, ashlar_idx, channel_idx = headers.index('cycle'), headers.index('ashlar'), headers.index('channel')

    clean_headers = [h for h in headers if h.lower() != 'ashlar']
    clean_rows.append(clean_headers)

    current_cycle, relative_idx, new_channel_counter = -1, 0, 1

    for row in reader:
        if not row: continue
        c_cycle = int(row[cycle_idx])
        ashlar_val = row[ashlar_idx].strip().lower() if ashlar_idx < len(row) else ''

        if c_cycle != current_cycle:
            current_cycle = c_cycle
            relative_idx = 0
            cycle_to_kept[current_cycle] = []

        if ashlar_val != 'skip':
            cycle_to_kept[current_cycle].append(relative_idx)
            clean_row = []
            for i, val in enumerate(row):
                if i == ashlar_idx: continue
                if i == channel_idx:
                    clean_row.append(str(new_channel_counter))
                    new_channel_counter += 1
                else: clean_row.append(val)
            clean_rows.append(clean_row)
        relative_idx += 1

with open('markers.csv', 'w', newline='') as f:
    csv.writer(f).writerows(clean_rows)

original_init = ashlar.fileseries.FileSeriesReader.__init__
original_read = ashlar.fileseries.FileSeriesReader.read
reader_counter = 0

def patched_init(self, *args, **kwargs):
    global reader_counter
    original_init(self, *args, **kwargs)
    cycle_idx = reader_counter
    reader_counter += 1

    kept_raw_indices = cycle_to_kept.get(cycle_idx, list(range(self.metadata.num_channels)))

    self.my_raw_map = {{virtual_idx: raw_idx for virtual_idx, raw_idx in enumerate(kept_raw_indices)}}
    if kept_raw_indices: self.my_raw_map[-1] = kept_raw_indices[-1]

    class MockMetadata(self.metadata.__class__):
        @property
        def num_channels(self): return self._patched_num_channels
    self.metadata.__class__ = MockMetadata
    self.metadata._patched_num_channels = len(kept_raw_indices)

def patched_read(self, series, c):
    return original_read(self, series, self.my_raw_map.get(c, c))

ashlar.fileseries.FileSeriesReader.__init__ = patched_init
ashlar.fileseries.FileSeriesReader.read = patched_read

if __name__ == '__main__': sys.exit(ashlar_main())
EOF

c=()
for dir in {' '.join(dirs)}; do
    c+=( "fileseries|${{dir}}|pattern=F{{series}}{ext}|width=${{WIDTH}}|height=${{HEIGHT}}|pixel_size=${{PIXEL_SIZE}}|overlap=${{OVERLAP}}|layout=${{LAYOUT}}|direction=${{DIRECTION}}" )
done
mkdir -p registration
python3 ./ashlar_wrapper.py "${{c[@]}}" --flip-y --pyramid --maximum-shift $MAX_SHIFT --filter-sigma $SIGMA --align-channel=$ALIGN_CHANNEL -o "registration/stitched.ome.tiff"
"""

def get_params_preprocess_yml(is_tma):
    tma_str = 'true' if is_tma else 'false'
    stop_at = 'dearray' if is_tma else 'background'
    return f"""
workflow:
  background: true
  tma: {tma_str}
  start-at: background
  stop-at: {stop_at}
options:
  coreograph: --channel 0 --downsampleFactor 8
modules:
  dearray:
    name: coreograph
    version: 2.4.6
"""

def get_params_analysis_yml(is_tma, compartment, nuc_idx, memb_indices):
    tma_str = 'true' if is_tma else 'false'

    # Recyze uses 1-based extraction logic
    nuc_1b = nuc_idx + 1
    all_channels = [str(nuc_1b)] + [str(i + 1) for i in memb_indices]

    # Format strictly as a space-separated string enclosed in quotes: '1 13 15'
    yaml_list = "'" + " ".join(all_channels) + "'"

    return f"""
workflow:
  background: true
  tma: {tma_str}
  start-at: segmentation
  stop-at: segmentation
  segmentation: [mesmer]
  segmentation-channel: {yaml_list}
  segmentation-recyze: true
options:
  mesmer: --image-mpp {IMAGE_MPP} --compartment {compartment}
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('data_dir', type=str, help='Path to the raw data directory.')
    parser.add_argument('--dry-run', action='store_true', help='Generate scripts without submitting.')
    args = parser.parse_args()

    source_data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(source_data_dir):
        print(f"Error: Directory not found '{source_data_dir}'"); sys.exit(1)

    markers_ashlar_path = os.path.join(source_data_dir, 'markers_ashlar.csv')
    if not os.path.exists(markers_ashlar_path):
        print("Error: markers_ashlar.csv not found in data directory."); sys.exit(1)

    # Resolve Dynamic Channel Names
    print(f"\n--- Resolving Marker Names ---")
    try:
        nuc_idx, memb_indices = get_final_indices(markers_ashlar_path, NUC_MARKER_NAME, MEMBRANE_MARKER_NAMES)
        if nuc_idx == -1:
            print(f"ERROR: Nuclear marker '{NUC_MARKER_NAME}' not found or is marked to skip/remove in markers_ashlar.csv."); sys.exit(1)
        if len(memb_indices) != len(MEMBRANE_MARKER_NAMES.split()):
            print(f"ERROR: One or more membrane markers from '{MEMBRANE_MARKER_NAMES}' not found or marked to skip/remove."); sys.exit(1)

        print(f"  [Nuclear Marker]  {NUC_MARKER_NAME} -> Mapped to post-backsub Index {nuc_idx}")
        print(f"  [Membrane Marker] {MEMBRANE_MARKER_NAMES} -> Mapped to post-backsub Indices {memb_indices}")
    except Exception as e:
        print(f"Error reading markers_ashlar.csv: {e}"); sys.exit(1)

    xml_path = next((os.path.join(root, f) for root, _, files in os.walk(source_data_dir) for f in files if f.endswith('.xml')), None)
    if not xml_path: print("Error: No XML metadata file found."); sys.exit(1)

    print("\n--- Scanning Source Data ---")
    num_cores, total_gb, prospective_dirs, file_extension, dirs_to_rename_map = get_data_info(source_data_dir)
    resources = estimate_resources(total_gb, num_cores, IS_TMA_WORKFLOW)

    run_name = os.path.basename(source_data_dir)
    scratch_run_dir = os.path.join(SCRATCH_BASE_DIR, run_name)
    staged_scripts_dir = os.path.join(scratch_run_dir, "scripts")
    home_dir = os.environ.get('HOME', os.path.expanduser('~')).rstrip('/')

    print("\n" + "="*60)
    print("--- MCMICRO SUBMISSION SUMMARY (V53) ---")
    print("="*60)
    print(f"Data Directory:         {source_data_dir}")
    print(f"Run Name:               {run_name}")
    print(f"Total Cycles/Fields:    {num_cores}")
    print(f"Total Raw Data Size:    {total_gb:.2f} GB")
    print(f"Scratch Directory:      {scratch_run_dir}")
    print(f"Slurm .out/.err dir:    {source_data_dir}")
    print(f"Diagnostic logs:        {scratch_run_dir}/diag_<job>.<jobid>.log")
    print("\n--- Estimated Resources per Job ---")
    for job_name, res in resources.items():
        print(f"  - {job_name:<15} Time={res['time']}, Mem={res['mem']}, CPUs={res['cpu']}")
    print("="*60 + "\n")

    if not args.dry_run:
        os.makedirs(scratch_run_dir, exist_ok=True)
        os.makedirs(staged_scripts_dir, exist_ok=True)
        for path in QC_SCRIPT_PATHS.values(): subprocess.run(['cp', path, staged_scripts_dir], check=True)
        subprocess.run(['cp', markers_ashlar_path, scratch_run_dir], check=True)

    scripts_to_write = {
        os.path.join(staged_scripts_dir, '_renamer.py'): generate_rename_script(scratch_run_dir, dirs_to_rename_map, os.path.join(scratch_run_dir, os.path.basename(xml_path))),
        os.path.join(scratch_run_dir, 'ashlar.sh'): create_ashlar_script(prospective_dirs, file_extension, xml_path),
        os.path.join(scratch_run_dir, 'params-preprocess.yml'): get_params_preprocess_yml(IS_TMA_WORKFLOW),
        os.path.join(scratch_run_dir, 'params-wc.yml'): get_params_analysis_yml(IS_TMA_WORKFLOW, 'whole-cell', nuc_idx, memb_indices),
        os.path.join(scratch_run_dir, 'params-nuc.yml'): get_params_analysis_yml(IS_TMA_WORKFLOW, 'nuclear', nuc_idx, memb_indices),
        os.path.join(scratch_run_dir, 'params-quant.yml'): f"workflow:\n  background: true\n  tma: {'true' if IS_TMA_WORKFLOW else 'false'}\n  start-at: quantification\n  stop-at: quantification\noptions:\n  mcquant: --masks cell.tif nuclear.tif --intensity_props intensity_median"
    }

    if not args.dry_run:
        for path, content in scripts_to_write.items():
            with open(path, 'w') as f: f.write(content)

    # =========================================================================
    # generate_slurm_header uses the resolved absolute path of the user's home directory.
    # The scratch filesystem has different mount points on login vs compute nodes
    # (/scratch/ on login, /ess/scratch/ on compute). The home directory (/home/... or
    # /ess/home/...) is confirmed reachable and writable from ALL compute nodes.
    # Python resolves $HOME on the login node and writes absolute paths to prevent
    # Slurm from failing to expand the literal '$HOME' string.
    # =========================================================================
    def generate_slurm_header(job_key, res):
        # Automatically switch to tier2q if memory is > 128G
        mem_gb = int(res['mem'].replace('G', ''))
        partition = "tier2q" if mem_gb > 128 else "tier1q"

        out_path = os.path.join(source_data_dir, f"{job_key}.%j.out")
        err_path = os.path.join(source_data_dir, f"{job_key}.%j.err")
        exclude_line = f"#SBATCH --exclude={EXCLUDE_NODES}\n" if EXCLUDE_NODES else ""
        return (
            f"#!/bin/bash -l\n"
            f"#SBATCH --job-name={job_key}\n"
            f"#SBATCH --account=hdid-share\n"
            f"#SBATCH --partition={partition}\n"
            f"#SBATCH --time={res['time']}\n"
            f"#SBATCH --cpus-per-task={res['cpu']}\n"
            f"#SBATCH --mem={res['mem']}\n"
            f"#SBATCH --output={out_path}\n"
            f"#SBATCH --error={err_path}\n"
            {exclude_line}\n"
        )

    # Grab the first membrane index for DeepCell QC overlays
    cyto_idx = memb_indices[0] if len(memb_indices) > 0 else nuc_idx
    tma_flag = " --is-tma" if IS_TMA_WORKFLOW else ""

    # =========================================================================
    # Job payloads — each critical command captures its exit code into PAYLOAD_EXIT
    # =========================================================================

    # --- Job 1: Data Staging ---
    j1_payload = f"""
mkdir -p {scratch_run_dir}
sleep 5
diag_log "Starting rsync: {source_data_dir}/ -> {scratch_run_dir}/"
rsync -ah --info=progress2 {source_data_dir}/ {scratch_run_dir}/
RSYNC_EXIT=$?
diag_log "EXIT [rsync staging]: $RSYNC_EXIT"
[ $RSYNC_EXIT -ne 0 ] && PAYLOAD_EXIT=$RSYNC_EXIT && diag_log "*** FAILED: rsync staging ***"
diag_log "--- Post-rsync scratch disk usage ---"
du -sh "{scratch_run_dir}" >> "${{DIAG_LOG}}" 2>&1
diag_log "--- Scratch root listing ---"
ls -lh "{scratch_run_dir}" >> "${{DIAG_LOG}}" 2>&1
"""

    # --- Job 2: Renaming + File QC ---
    j2_payload = f"""
module load go/1.20.1 miniconda3
diag_log "Running renamer script: {staged_scripts_dir}/_renamer.py"
python {staged_scripts_dir}/_renamer.py
RENAME_EXIT=$?
diag_log "EXIT [_renamer.py]: $RENAME_EXIT"
[ $RENAME_EXIT -ne 0 ] && PAYLOAD_EXIT=$RENAME_EXIT && diag_log "*** FAILED: _renamer.py ***"
diag_log "--- QA_QC_Report.txt contents ---"
cat "{scratch_run_dir}/QA_QC_Report.txt" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: QA_QC_Report.txt not found"
diag_log "--- Renamed cycle dirs in scratch ---"
ls -lh "{scratch_run_dir}" >> "${{DIAG_LOG}}" 2>&1
"""

    # --- Job 3: Ashlar Stitching ---
    j3_payload = f"""
export PYTHONNOUSERSITE=1
cd {scratch_run_dir}
module load go/1.20.1 miniconda3 openjdk/17.0.2
source activate {CONDA_BASE_PATH}/ashlar_group
diag_log "Starting Ashlar stitching (bash ./ashlar.sh)"
bash ./ashlar.sh
ASHLAR_EXIT=$?
diag_log "EXIT [ashlar.sh]: $ASHLAR_EXIT"
[ $ASHLAR_EXIT -ne 0 ] && PAYLOAD_EXIT=$ASHLAR_EXIT && diag_log "*** FAILED: ashlar.sh ***"
diag_log "--- Registration output listing ---"
ls -lh "{scratch_run_dir}/registration/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: registration/ not found after Ashlar"
"""

    # --- Job 4: Background Subtraction (MCMICRO/Backsub via Nextflow) ---
    j4_payload = f"""
cd {scratch_run_dir}
module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow
diag_log "--- Pre-nextflow: registration/ listing ---"
ls -lh "{scratch_run_dir}/registration/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: registration/ not found before nextflow"
diag_log "Starting Nextflow preprocess (background subtraction) with params-preprocess.yml"
nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-preprocess.yml
NF_PRE_EXIT=$?
diag_log "EXIT [nextflow preprocess]: $NF_PRE_EXIT"
[ $NF_PRE_EXIT -ne 0 ] && PAYLOAD_EXIT=$NF_PRE_EXIT && diag_log "*** FAILED: nextflow preprocess ***"
diag_log "--- Post-nextflow: background/ listing ---"
ls -lh "{scratch_run_dir}/background/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: background/ not found after nextflow"
diag_log "--- Nextflow work dir disk usage ---"
du -sh "{scratch_run_dir}/work/" >> "${{DIAG_LOG}}" 2>&1 || true
"""

    # --- Job 5: Segmentation (WC + Nuclear) + QC Overlays ---
    j5_payload = f"""
cd {scratch_run_dir}
module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3

diag_log "=== PHASE 1: Whole-cell segmentation (params-wc.yml) ==="
nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-wc.yml
NF_WC_EXIT=$?
diag_log "EXIT [nextflow WC segmentation]: $NF_WC_EXIT"
[ $NF_WC_EXIT -ne 0 ] && PAYLOAD_EXIT=$NF_WC_EXIT && diag_log "*** FAILED: nextflow WC segmentation ***"

diag_log "Activating DeepCell env: {QC_ENV_DEEPCELL}"
source activate {QC_ENV_DEEPCELL}
bash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel {cyto_idx} --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer{tma_flag}
QC_WC_EXIT=$?
diag_log "EXIT [run_overlay_qc.sh WC]: $QC_WC_EXIT"
[ $QC_WC_EXIT -ne 0 ] && PAYLOAD_EXIT=$QC_WC_EXIT && diag_log "*** FAILED: WC QC overlay ***"
conda deactivate

diag_log "Moving WC outputs: segmentation -> segmentation_wc"
mv segmentation segmentation_wc
mv qc_overlays_segmentation-mesmer qc_overlays_wc
diag_log "--- segmentation_wc/ listing ---"
ls -lh "{scratch_run_dir}/segmentation_wc/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: segmentation_wc/ not found"

diag_log "=== PHASE 2: Nuclear segmentation (params-nuc.yml) ==="
nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-nuc.yml
NF_NUC_EXIT=$?
diag_log "EXIT [nextflow nuclear segmentation]: $NF_NUC_EXIT"
[ $NF_NUC_EXIT -ne 0 ] && PAYLOAD_EXIT=$NF_NUC_EXIT && diag_log "*** FAILED: nextflow nuclear segmentation ***"

diag_log "Activating DeepCell env: {QC_ENV_DEEPCELL}"
source activate {QC_ENV_DEEPCELL}
bash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel {cyto_idx} --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer{tma_flag}
QC_NUC_EXIT=$?
diag_log "EXIT [run_overlay_qc.sh nuclear]: $QC_NUC_EXIT"
[ $QC_NUC_EXIT -ne 0 ] && PAYLOAD_EXIT=$QC_NUC_EXIT && diag_log "*** FAILED: nuclear QC overlay ***"
conda deactivate

diag_log "Renaming nuclear masks: cell.tif -> nuclear.tif"
for dir in segmentation/mesmer-*; do mv "$dir/cell.tif" "$dir/nuclear.tif" 2>/dev/null; done
mv segmentation segmentation_nuc
mv qc_overlays_segmentation-mesmer qc_overlays_nuc

diag_log "=== PHASE 3: Merging WC + nuclear into final segmentation/ ==="
mkdir -p segmentation
for dir in segmentation_nuc/mesmer-*; do core_id=$(basename $dir); mkdir -p segmentation/$core_id; cp segmentation_wc/$core_id/cell.tif segmentation/$core_id/ 2>/dev/null; cp $dir/nuclear.tif segmentation/$core_id/ 2>/dev/null; done
diag_log "--- Final merged segmentation/ listing ---"
ls -lhR "{scratch_run_dir}/segmentation/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: segmentation/ not found after merge"
"""

    # Define TMA core merging vs WSI direct copying
    if IS_TMA_WORKFLOW:
        merge_block = f"""
diag_log "Running batch merge: batch_merge_all_cores_v10.py"
source activate {QC_ENV_REPORTS}
python ./scripts/batch_merge_all_cores_v10.py --project-dir .
MERGE_EXIT=$?
diag_log "EXIT [batch_merge_all_cores_v10.py]: $MERGE_EXIT"
[ $MERGE_EXIT -ne 0 ] && PAYLOAD_EXIT=$MERGE_EXIT && diag_log "*** FAILED: batch_merge_all_cores_v10.py ***"
conda deactivate
"""
    else:
        merge_block = f"""
diag_log "WSI Workflow: Skipping batch merge and copying single whole-slide CSVs to quantification_merged/"
mkdir -p quantification_merged
cp quantification/*cell.csv quantification_merged/cell.csv 2>/dev/null
cp quantification/*nuclear.csv quantification_merged/nuclear.csv 2>/dev/null
MERGE_EXIT=$?
diag_log "EXIT [WSI mock merge]: $MERGE_EXIT"
"""

    # --- Job 6: Quantification + QC Reports + Batch Merge ---
    j6_payload = f"""
cd {scratch_run_dir}
module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3

diag_log "Starting Nextflow quantification (params-quant.yml)"
nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-quant.yml
NF_QUANT_EXIT=$?
diag_log "EXIT [nextflow quantification]: $NF_QUANT_EXIT"
[ $NF_QUANT_EXIT -ne 0 ] && PAYLOAD_EXIT=$NF_QUANT_EXIT && diag_log "*** FAILED: nextflow quantification ***"

diag_log "Splitting WC/nuclear quantification CSVs"
mkdir -p quantification_wc quantification_nuc
mv quantification/*cell.csv quantification_wc/
mv quantification/*nuclear.csv quantification_nuc/
diag_log "--- quantification_wc/ listing ---"
ls -lh "{scratch_run_dir}/quantification_wc/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: quantification_wc/ not found"

diag_log "Activating QC env: {QC_ENV_REPORTS}"
source activate {QC_ENV_REPORTS}

diag_log "Running QC analyzer (whole-cell)"
python ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_wc --output-path qc_report_wc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}
QC_ANA_WC_EXIT=$?
diag_log "EXIT [mcmicro_qc_analyzer.py WC]: $QC_ANA_WC_EXIT"
[ $QC_ANA_WC_EXIT -ne 0 ] && PAYLOAD_EXIT=$QC_ANA_WC_EXIT && diag_log "*** FAILED: QC analyzer WC ***"

diag_log "Running QC analyzer (nuclear)"
python ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_nuc --output-path qc_report_nuc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}
QC_ANA_NUC_EXIT=$?
diag_log "EXIT [mcmicro_qc_analyzer.py nuclear]: $QC_ANA_NUC_EXIT"
[ $QC_ANA_NUC_EXIT -ne 0 ] && PAYLOAD_EXIT=$QC_ANA_NUC_EXIT && diag_log "*** FAILED: QC analyzer nuclear ***"
conda deactivate

diag_log "Consolidating quantification CSVs back into quantification/"
mv quantification_wc/* quantification/
mv quantification_nuc/* quantification/
rmdir quantification_wc quantification_nuc

{merge_block}

diag_log "--- Final quantification/ listing ---"
ls -lh "{scratch_run_dir}/quantification/" >> "${{DIAG_LOG}}" 2>&1 || diag_log "  WARNING: quantification/ not found"
"""

    # --- Job 7: Archive Results to Source ---
    j7_payload = f"""
diag_log "Starting archive rsync from scratch to source data dir"
for item in registration background dearray segmentation quantification_merged qc_overlays_wc qc_overlays_nuc qc_report_wc.html qc_report_nuc.html QA_QC_Report.txt; do
    if [ -e "{scratch_run_dir}/$item" ]; then
        diag_log "  Archiving item: $item"
        rsync -ah --info=progress2 "{scratch_run_dir}/$item" "{source_data_dir}/"
        RSYNC_ITEM_EXIT=$?
        diag_log "  EXIT [rsync $item]: $RSYNC_ITEM_EXIT"
        [ $RSYNC_ITEM_EXIT -ne 0 ] && PAYLOAD_EXIT=$RSYNC_ITEM_EXIT && diag_log "*** FAILED: rsync archiving of $item ***"
    else
        diag_log "  Skipping archive for missing item: $item"
    fi
done
diag_log "Archive rsync loop complete"
"""

    # =========================================================================
    # Assemble final job scripts
    # =========================================================================
    job_scripts = {
        '1_staging': (
            generate_slurm_header('1_staging', resources['1_staging']) +
            get_diag_start_block('1_staging', scratch_run_dir,
                check_paths=[source_data_dir]) +
            j1_payload +
            get_diag_end_block('1_staging')
        ),
        '2_renaming': (
            generate_slurm_header('2_renaming', resources['2_renaming']) +
            get_diag_start_block('2_renaming', scratch_run_dir,
                check_paths=[f"{staged_scripts_dir}/_renamer.py", scratch_run_dir]) +
            j2_payload +
            get_diag_end_block('2_renaming')
        ),
        '3_ashlar': (
            generate_slurm_header('3_ashlar', resources['3_ashlar']) +
            get_diag_start_block('3_ashlar', scratch_run_dir,
                check_paths=[f"{scratch_run_dir}/ashlar.sh",
                             f"{scratch_run_dir}/markers_ashlar.csv"]) +
            get_monitoring_script_block(os.path.join(scratch_run_dir, '3_ashlar')) +
            j3_payload +
            get_diag_end_block('3_ashlar')
        ),
        '4_preprocess': (
            generate_slurm_header('4_preprocess', resources['4_preprocess']) +
            get_diag_start_block('4_preprocess', scratch_run_dir,
                check_paths=[f"{scratch_run_dir}/registration/stitched.ome.tiff",
                             f"{scratch_run_dir}/params-preprocess.yml"]) +
            get_monitoring_script_block(os.path.join(scratch_run_dir, '4_preprocess')) +
            j4_payload +
            get_diag_end_block('4_preprocess')
        ),
        '5_segment': (
            generate_slurm_header('5_segment', resources['5_segment']) +
            get_diag_start_block('5_segment', scratch_run_dir,
                check_paths=[f"{scratch_run_dir}/background",
                             f"{scratch_run_dir}/params-wc.yml",
                             f"{scratch_run_dir}/params-nuc.yml"]) +
            get_monitoring_script_block(os.path.join(scratch_run_dir, '5_segment')) +
            j5_payload +
            get_diag_end_block('5_segment')
        ),
        '6_quant_merge': (
            generate_slurm_header('6_quant_merge', resources['6_quant_merge']) +
            get_diag_start_block('6_quant_merge', scratch_run_dir,
                check_paths=[f"{scratch_run_dir}/segmentation",
                             f"{scratch_run_dir}/params-quant.yml"]) +
            get_monitoring_script_block(os.path.join(scratch_run_dir, '6_quant_merge')) +
            j6_payload +
            get_diag_end_block('6_quant_merge')
        ),
        '7_archive': (
            generate_slurm_header('7_archive', resources['7_archive']) +
            get_diag_start_block('7_archive', scratch_run_dir,
                check_paths=[f"{scratch_run_dir}/registration",
                             f"{scratch_run_dir}/quantification_merged"]) +
            j7_payload +
            get_diag_end_block('7_archive')
        ),
    }

    slurm_paths = {}
    if not args.dry_run:
        for k, v in job_scripts.items():
            p = os.path.join(scratch_run_dir, f"{k}.slurm")
            with open(p, 'w') as f: f.write(v)
            slurm_paths[k] = p

    if args.dry_run: sys.exit(0)

    # =========================================================================
    # V049 FIX: submit_job() validates the sbatch-returned Job ID.
    # An empty or non-numeric ID causes an immediate RuntimeError, preventing
    # a corrupted dependency chain (--dependency=afterok:) from being submitted.
    # =========================================================================
    def submit_job(cmd):
        """Submits a Slurm job via sbatch and validates the returned Job ID."""
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        jid = result.stdout.strip()
        if not jid.isdigit():
            stderr_msg = result.stderr.strip()
            raise RuntimeError(
                f"sbatch failed.\n"
                f"  Command: {' '.join(cmd)}\n"
                f"  stdout : '{jid}'\n"
                f"  stderr : '{stderr_msg}'"
            )
        script_name = os.path.basename(cmd[-1])
        print(f"  Submitted Job {jid}  ({script_name})")
        return jid

    try:
        print("\n--- Submitting 7-Job Chain ---")
        j1 = submit_job(["sbatch", "--parsable", slurm_paths['1_staging']])
        j2 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j1}", slurm_paths['2_renaming']])
        j3 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j2}", slurm_paths['3_ashlar']])
        j4 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j3}", slurm_paths['4_preprocess']])
        j5 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j4}", slurm_paths['5_segment']])
        j6 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j5}", slurm_paths['6_quant_merge']])
        j7 = submit_job(["sbatch", "--parsable", f"--dependency=afterok:{j6}", slurm_paths['7_archive']])
        print(f"\nSuccess! 7-Job Chain Submitted.")
        print(f"  Chain     : {j1} -> {j2} -> {j3} -> {j4} -> {j5} -> {j6} -> {j7}")
        print(f"  Diag logs : {scratch_run_dir}/diag_<job>.<jobid>.log")
        print(f"  Slurm out : {source_data_dir}/<job>.<jobid>.out/.err")
    except RuntimeError as e:
        print(f"\nFATAL: Job chain submission aborted.\n{e}")
        sys.exit(1)
