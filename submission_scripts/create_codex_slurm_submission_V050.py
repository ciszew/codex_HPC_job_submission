#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
MCMICRO Pipeline Submission Script for Slurm HPC
================================================================================
Version: 50 (Inline Rescue Architecture)

Purpose:
  Automates the MCMICRO pipeline. Users specify target channel NAMES instead
  of indices. The script calculates the exact final indices by simulating
  the Ashlar drop (ashlar=skip) and Backsub drop (remove=TRUE).

Changes from V49:
  - Replaced legacy shutil.copy2() corruption handler in the generated
    _renamer.py with an inline h5py/tifffile zero-padding rescue.
  - When a corrupted .ims file (file size < 85% of median) is detected,
    the renamer now:
      1. Opens the file via h5py (read-only).
      2. Extracts surviving channels from DataSet/ResolutionLevel 0/TimePoint 0.
      3. Zero-pads missing channels with np.zeros(dtype=uint16).
      4. Overwrites the file via tifffile.imwrite() keeping the .ims extension
         ("Trojan Horse" technique for downstream compatibility).
  - Job 2 (2_renaming) now activates the ashlar_group Conda environment
    to provide h5py, numpy, tifffile at runtime.
  - Job 2 resource allocation increased to 16G RAM to support HDF5 reads.
  - New user-configurable variables: EXPECTED_CHANNELS, TRUE_HEIGHT, TRUE_WIDTH.
  - All V49 security hardening preserved (shlex.quote, json.dumps, submit_chain).
"""

import os
import re
import csv
import sys
import json
import shlex
import argparse
import subprocess
from datetime import datetime
import xml.etree.ElementTree as ET

# ==============================================================================
# --- USER-CONFIGURABLE VARIABLES ---
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. PIPELINE MODE & MICROSCOPE SETTINGS
# ------------------------------------------------------------------------------
# IS_TMA_WORKFLOW: Set to True for Tissue Microarrays, False for Whole Slide Images (WSI).
# This dynamically adjusts Nextflow routing, DeepCell QC flags, and Slurm resource math.
IS_TMA_WORKFLOW = False

# IMAGE_MPP: Microns per pixel. Critical for DeepCell/Mesmer segmentation scale.
# (e.g., 0.16 for 40x magnification, 0.32 for standard 20x).
IMAGE_MPP = 0.16

# ARCSINH_COFACTOR: Used for arcsinh data transformation in the final QC reports. 
# Standard Codex default is 5.
ARCSINH_COFACTOR = 5

# ------------------------------------------------------------------------------
# 2. MARKER CONFIGURATION
# ------------------------------------------------------------------------------
# IMPORTANT: These names must EXACTLY match the 'marker_name' column in your 
# markers_ashlar.csv file (case-insensitive). The script will automatically calculate 
# the correct indices by dropping 'ashlar=skip' and 'remove=true' channels.

# NUC_MARKER_NAME: The exact string name of the nuclear marker (e.g., "DAPI", "UV_high").
NUC_MARKER_NAME = "UV_high"

# MEMBRANE_MARKER_NAMES: Space-separated string of membrane/cytoplasm markers 
# used for whole-cell segmentation. 
# Example: "CD45_Atto550 PanCK_AF750 CD3e_AF488"
MEMBRANE_MARKER_NAMES = "CD45_Atto550 PanCK_AF488"

# ------------------------------------------------------------------------------
# 3. ASHLAR STITCHING PARAMETERS
# ------------------------------------------------------------------------------
ASHLAR_OVERLAP = 0.1          # Expected overlap between microscope tiles (10% = 0.1)
ASHLAR_LAYOUT = 'snake'       # Tile acquisition layout ('snake' or 'raster')
ASHLAR_DIRECTION = 'vertical' # Acquisition direction ('vertical' or 'horizontal')
ASHLAR_MAX_SHIFT = 30         # Maximum allowed shift in pixels during alignment
ASHLAR_FILTER_SIGMA = 2       # Smoothing filter for noisy channels (reduces stitch errors)

# ------------------------------------------------------------------------------
# 4. HPC RESOURCE LIMITS
# ------------------------------------------------------------------------------
# MAX_WALL_HOURS: The absolute maximum time limit for your Slurm partition (168h = 7 days).
# If dataset estimates exceed this, the script will clamp to this value and print a warning.
MAX_WALL_HOURS = 168

# ------------------------------------------------------------------------------
# 5. INLINE RESCUE — CORRUPTED TILE RECOVERY
# ------------------------------------------------------------------------------
# These parameters control the zero-padding rescue for corrupted .ims tiles.
# When the renamer detects a file smaller than 85% of the cycle's median size,
# it reads surviving channels via h5py and pads missing channels with zeros.
#
# EXPECTED_CHANNELS: Total number of channels expected per tile (e.g., 5 for CODEX).
# TRUE_HEIGHT: Pixel height to crop each channel slice to. Must match acquisition ROI.
# TRUE_WIDTH:  Pixel width to crop each channel slice to. Must match acquisition ROI.
EXPECTED_CHANNELS = 5
TRUE_HEIGHT = 1817
TRUE_WIDTH = 1979

# ==============================================================================
# --- CONFIGURATION RESOLUTION ---
# ==============================================================================

# --- Script Defaults (Priority 4 — lowest) ---
_DEFAULTS = {
    "scratch_base_dir": "/scratch/cciszews/nextflow_runs/04142026/",
    "conda_base_path": "/gpfs/data/hdid-share/conda",
    "qc_scripts_dir": "/gpfs/data/hdid-share/Codex/HDID/scripts/current_working_scripts/",
}

# --- Env-var-to-config key mapping (Priority 2) ---
_ENV_MAP = {
    "CODEX_SCRATCH_DIR": "scratch_base_dir",
    "CODEX_CONDA_BASE": "conda_base_path",
    "CODEX_QC_SCRIPTS_DIR": "qc_scripts_dir",
}

def load_config(cli_args):
    """Resolve configuration with cascading priority: CLI > Env > File > Defaults."""
    # Priority 4: defaults
    config = dict(_DEFAULTS)

    # Priority 3: config file in data directory
    config_path = os.path.join(os.path.abspath(cli_args.data_dir), '.codex_pipeline.json')
    if os.path.isfile(config_path):
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                file_cfg = json.load(f)
            if not isinstance(file_cfg, dict):
                print(f"WARNING: {config_path} does not contain a JSON object. Ignoring.")
            else:
                config.update(file_cfg)
        except (json.JSONDecodeError, OSError) as e:
            print(f"WARNING: Could not read config file {config_path}: {e}")

    # Priority 2: environment variables
    for env_key, config_key in _ENV_MAP.items():
        val = os.environ.get(env_key)
        if val:
            config[config_key] = val

    # Priority 1: CLI overrides
    if getattr(cli_args, 'scratch_dir', None):
        config["scratch_base_dir"] = cli_args.scratch_dir
    if getattr(cli_args, 'conda_base', None):
        config["conda_base_path"] = cli_args.conda_base
    if getattr(cli_args, 'qc_scripts_dir', None):
        config["qc_scripts_dir"] = cli_args.qc_scripts_dir

    return config


# ==============================================================================
# --- HELPER: CASE-INSENSITIVE COLUMN LOOKUP ---
# ==============================================================================

def _ci_lookup(fieldnames, target):
    """Return the actual column name matching *target* (case-insensitive). Raises SystemExit on miss."""
    target_lower = target.lower()
    for name in fieldnames:
        if name.lower() == target_lower:
            return name
    raise SystemExit(
        f"FATAL: Required CSV column '{target}' not found. "
        f"Available columns: {fieldnames}"
    )


def _strip_bom(fieldnames):
    """Strip UTF-8 BOM from the first fieldname and whitespace from all."""
    if not fieldnames:
        return fieldnames
    cleaned = list(fieldnames)
    cleaned[0] = cleaned[0].lstrip('\ufeff')
    return [h.strip() for h in cleaned]


# ==============================================================================
# --- CORE PIPELINE FUNCTIONS ---
# ==============================================================================

def get_final_indices(markers_path, nuc_name, memb_names_str):
    """Calculates the final 0-based image indices by simulating dropped channels.

    Raises SystemExit with a diagnostic dump of surviving markers if any
    requested marker name is not found.
    """
    memb_names = memb_names_str.split()
    surviving_markers = []

    with open(markers_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        raw_fields = _strip_bom(reader.fieldnames or [])

        # Build case-insensitive column mapping
        col_marker = _ci_lookup(raw_fields, 'marker_name')
        col_ashlar = _ci_lookup(raw_fields, 'ashlar')
        col_remove = _ci_lookup(raw_fields, 'remove')

        for row in reader:
            ashlar_val = row.get(col_ashlar, '').strip().lower()
            remove_val = row.get(col_remove, '').strip().lower()

            # Channel survives if it isn't skipped by Ashlar AND isn't removed by Backsub
            if ashlar_val != 'skip' and remove_val not in ('true', 't', '1', 'yes'):
                surviving_markers.append(row.get(col_marker, '').strip())

    # --- Resolve indices with diagnostic output on failure ---
    if nuc_name not in surviving_markers:
        raise SystemExit(
            f"FATAL: Nuclear marker '{nuc_name}' not found in surviving channels.\n"
            f"  Surviving markers ({len(surviving_markers)}): {surviving_markers}"
        )
    nuc_idx = surviving_markers.index(nuc_name)

    memb_indices = []
    for name in memb_names:
        if name not in surviving_markers:
            raise SystemExit(
                f"FATAL: Membrane marker '{name}' not found in surviving channels.\n"
                f"  Surviving markers ({len(surviving_markers)}): {surviving_markers}"
            )
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
    days, rem_hours = divmod(int(hours), 24)
    return f"{days}-{rem_hours:02d}:00:00" if days > 0 else f"{int(hours):02d}:00:00"


def _clamp_wall_time(estimated_hours, job_label):
    """Clamp an estimated wall time to MAX_WALL_HOURS and warn if exceeded."""
    if estimated_hours > MAX_WALL_HOURS:
        print(
            f"WARNING: Estimated {job_label} time ({estimated_hours:.0f}h) exceeds "
            f"{MAX_WALL_HOURS}h partition limit. Clamping to {MAX_WALL_HOURS}h. "
            f"Consider splitting into array jobs or reducing core count.",
            file=sys.stderr,
        )
        return MAX_WALL_HOURS
    return estimated_hours


def estimate_resources(total_gb, num_cores, is_tma):
    ref_gb, ref_ashlar_hours = 326.0, 10.25
    ashlar_hours_per_gb = ref_ashlar_hours / ref_gb
    est_ashlar_hours = (total_gb * ashlar_hours_per_gb) * 1.5 + 2
    est_ashlar_hours = _clamp_wall_time(est_ashlar_hours, "ashlar")

    # Dynamic scaling for Job 4 (Preprocess)
    est_prep_mem = max(48, int(32 + (total_gb * 0.02)))
    est_prep_hours = max(4, int(3 + (total_gb * 0.003)))
    est_prep_hours = _clamp_wall_time(est_prep_hours, "preprocess")

    if is_tma:
        # TMAs are cropped into small pieces, memory footprint is low
        est_segmentation_hours = (num_cores * 0.8) + 3
        est_quant_hours = (num_cores * 0.6) + 3
        seg_mem = "96G"
    else:
        # 40x Whole Slide Imaging loads massive arrays into RAM
        est_segmentation_hours = (total_gb * 0.05) + 6
        est_quant_hours = (total_gb * 0.02) + 4
        seg_mem = "160G"

    est_segmentation_hours = _clamp_wall_time(est_segmentation_hours, "segmentation")
    est_quant_hours = _clamp_wall_time(est_quant_hours, "quantification")

    return {
        "1_staging":     {"time": "04:00:00", "mem": "8G", "cpu": "1"},
        "2_renaming":    {"time": "01:00:00", "mem": "16G", "cpu": "1"},
        "3_ashlar":      {"time": format_slurm_time(est_ashlar_hours), "mem": "64G", "cpu": "4"},
        "4_preprocess":  {"time": format_slurm_time(est_prep_hours), "mem": f"{est_prep_mem}G", "cpu": "4"},
        "5_segment":     {"time": format_slurm_time(est_segmentation_hours), "mem": seg_mem, "cpu": "4"},
        "6_quant_merge": {"time": format_slurm_time(est_quant_hours), "mem": "64G", "cpu": "4"},
        "7_archive":     {"time": "04:00:00", "mem": "8G", "cpu": "1"}
    }


def get_monitoring_script_block(log_prefix):
    safe_prefix = shlex.quote(log_prefix)
    return f"""
LOG_FILE={safe_prefix}"_resource_log.${{SLURM_JOB_ID}}.txt"
echo "TimeElapsed|AveCPU|AveRSS|MaxRSS|AveDiskRead|MaxDiskRead|AveDiskWrite|MaxDiskWrite" > "${{LOG_FILE}}"
(while true; do
    ELAPSED=$(squeue -h -j $SLURM_JOB_ID -o %M);
    STATS=$(sstat --format=AveCPU,AveRSS,MaxRSS,AveDiskRead,MaxDiskRead,AveDiskWrite,MaxDiskWrite -P -n -j "${{SLURM_JOB_ID}}.batch" | tail -n1);
    echo "${{ELAPSED}}|${{STATS}}" >> "$LOG_FILE";
    sleep 900;
done) &
MONITOR_PID=$!
trap "echo '>>> Cleaning up monitor process PID $MONITOR_PID...'; kill $MONITOR_PID" EXIT
"""


# ==============================================================================
# --- INLINE RESCUE: Generated Python Source for _renamer.py ---
# ==============================================================================
# This string is NOT an f-string. It is plain source code that gets prepended
# to the generated renamer script. No brace-escaping needed.

_RESCUE_FUNC_SRC = r'''
def rescue_corrupted_tile(file_path, log_file):
    """Inline zero-padding rescue for corrupted .ims tiles.

    Reads surviving channels from the HDF5 structure inside the .ims file,
    pads any missing channels with np.zeros(dtype=uint16), and overwrites
    the file as a TIFF using tifffile.imwrite() while keeping the .ims
    extension for downstream compatibility ('Trojan Horse' technique).
    """
    base_name = os.path.basename(file_path)
    missing_channels = []
    dead_channels = []

    try:
        with h5py.File(file_path, 'r') as f:
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
                    missing_channels.append(i)
                    blank_layer = np.zeros((TRUE_HEIGHT, TRUE_WIDTH), dtype=np.uint16)
                    image_stack.append(blank_layer)

            final_image = np.stack(image_stack)

            # Overwrite the corrupted .ims file with valid TIFF data
            tifffile.imwrite(file_path, final_image, imagej=True)

        # --- Build detailed log entry ---
        ts = datetime.now().strftime('%H:%M:%S')
        status_msg = f"{ts} | RESCUED: {base_name} "
        if missing_channels:
            status_msg += f"| Missing (Injected Blanks): {missing_channels} "
        if dead_channels:
            status_msg += f"| Dead (All Zeros): {dead_channels}"
        if not missing_channels and not dead_channels:
            status_msg += "| All channels present (partial data recovered)"

        print(status_msg)
        log_file.write(status_msg + "\n")
        return True

    except Exception as e:
        ts = datetime.now().strftime('%H:%M:%S')
        err_msg = (
            f"{ts} | RESCUE FAILED: {base_name}: {e}\n"
            f"  -> HALTING: This tile could not be rescued. "
            f"Manual intervention required."
        )
        print(err_msg)
        log_file.write(err_msg + "\n")
        return False
'''


def generate_rename_script(path, dirs_to_rename_map, xml_path,
                           expected_channels, true_height, true_width):
    """Generate the _renamer.py source code with inline rescue capability.

    The rescue configuration values are injected via repr() at the top of the
    generated script, avoiding any f-string brace-escaping issues. The rescue
    function itself is prepended as a raw string constant (_RESCUE_FUNC_SRC).
    """
    # Build the main body as an f-string (same pattern as V049)
    main_body = f'''
import os, re, statistics, json
import h5py
import numpy as np
import tifffile
from datetime import datetime

path = {repr(path)}
dirs_map = json.loads({repr(json.dumps(dirs_to_rename_map))})

# --- Inline Rescue Configuration ---
EXPECTED_CHANNELS = {repr(expected_channels)}
TRUE_HEIGHT = {repr(true_height)}
TRUE_WIDTH = {repr(true_width)}

for old, new in dirs_map.items():
    try: os.rename(os.path.join(path, old), os.path.join(path, new))
    except OSError as e: raise RuntimeError(f"FATAL: Directory rename failed {{old}} -> {{new}}: {{e}}")

all_cycles = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and re.search(r'_\\d+_Field', d)]
all_cycles.sort(key=lambda x: int(re.search(r'_(\\d+)_Field', x).group(1)))

log_file = open(os.path.join(path, "QA_QC_Report.txt"), "w")
log_file.write("--- Data Integrity Check ---\\n")

rescue_failure_count = 0

for dir_name in all_cycles:
    dir_path = os.path.join(path, dir_name)
    
    for f in os.listdir(dir_path):
        match = re.search(r'_F(\\d+)(\\.ims|\\.ome\\.tif)', f, re.IGNORECASE)
        if match: 
            new_name = f"F{{match.group(1)}}{{match.group(2)}}"
            try: os.rename(os.path.join(dir_path, f), os.path.join(dir_path, new_name))
            except OSError as e: raise RuntimeError(f"FATAL: File rename failed {{f}} -> {{new_name}}: {{e}}")
            
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
            success = rescue_corrupted_tile(file_path, log_file)
            if not success:
                rescue_failure_count += 1

if rescue_failure_count > 0:
    log_file.write(f"\\n!!! CRITICAL: {{rescue_failure_count}} tile(s) failed rescue. Review log above. !!!\\n")
    print(f"CRITICAL: {{rescue_failure_count}} tile(s) could not be rescued. See QA_QC_Report.txt.")

log_file.close()
'''
    # Concatenate: rescue function definition FIRST, then the main body
    return _RESCUE_FUNC_SRC + main_body


def create_ashlar_script(dirs, ext, xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    width = root.find('dimensions').attrib['stack_columns']
    height = root.find('dimensions').attrib['stack_rows']
    pixel_size = root.find('voxel_dims').attrib['H']
    align_channel = -1

    # Shell-quote every directory name for safe bash interpolation
    quoted_dirs = ' '.join(shlex.quote(d) for d in dirs)

    return f"""
WIDTH={shlex.quote(width)}
HEIGHT={shlex.quote(height)}
PIXEL_SIZE={shlex.quote(pixel_size)}
OVERLAP={ASHLAR_OVERLAP}
LAYOUT={ASHLAR_LAYOUT}
DIRECTION={ASHLAR_DIRECTION}
MAX_SHIFT={ASHLAR_MAX_SHIFT}
SIGMA={ASHLAR_FILTER_SIGMA}
ALIGN_CHANNEL={align_channel}

cat << 'EOF' > ashlar_wrapper.py
import sys, csv, json, ashlar.fileseries
from ashlar.scripts.ashlar import main as ashlar_main

def get_col_idx(headers, name):
    name_lower = name.lower()
    for i, h in enumerate(headers):
        if h.lower() == name_lower:
            return i
    raise SystemExit(
        f"FATAL: Column '{{name}}' not found in markers_ashlar.csv. "
        f"Available columns: {{headers}}"
    )

cycle_to_kept = {{}}
clean_rows = []

with open('markers_ashlar.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    headers = [h.strip().lstrip('\\ufeff') for h in next(reader)]
    cycle_idx = get_col_idx(headers, 'cycle')
    ashlar_idx = get_col_idx(headers, 'ashlar')
    channel_idx = get_col_idx(headers, 'channel')

    clean_headers = [h for h in headers if h.lower() != 'ashlar']
    clean_rows.append(clean_headers)

    current_cycle, relative_idx, new_channel_counter = -1, 0, 1

    for row in reader:
        if not row or all(cell.strip() == '' for cell in row):
            continue
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
for dir in {quoted_dirs}; do
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
    
    if compartment == 'nuclear':
        all_channels = [str(nuc_1b)]
    else:
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


# ==============================================================================
# --- SLURM ORCHESTRATION ---
# ==============================================================================

def _rollback_jobs(job_ids):
    """Cancel all previously submitted Slurm jobs."""
    if job_ids:
        subprocess.run(["scancel"] + job_ids)
        print(f"  Rolled back jobs: {', '.join(job_ids)}")


def submit_chain(slurm_paths, job_order):
    """Submit a Slurm job chain with validation and automatic rollback on failure.

    For each job:
      1. Builds the sbatch command with --parsable.
      2. Appends --dependency=afterok:<prev_id> if not the first job.
      3. Validates that sbatch returns exit code 0.
      4. Validates that the returned job ID is numeric.
      5. On any failure: cancels all previously submitted jobs and exits.

    Returns the list of submitted job IDs.
    """
    submitted_ids = []

    for job_key in job_order:
        if job_key not in slurm_paths:
            print(f"FATAL: No slurm script path for job '{job_key}'.")
            _rollback_jobs(submitted_ids)
            sys.exit(1)

        cmd = ["sbatch", "--parsable"]
        if submitted_ids:
            cmd.append(f"--dependency=afterok:{submitted_ids[-1]}")
        cmd.append(slurm_paths[job_key])

        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        if result.returncode != 0:
            print(f"FATAL: sbatch failed for {job_key}: {result.stderr.strip()}")
            _rollback_jobs(submitted_ids)
            sys.exit(1)

        job_id = result.stdout.strip()
        if not re.match(r'^\d+$', job_id):
            print(f"FATAL: sbatch returned non-numeric job ID for {job_key}: '{job_id}'")
            _rollback_jobs(submitted_ids)
            sys.exit(1)

        submitted_ids.append(job_id)
        print(f"  Submitted {job_key}: Job ID {job_id}")

    return submitted_ids


# ==============================================================================
# --- MAIN ---
# ==============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('data_dir', type=str, help='Path to the raw data directory.')
    parser.add_argument('--dry-run', action='store_true', help='Generate scripts without submitting.')
    parser.add_argument('--scratch-dir', type=str, default=None,
                        help='Override scratch base directory (Priority 1).')
    parser.add_argument('--conda-base', type=str, default=None,
                        help='Override conda base path (Priority 1).')
    parser.add_argument('--qc-scripts-dir', type=str, default=None,
                        help='Override QC scripts directory (Priority 1).')
    args = parser.parse_args()

    # --- Resolve configuration ---
    cfg = load_config(args)
    SCRATCH_BASE_DIR = cfg["scratch_base_dir"]
    CONDA_BASE_PATH = cfg["conda_base_path"]
    QC_SCRIPTS_DIR = cfg["qc_scripts_dir"]

    QC_ENV_DEEPCELL = os.path.join(CONDA_BASE_PATH, 'deepcell')
    QC_ENV_REPORTS = os.path.join(CONDA_BASE_PATH, 'qc')

    QC_SCRIPT_PATHS = {
        'overlay_sh': os.path.join(QC_SCRIPTS_DIR, 'run_overlay_qc.sh'),
        'overlay_py': os.path.join(QC_SCRIPTS_DIR, 'create_overlay_final.py'),
        'analyzer_py': os.path.join(QC_SCRIPTS_DIR, 'mcmicro_qc_analyzer.py'),
        'merger_py': os.path.join(QC_SCRIPTS_DIR, 'batch_merge_all_cores_v10.py')
    }

    source_data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(source_data_dir):
        print(f"Error: Directory not found '{source_data_dir}'"); sys.exit(1)

    markers_ashlar_path = os.path.join(source_data_dir, 'markers_ashlar.csv')
    if not os.path.exists(markers_ashlar_path):
        print("Error: markers_ashlar.csv not found in data directory."); sys.exit(1)

    # Resolve Dynamic Channel Names
    print(f"\n--- Resolving Marker Names ---")
    nuc_idx, memb_indices = get_final_indices(markers_ashlar_path, NUC_MARKER_NAME, MEMBRANE_MARKER_NAMES)
    # get_final_indices now raises SystemExit on failure — no sentinel checks needed
    print(f"  [Nuclear Marker]  {NUC_MARKER_NAME} -> Mapped to post-backsub Index {nuc_idx}")
    print(f"  [Membrane Marker] {MEMBRANE_MARKER_NAMES} -> Mapped to post-backsub Indices {memb_indices}")

    xml_path = next((os.path.join(root, f) for root, _, files in os.walk(source_data_dir) for f in files if f.endswith('.xml')), None)
    if not xml_path: print("Error: No XML metadata file found."); sys.exit(1)

    print("\n--- Scanning Source Data ---")
    num_cores, total_gb, prospective_dirs, file_extension, dirs_to_rename_map = get_data_info(source_data_dir)
    resources = estimate_resources(total_gb, num_cores, IS_TMA_WORKFLOW)

    run_name = os.path.basename(source_data_dir)
    scratch_run_dir = os.path.join(SCRATCH_BASE_DIR, run_name)
    staged_scripts_dir = os.path.join(scratch_run_dir, "scripts")

    print("\n" + "="*60)
    print("--- MCMICRO SUBMISSION SUMMARY (V50) ---")
    print("="*60)
    print(f"Data Directory:         {source_data_dir}")
    print(f"Run Name:               {run_name}")
    print(f"Total Cycles/Fields:    {num_cores}")
    print(f"Total Raw Data Size:    {total_gb:.2f} GB")
    print(f"Scratch Directory:      {scratch_run_dir}")
    print(f"Conda Base Path:        {CONDA_BASE_PATH}")
    print(f"QC Scripts Directory:   {QC_SCRIPTS_DIR}")
    print(f"Inline Rescue Config:   {EXPECTED_CHANNELS}ch, {TRUE_HEIGHT}x{TRUE_WIDTH}px")
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
        os.path.join(staged_scripts_dir, '_renamer.py'): generate_rename_script(
            scratch_run_dir, dirs_to_rename_map,
            os.path.join(scratch_run_dir, os.path.basename(xml_path)),
            EXPECTED_CHANNELS, TRUE_HEIGHT, TRUE_WIDTH
        ),
        os.path.join(scratch_run_dir, 'ashlar.sh'): create_ashlar_script(prospective_dirs, file_extension, xml_path),
        os.path.join(scratch_run_dir, 'params-preprocess.yml'): get_params_preprocess_yml(IS_TMA_WORKFLOW),
        os.path.join(scratch_run_dir, 'params-wc.yml'): get_params_analysis_yml(IS_TMA_WORKFLOW, 'whole-cell', nuc_idx, memb_indices),
        os.path.join(scratch_run_dir, 'params-nuc.yml'): get_params_analysis_yml(IS_TMA_WORKFLOW, 'nuclear', nuc_idx, memb_indices),
        os.path.join(scratch_run_dir, 'params-quant.yml'): "workflow:\n  background: true\n  tma: true\n  start-at: quantification\n  stop-at: quantification\noptions:\n  mcquant: --masks cell.tif nuclear.tif --intensity_props intensity_median"
    }

    if not args.dry_run:
        for path, content in scripts_to_write.items():
            with open(path, 'w') as f: f.write(content)

    def generate_slurm_header(job_key, res):
        # Automatically switch to tier2q if memory is > 128G
        mem_gb = int(res['mem'].replace('G', ''))
        partition = "tier2q" if mem_gb > 128 else "tier1q"

        return f"#!/bin/bash -l\n#SBATCH --job-name={job_key}\n#SBATCH --account=hdid-share\n#SBATCH --partition={partition}\n#SBATCH --time={res['time']}\n#SBATCH --cpus-per-task={res['cpu']}\n#SBATCH --mem={res['mem']}\n#SBATCH --output={job_key}.%J.out\n#SBATCH --error={job_key}.%J.err\n"

    # Grab the first membrane index for DeepCell QC overlays
    cyto_idx_str = " ".join(map(str, memb_indices)) if len(memb_indices) > 0 else str(nuc_idx)

    # --- Build bash variable blocks for safe path interpolation ---
    def _path_vars_block(src, scratch, scripts, conda, qc_dc, qc_rpt):
        """Generate a bash block that assigns paths to shell variables (double-quoted)."""
        return (
            f'SOURCE_DIR={shlex.quote(src)}\n'
            f'SCRATCH_DIR={shlex.quote(scratch)}\n'
            f'SCRIPTS_DIR={shlex.quote(scripts)}\n'
            f'CONDA_BASE={shlex.quote(conda)}\n'
            f'QC_ENV_DC={shlex.quote(qc_dc)}\n'
            f'QC_ENV_RPT={shlex.quote(qc_rpt)}\n'
        )

    path_block = _path_vars_block(
        source_data_dir, scratch_run_dir, staged_scripts_dir,
        CONDA_BASE_PATH, QC_ENV_DEEPCELL, QC_ENV_REPORTS
    )

    tma_flag = " --is-tma" if IS_TMA_WORKFLOW else ""

    job_scripts = {
        '1_staging': (
            generate_slurm_header('1_staging', resources['1_staging'])
            + path_block
            + 'rsync -ah --info=progress2 "${SOURCE_DIR}/" "${SCRATCH_DIR}/"\n'
        ),
        '2_renaming': (
            generate_slurm_header('2_renaming', resources['2_renaming'])
            + path_block
            + 'module load go/1.20.1 miniconda3\n'
            + f'source activate "${{QC_ENV_DC}}"\n'
            + 'python "${SCRIPTS_DIR}/_renamer.py"\n'
        ),
        '3_ashlar': (
            generate_slurm_header('3_ashlar', resources['3_ashlar'])
            + path_block
            + get_monitoring_script_block('3_ashlar')
            + '\nexport PYTHONNOUSERSITE=1\n'
            + 'cd "${SCRATCH_DIR}"\n'
            + 'module load go/1.20.1 miniconda3 openjdk/17.0.2\n'
            + f'source activate "${{CONDA_BASE}}/ashlar_group"\n'
            + 'bash ./ashlar.sh\n'
        ),
        '4_preprocess': (
            generate_slurm_header('4_preprocess', resources['4_preprocess'])
            + path_block
            + get_monitoring_script_block('4_preprocess')
            + '\ncd "${SCRATCH_DIR}"\n'
            + 'module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow\n'
            + 'nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-preprocess.yml\n'
        ),
        '5_segment': (
            generate_slurm_header('5_segment', resources['5_segment'])
            + path_block
            + get_monitoring_script_block('5_segment')
            + '\ncd "${SCRATCH_DIR}"\n'
            + 'module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3\n'
            + 'nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-wc.yml\n'
            + f'source activate "${{QC_ENV_DC}}"\n'
            + f'bash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel "{cyto_idx_str}"    --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer{tma_flag}\n'
            + 'conda deactivate\n'
            + 'mv segmentation segmentation_wc; mv qc_overlays_segmentation-mesmer qc_overlays_wc\n'
            + 'nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-nuc.yml\n'
            + f'source activate "${{QC_ENV_DC}}"\n'
            + f'bash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel "{cyto_idx_str}" --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer{tma_flag}\n'
            + 'conda deactivate\n'
            + 'for dir in segmentation/mesmer-*; do mv "$dir/cell.tif" "$dir/nuclear.tif" 2>/dev/null; done\n'
            + 'mv segmentation segmentation_nuc; mv qc_overlays_segmentation-mesmer qc_overlays_nuc\n'
            + 'mkdir -p segmentation\n'
            + 'for dir in segmentation_nuc/mesmer-*; do core_id=$(basename "$dir"); mkdir -p "segmentation/$core_id"; cp "segmentation_wc/$core_id/cell.tif" "segmentation/$core_id/" 2>/dev/null; cp "$dir/nuclear.tif" "segmentation/$core_id/" 2>/dev/null; done\n'
        ),
        '6_quant_merge': (
            generate_slurm_header('6_quant_merge', resources['6_quant_merge'])
            + path_block
            + get_monitoring_script_block('6_quant_merge')
            + '\ncd "${SCRATCH_DIR}"\n'
            + 'module load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3\n'
            + 'nextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-quant.yml\n'
            + 'mkdir -p quantification_wc quantification_nuc\n'
            + 'mv quantification/*cell.csv quantification_wc/\n'
            + 'mv quantification/*nuclear.csv quantification_nuc/\n'
            + f'source activate "${{QC_ENV_RPT}}"\n'
            + f'python ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_wc --output-path qc_report_wc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}\n'
            + f'python ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_nuc --output-path qc_report_nuc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}\n'
            + 'conda deactivate\n'
            + 'mv quantification_wc/* quantification/; mv quantification_nuc/* quantification/\n'
            + 'rmdir quantification_wc quantification_nuc\n'
            + f'source activate "${{QC_ENV_RPT}}"\n'
            + 'python ./scripts/batch_merge_all_cores_v10.py --project-dir .\n'
            + 'conda deactivate\n'
        ),
        '7_archive': (
            generate_slurm_header('7_archive', resources['7_archive'])
            + path_block
            + '\nfor item in registration background dearray segmentation quantification_merged qc_overlays_wc qc_overlays_nuc qc_report_wc.html qc_report_nuc.html QA_QC_Report.txt; do\n'
            + '    rsync -ah --info=progress2 "${SCRATCH_DIR}/$item" "${SOURCE_DIR}/" 2>/dev/null\n'
            + 'done\n'
        ),
    }

    slurm_paths = {}
    if not args.dry_run:
        for k, v in job_scripts.items():
            p = os.path.join(scratch_run_dir, f"{k}.slurm")
            with open(p, 'w') as f: f.write(v)
            slurm_paths[k] = p

    if args.dry_run: sys.exit(0)

    # --- Submit the chain with validation and rollback ---
    JOB_ORDER = ['1_staging', '2_renaming', '3_ashlar', '4_preprocess',
                 '5_segment', '6_quant_merge', '7_archive']
    submitted = submit_chain(slurm_paths, JOB_ORDER)
    print(f"\nSuccess! 7-Job Chain Submitted. IDs: {', '.join(submitted)}")
