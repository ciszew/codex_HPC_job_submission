#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
MCMICRO Pipeline Submission Script for Slurm HPC
================================================================================
Version: 48 (Dynamic Marker Name Resolution)

Purpose:
  Automates the MCMICRO pipeline. Users specify target channel NAMES instead 
  of indices. The script calculates the exact final indices by simulating 
  the Ashlar drop (ashlar=skip) and Backsub drop (remove=TRUE).
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
SCRATCH_BASE_DIR = "/scratch/cciszews/nextflow_runs/"
CONDA_BASE_PATH = "/gpfs/data/hdid-share/conda"
QC_SCRIPTS_DIR = "/gpfs/data/hdid-share/Codex/HDID/scripts/current_working_scripts/"

IS_TMA_WORKFLOW = False
IMAGE_MPP = 0.32
ARCSINH_COFACTOR = 5

# Set Segmentation Channels by EXACT MARKER NAME (as written in markers_ashlar.csv)
NUC_MARKER_NAME = "UV_high"
MEMBRANE_MARKER_NAMES = "CD45_Atto550 PanCK_AF750"

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
    days, rem_hours = divmod(int(hours), 24)
    return f"{days}-{rem_hours:02d}:00:00" if days > 0 else f"{int(hours):02d}:00:00"

def estimate_resources(total_gb, num_cores, is_tma):
    ref_gb, ref_ashlar_hours = 326.0, 10.25
    ashlar_hours_per_gb = ref_ashlar_hours / ref_gb
    est_ashlar_hours = (total_gb * ashlar_hours_per_gb) * 1.5 + 2
    
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
        
    return {
        "1_staging":     {"time": "04:00:00", "mem": "8G", "cpu": "1"},
        "2_renaming":    {"time": "01:00:00", "mem": "8G", "cpu": "1"},
        "3_ashlar":      {"time": format_slurm_time(est_ashlar_hours), "mem": "64G", "cpu": "4"},
        "4_preprocess":  {"time": "04:00:00", "mem": "48G", "cpu": "4"},
        "5_segment":     {"time": format_slurm_time(est_segmentation_hours), "mem": seg_mem, "cpu": "4"},
        "6_quant_merge": {"time": format_slurm_time(est_quant_hours), "mem": "64G", "cpu": "4"},
        "7_archive":     {"time": "04:00:00", "mem": "8G", "cpu": "1"}
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

def generate_rename_script(path, dirs_to_rename_map, xml_path):
    return f'''
import os, re, statistics, shutil

path = "{path}"
dirs_map = {dirs_to_rename_map}

for old, new in dirs_map.items():
    try: os.rename(os.path.join(path, old), os.path.join(path, new))
    except OSError: pass

all_cycles = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and re.search(r'_\d+_Field', d)]
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
    
    # Format strictly as a YAML list: [1, 13, 15]
    yaml_list = "[" + ", ".join(all_channels) + "]"

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
    
    print(f"\n--- MCMICRO SUBMISSION V48 ---\nScratch Dir: {scratch_run_dir}")
    
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
    cyto_idx = memb_indices[0] if len(memb_indices) > 0 else nuc_idx

    job_scripts = {
        '1_staging': generate_slurm_header('1_staging', resources['1_staging']) + f"rsync -ah --info=progress2 {source_data_dir}/ {scratch_run_dir}/",
        '2_renaming': generate_slurm_header('2_renaming', resources['2_renaming']) + f"module load go/1.20.1 miniconda3\npython {staged_scripts_dir}/_renamer.py",
        '3_ashlar': generate_slurm_header('3_ashlar', resources['3_ashlar']) + get_monitoring_script_block('3_ashlar') + f"\nexport PYTHONNOUSERSITE=1\ncd {scratch_run_dir}\nmodule load go/1.20.1 miniconda3 openjdk/17.0.2\nsource activate {CONDA_BASE_PATH}/ashlar_group\nbash ./ashlar.sh",
        '4_preprocess': generate_slurm_header('4_preprocess', resources['4_preprocess']) + get_monitoring_script_block('4_preprocess') + f"\ncd {scratch_run_dir}\nmodule load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow\nnextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-preprocess.yml",
        '5_segment': generate_slurm_header('5_segment', resources['5_segment']) + get_monitoring_script_block('5_segment') + f"\ncd {scratch_run_dir}\nmodule load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3\nnextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-wc.yml\nsource activate {QC_ENV_DEEPCELL}\nbash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel {cyto_idx} --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer --is-tma\nconda deactivate\nmv segmentation segmentation_wc; mv qc_overlays_segmentation-mesmer qc_overlays_wc\nnextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-nuc.yml\nsource activate {QC_ENV_DEEPCELL}\nbash ./scripts/run_overlay_qc.sh --base-dir . --nuc-channel {nuc_idx} --cyto-channel {cyto_idx} --overlay-script ./scripts/create_overlay_final.py --seg-base segmentation/mesmer --is-tma\nconda deactivate\nfor dir in segmentation/mesmer-*; do mv \"$dir/cell.tif\" \"$dir/nuclear.tif\" 2>/dev/null; done\nmv segmentation segmentation_nuc; mv qc_overlays_segmentation-mesmer qc_overlays_nuc\nmkdir -p segmentation\nfor dir in segmentation_nuc/mesmer-*; do core_id=$(basename $dir); mkdir -p segmentation/$core_id; cp segmentation_wc/$core_id/cell.tif segmentation/$core_id/ 2>/dev/null; cp $dir/nuclear.tif segmentation/$core_id/ 2>/dev/null; done\n",
        '6_quant_merge': generate_slurm_header('6_quant_merge', resources['6_quant_merge']) + get_monitoring_script_block('6_quant_merge') + f"\ncd {scratch_run_dir}\nmodule load go/1.20.1 openjdk/17.0.2 singularity/3.8.7 nextflow miniconda3\nnextflow run labsyspharm/mcmicro --in . -profile singularity -params-file params-quant.yml\nmkdir -p quantification_wc quantification_nuc\nmv quantification/*cell.csv quantification_wc/\nmv quantification/*nuclear.csv quantification_nuc/\nsource activate {QC_ENV_REPORTS}\npython ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_wc --output-path qc_report_wc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}\npython ./scripts/mcmicro_qc_analyzer.py --input-dir ./quantification_nuc --output-path qc_report_nuc.html --pixel-size {IMAGE_MPP} --cofactor {ARCSINH_COFACTOR}\nconda deactivate\nmv quantification_wc/* quantification/; mv quantification_nuc/* quantification/\nrmdir quantification_wc quantification_nuc\nsource activate {QC_ENV_REPORTS}\npython ./scripts/batch_merge_all_cores_v10.py --project-dir .\nconda deactivate\n",
        '7_archive': generate_slurm_header('7_archive', resources['7_archive']) + f"\nfor item in registration background dearray segmentation quantification_merged qc_overlays_wc qc_overlays_nuc qc_report_wc.html qc_report_nuc.html QA_QC_Report.txt; do rsync -ah --info=progress2 {scratch_run_dir}/$item {source_data_dir}/ 2>/dev/null; done\n"
    }
    
    slurm_paths = {}
    if not args.dry_run:
        for k, v in job_scripts.items():
            p = os.path.join(scratch_run_dir, f"{k}.slurm")
            with open(p, 'w') as f: f.write(v)
            slurm_paths[k] = p

    if args.dry_run: sys.exit(0)

    try:
        j1 = subprocess.run(["sbatch", "--parsable", slurm_paths['1_staging']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        j2 = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{j1}", slurm_paths['2_renaming']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        j3 = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{j2}", slurm_paths['3_ashlar']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        j4 = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{j3}", slurm_paths['4_preprocess']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        j5 = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{j4}", slurm_paths['5_segment']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        j6 = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{j5}", slurm_paths['6_quant_merge']], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        subprocess.run(["sbatch", f"--dependency=afterok:{j6}", slurm_paths['7_archive']], stdout=subprocess.PIPE, universal_newlines=True)
        print("Success! 7-Job Chain Submitted.")
    except Exception as e: 
        print(f"Error submitting jobs: {e}")