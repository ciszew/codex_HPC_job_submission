# Standard Operating Procedure: MCMICRO SLURM Automation Pipeline (V048)

## 1. Introduction
The `create_codex_slurm_submission_V048.py` script is a master automation wrapper designed to deploy modular components (Ashlar, Background Subtraction, Deepcell Segmentation, and Quantification) , some of them part of the MCMICRO pipeline on a SLURM-managed High-Performance Computing (HPC) cluster Randi provided by CRI at University of Chicago. This pipeline was created for the Codex spatial proteomics platform at HDID core facility at the University of Chicago (https://voices.uchicago.edu/hdid/codex-platform/).

Version 48 introduces a major architectural upgrade: **Dynamic Marker Name Resolution**. Instead of relying on hardcoded channel indices—which frequently break when upstream channels are dropped—the script accepts exact marker names (e.g., `UV_high`, `CD45_Atto550`). It simulates both the Ashlar drop phase (`ashlar=skip`) and the Background Subtraction drop phase (`remove=TRUE`) to dynamically calculate the precise 0-based indices for downstream segmentation and QC overlay generation.

The script provisions a highly robust, seven-job dependency chain that handles everything from initial data staging to final archival, minimizing manual intervention and ensuring safe, sequential execution.

---

## 2. Dependencies and Requirements

Before executing the master script, ensure the following environment constraints and file dependencies are met.

### A. Input Data Requirements
The target raw data directory must contain:
* **Raw Image Tiles:** `.ims` or `.ome.tif` files organized into per cycle subdirectories.
* **XML Metadata:** At least one `.xml` file containing voxel dimensions and stack layouts.
* **Marker Configuration (`markers_ashlar.csv`):** This file is strictly required in the root of the data directory. It must contain `marker_name`, `ashlar`, and `remove` columns to dictate channel inclusion.

### B. HPC Environment & Helper Scripts
The script expects specific Python scripts and Conda environments to exist on the shared HPC drive:
* **Conda Environments:** 
  * `ashlar_group` (Stitching)
  * `deepcell` (QC Overlays)
  * `qc` (HTML Reports)
* **QC Helper Scripts (`QC_SCRIPTS_DIR`):**
  * `run_overlay_qc.sh`: Automates mask/image overlays.
  * `create_overlay_final.py`: The python logic for generating the OME-TIFF overlays.
  * `mcmicro_qc_analyzer.py`: Generates the interactive HTML QC reports.
  * `batch_merge_all_cores_v10.py`: Concatenates all core-level CSVs into project-level summary tables.

---

## 3. How to Run the Master Script

### Step 1: Configure Pipeline Variables
Open `create_codex_slurm_submission_V048.py` in a text editor and adjust the `USER-CONFIGURABLE VARIABLES` block. 
* Toggle `IS_TMA_WORKFLOW` (`True` for TMAs, `False` for Whole Slide).
* Set the `NUC_MARKER_NAME` and `MEMBRANE_MARKER_NAMES` to exactly match the strings in your `markers_ashlar.csv`.

### Step 2: Dry-Run (Recommended)
Always run a pre-flight check to ensure the dynamic index mapping resolves correctly and scripts generate without syntax errors.

```bash
python3 create_codex_slurm_submission_V048.py /path/to/raw_data_directory --dry-run
```
*Note: This will print the calculated resource allocations and marker mappings to the console without submitting anything to the SLURM queue.*

### Step 3: Execute Full Pipeline
Submit the master script. It will generate all necessary YAMLs, Python wrappers, and SLURM `.slurm` scripts in your scratch directory, and submit the 7-job dependency chain.

```bash
python3 create_codex_slurm_submission_V048.py /path/to/raw_data_directory
```

---

## 4. Detailed Job Functionality (The 7-Job Chain)

The pipeline submits seven dependent SLURM jobs. If one fails, the subsequent jobs remain in a pending state or abort, preserving the data state.

### Job 1: Staging (`1_staging`)
* **Function:** Safely transfers raw data from the persistent storage drive to the high-speed scratch directory using `rsync`.
* **Benefit:** Prevents network bottlenecking during heavy I/O operations (like Ashlar stitching).

### Job 2: Renaming & Integrity Check (`2_renaming`)
* **Function:** Sorts raw FOV directories chronologically using date/time stamps parsed from folder names, renaming them to a clean sequential format (e.g., `Core_Run_1_Field_1`).
* **QA/QC Check:** Scans all `.ims` files and compares their byte size against the cycle median. If a file is significantly smaller (corrupted), it automatically patches the missing data by duplicating a healthy adjacent tile, preventing Ashlar from crashing.

### Job 3: Stitching (`3_ashlar`)
* **Function:** Dynamically generates an `ashlar_wrapper.py` script. This wrapper intercepts the standard Ashlar fileseries reader, hiding any channels marked `skip` in `markers_ashlar.csv`. 
* **Output:** Produces a fully registered `stitched.ome.tiff` and writes a cleaned `markers.csv` (with updated, contiguous channel numbers) for downstream MCMICRO modules.

### Job 4: Preprocessing (`4_preprocess`)
* **Function:** Invokes the Nextflow MCMICRO pipeline for early-stage processing. 
* **TMA vs. WSI Logic:** If `IS_TMA_WORKFLOW` is True, it runs background subtraction followed by `coreograph` (dearraying). If False, it stops after background subtraction.

### Job 5: Dual-Compartment Segmentation (`5_segment`)
* **Function:** Orchestrates two separate runs of the Deepcell/Mesmer segmentation module.
  * **Pass 1 (Whole-Cell):** Segments the whole cell. Generates deepcell QC overlays using the dynamically calculated membrane and nuclear channel indices.
  * **Pass 2 (Nuclear):** Segments only the nucleus. Generates nuclear QC overlays. Renames the output mask to `nuclear.tif`.
* **Consolidation:** Merges both `cell.tif` and `nuclear.tif` masks into a unified final segmentation directory.

### Job 6: Quantification & QC Reporting (`6_quant_merge`)
* **Function:** Runs `mcquant` on both the whole-cell and nuclear masks to extract single-cell spatial and intensity metrics.
* **QC Generation:** Separates the outputs and passes them to `mcmicro_qc_analyzer.py`. This generates two interactive HTML reports (`qc_report_wc.html` and `qc_report_nuc.html`) containing marker correlations, arcsinh-transformed intensity distributions, and cell shape analysis.
* **Merging:** Runs the batch merge script to concatenate all core-level data into a single master project CSV.

### Job 7: Archive (`7_archive`)
* **Function:** The cleanup phase. Uses `rsync` to push all completed pipeline artifacts (registration files, background subtracted images, dearrayed TIFFs, segmentation masks, merged quantification CSVs, QC overlays, and HTML reports) back to the original source data directory.