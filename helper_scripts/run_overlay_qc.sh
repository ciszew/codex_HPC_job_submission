#!/bin/bash

# =============================================================================
# MCMICRO Segmentation QC Overlay Generator (Dynamic Launcher V4)
# =============================================================================
#
# Purpose:
#   This script is a generic tool to automate the creation of QC overlays. It
#   receives all configuration via command-line arguments, making it reusable
#   for any TMA or whole-slide workflow without modification.
#
# =============================================================================

# --- ARGUMENT PARSING ---
IS_TMA=false
CYTO_CHANNEL_CMD=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --base-dir) BASE_DIR="$2"; shift ;;
        --nuc-channel) NUC_CHANNEL="$2"; shift ;;
        --cyto-channel) CYTO_CHANNEL_CMD="--cyto-channel-index $2"; shift ;;
        --overlay-script) OVERLAY_SCRIPT="$2"; shift ;;
        --seg-base) SEG_FOLDER_BASE="$2"; shift ;;
        --is-tma) IS_TMA=true ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# --- WORKFLOW EXECUTION ---

if [ "$IS_TMA" = true ]; then
    echo "--- Starting TMA Overlay Generation ---"
    DEARRAY_DIR="${BASE_DIR}/dearray"
    OUTPUT_DIR="${BASE_DIR}/qc_overlays_${SEG_FOLDER_BASE//\//-}" 

    mkdir -p "$OUTPUT_DIR"

    for image_file in ${DEARRAY_DIR}/*.tif*; do
        base_name=$(basename "$image_file")
        core_id=$(echo "$base_name" | cut -d'.' -f1)
        
        # TMA Mesmer mask path
        mask_file="${BASE_DIR}/${SEG_FOLDER_BASE}-${core_id}/cell.tif"
        output_file="${OUTPUT_DIR}/${core_id}_overlay.ome.tif"
        
        if [ -f "$mask_file" ]; then
            # CYTO_CHANNEL_CMD is unquoted here so Python receives it as distinct arguments if there are multiple indices
            CMD="python $OVERLAY_SCRIPT --image $image_file --mask $mask_file --output $output_file --nuclear-channel-index $NUC_CHANNEL $CYTO_CHANNEL_CMD"
            echo "Executing: $CMD"
            eval "$CMD"
        else
            echo "Warning: Mask not found for core ${core_id} at expected path ${mask_file}. Skipping."
        fi
    done
    echo "--- TMA Overlay Generation Complete ---"
    
else
    echo "--- Starting Whole Slide Overlay Generation ---"
    
    # 1. Find the stitched image
    IMAGE_FILE=$(find "${BASE_DIR}/background" -name "*.ome.tif" -print -quit)
    
    # 2. Dynamically resolve the image base name (e.g., 'stitched_backsub')
    BASE_NAME=$(basename "$IMAGE_FILE" .ome.tif)
    
    # 3. Construct the exact path MCMICRO uses for WSI masks
    MASK_FILE="${BASE_DIR}/${SEG_FOLDER_BASE}-${BASE_NAME}/cell.tif"
    
    OUTPUT_DIR="${BASE_DIR}/qc_overlays_${SEG_FOLDER_BASE//\//-}"
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}_overlay.ome.tif"

    mkdir -p "$OUTPUT_DIR"
    
    if [ -f "$IMAGE_FILE" ] && [ -f "$MASK_FILE" ]; then
        CMD="python $OVERLAY_SCRIPT --image $IMAGE_FILE --mask $MASK_FILE --output $OUTPUT_FILE --nuclear-channel-index $NUC_CHANNEL $CYTO_CHANNEL_CMD"
        echo "Executing: $CMD"
        eval "$CMD"
    else
        echo "Warning: Could not find required files for whole-slide overlay."
        echo "Image: $IMAGE_FILE"
        echo "Mask: $MASK_FILE"
    fi
    echo "--- Whole Slide Overlay Generation Complete ---"
fi