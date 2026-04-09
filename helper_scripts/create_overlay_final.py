import argparse
import numpy as np
from tifffile import imread, imwrite
import os
import sys

# This function is only imported if deepcell is installed.
from deepcell.utils.plot_utils import make_outline_overlay

def create_qc_overlay(image_path, mask_path, output_path, nuclear_idx, cyto_idx=None):
    """
    Creates a multi-channel QC TIFF for a single image/mask pair.
    """
    print(f"\n--- Processing Image: {os.path.basename(image_path)} ---")
    try:
        image_stack = imread(image_path)
        mask_img = imread(mask_path)
        nuclear_img = image_stack[nuclear_idx, :, :]
        
        cyto_img = None
        if cyto_idx is not None:
            cyto_img = image_stack[cyto_idx, :, :]

        def normalize(img_array):
            img_min, img_max = np.min(img_array), np.max(img_array)
            if img_max > img_min:
                return ((img_array - img_min) / (img_max - img_min) * 255).astype(np.uint8)
            return np.zeros_like(img_array, dtype=np.uint8)

        nuclear_norm_8bit = normalize(nuclear_img)
        rgb_image = np.zeros((*nuclear_img.shape, 3), dtype=np.uint8)
        output_stack_channels = [nuclear_img]
        
        if cyto_img is not None:
            cyto_norm_8bit = normalize(cyto_img)
            rgb_image[..., 1] = cyto_norm_8bit
            rgb_image[..., 2] = nuclear_norm_8bit
            output_stack_channels.append(cyto_img)
        else:
            rgb_image[..., 2] = nuclear_norm_8bit

        rgb_image_4d = np.expand_dims(rgb_image, axis=0)
        mask_img_4d = mask_img[np.newaxis, ..., np.newaxis]
        
        overlay_4d = make_outline_overlay(rgb_data=rgb_image_4d, predictions=mask_img_4d)
        overlay_rgb = overlay_4d[0]

        is_outline = (overlay_rgb[..., 0] > 200) & (overlay_rgb[..., 1] > 200)
        overlay_rgb[is_outline] = [255, 255, 255]
        
        overlay_r, overlay_g, overlay_b = overlay_rgb[..., 0], overlay_rgb[..., 1], overlay_rgb[..., 2]
        final_stack_list = [overlay_r, overlay_g, overlay_b] + output_stack_channels
        
        imwrite(output_path, np.stack(final_stack_list, axis=0), imagej=True)
        print(f"Successfully saved {len(final_stack_list)}-channel overlay to: {output_path}")

    except Exception as e:
        print(f"Could not process {os.path.basename(image_path)}. Error: {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a flexible, multi-channel QC overlay TIFF.')
    parser.add_argument('--image', required=True, help='Path to the input multi-channel TIFF file.')
    parser.add_argument('--mask', required=True, help='Path to the segmentation label mask.')
    parser.add_argument('--output', required=True, help='Path for the output multi-channel TIFF file.')
    parser.add_argument('--nuclear-channel-index', type=int, required=True, help='0-based index of the nuclear channel.')
    parser.add_argument('--cyto-channel-index', type=int, help='(Optional) 0-based index of the cytoplasm channel.')
    args = parser.parse_args()
    create_qc_overlay(
        args.image, args.mask, args.output,
        args.nuclear_channel_index, args.cyto_channel_index
    )
