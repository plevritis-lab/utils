import argparse
import numpy as np
import os
import pandas as pd
from skimage.measure import regionprops
from tifffile import imread
from utils import extract_proteomic_panel

def quantify_expression(image, segmentation_mask, panel, save_path):
    """quantifies expression of every marker present in the panel in the provided image per cell

    args:
        image (array): loaded matrix representation of the image of shape (c, y, x)
        segmentation_mask (array): loaded matrix representation of the segmentation mask of shape (y, x)
        panel (list): list of protein channel names
        save_path (str): directory where output (cell expressions) will be written
    """

    cell_data = {}

    for channel_number in range(image.shape[0]):
        channel_image = image[channel_number, :, :]

        properties = regionprops(segmentation_mask, channel_image)

        for property in properties:
            cell_identifier = property.label

            cell_major_axis_length = property.axis_major_length
            cell_minor_axis_length = property.axis_minor_length

            cell_centroid_y, cell_centroid_x = property.centroid
            cell_eccentricity = property.eccentricity
            cell_intensity = property.intensity_mean
            cell_size = property.num_pixels
            cell_orientation = property.orientation

            if cell_identifier not in cell_data:
                cell_data[cell_identifier] = {"CELL_IDENTIFIER" : cell_identifier}

            cell_data[cell_identifier]["MAJOR_AXIS_LENGTH"] = cell_major_axis_length
            cell_data[cell_identifier]["MINOR_AXIS_LENGTH"] = cell_minor_axis_length
            
            cell_data[cell_identifier]["X"] = round(cell_centroid_x)
            cell_data[cell_identifier]["Y"] = round(cell_centroid_y)
            cell_data[cell_identifier]["SIZE"] = round(cell_size)
            cell_data[cell_identifier]["ECCENTRICITY"] = cell_eccentricity            
            cell_data[cell_identifier]["ORIENTATION"] = cell_orientation
            
            cell_data[cell_identifier][panel[channel_number].upper()] = cell_intensity
    
    cell_expressions = pd.DataFrame.from_dict(cell_data, orient = "index")
    cell_expressions.to_csv(os.path.join(save_path, "cell_expressions.csv"), index = False)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for quantifying marker expression from generated segmentation masks")

    parser.add_argument("-i", "--image_path", help = "file path that points to the underlying location of the QPTIFF or TIF")
    parser.add_argument("-m", "--mask_path", help = "file path that points to the underlying location of the segmentation mask, stored as a .npy file")
    parser.add_argument("-p", "--panel_path", help = "file path that points to the underlying location of the channel_names.txt file; \
                                                      this argument is optional and should only be used when metadata parsing fails")

    parser.add_argument("--use_mesmer", action = "store_true", help = "toggle only if using a mesmer-generated segmentation mask")
    parser.add_argument("--use_cellpose", action = "store_true", help = "toggle only if using a cellpose-generated segmentation mask")

    return parser.parse_args()

def main():
    arguments = parse_arguments()

    image_path = arguments.image_path
    mask_path = arguments.mask_path
    panel_path = arguments.panel_path

    print("\nprocessing sample", os.path.basename(os.path.splitext(image_path)[0]))

    try:
        panel = extract_proteomic_panel(image_path)
    
    except AttributeError:
        print("\nERROR: unable to extract marker metadata from the provided file; "
                  "please provide a channel_names.txt file instead")
        
        if panel_path:
            print(f"INFO: extracting protein panel from {panel_path}")
            panel = extract_proteomic_panel(image_path, panel_path)
        
        else:
            return

    image = imread(image_path) # (c, y, x)
    
    save_path = os.path.dirname(mask_path)

    if arguments.use_mesmer:
        mesmer_segmentation = np.load(mask_path, allow_pickle = True)
        quantify_expression(image, mesmer_segmentation, panel, save_path)
    
    if arguments.use_cellpose:
        cellpose_segmentation = np.load(mask_path, allow_pickle = True).item()["outlines"]
        quantify_expression(image, cellpose_segmentation, panel, save_path)

if __name__ == "__main__":
    main()