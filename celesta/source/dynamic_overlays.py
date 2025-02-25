import argparse
import json
import napari
import numpy as np
import os
import pandas as pd
from scipy.ndimage import binary_erosion
from tifffile import imread

def create_assignment_overlays(spatial_images, assignments, colormap, segmentation_mask = None):
    """creates interactive visualization overlays for cell assignments and spatial images using napari

    args:
        spatial_images (dictionary): dictionary containing image arrays for proteomic and / or histology data
        assignments (dataframe): dataframe containing cell assignments and their properties
        colormap (dictionary): dictionary mapping cell types and markers to their display colors
        segmentation_mask (array, optional): array containing cell segmentation masks; \
                                             defaults to None
    """

    viewer = napari.Viewer()

    exclude_columns = ["CELL_IDENTIFIER", "MAJOR_AXIS_LENGTH", "MINOR_AXIS_LENGTH", "X", "Y", 
                           "SIZE", "ECCENTRICITY", "ORIENTATION", "CELL_TYPE_NUMBER", "FINAL_CELL_TYPE"]
    panel = assignments.columns[~assignments.columns.str.contains("_PROBABILITY$") & 
                                    ~assignments.columns.isin(exclude_columns)].tolist()
    
    markers_to_colors = {info["marker"]: info["color"] for info in colormap.values() if info.get("marker")}
    
    for modality, image in spatial_images.items():
        spatial_images[modality] = (image - image.min()) / (image.max() - image.min())

        if modality == "proteomic":
            for i, protein in enumerate(panel):
                image_layer = viewer.add_image(spatial_images[modality][i], name = protein,
                                                   contrast_limits = [0, 1], blending = "additive")
                image_layer.visible = False
                
                if protein in markers_to_colors:
                    image_layer.colormap = napari.utils.colormaps.Colormap(["black", 
                                               markers_to_colors[protein]])

        elif modality == "histology":
            image_layer = viewer.add_image(np.transpose(spatial_images[modality], (1, 2, 0)), name = modality, 
                                               contrast_limits = [0, 1], blending = "additive")
            image_layer.visible = False

    unique_cell_types = assignments["FINAL_CELL_TYPE"].unique()

    if segmentation_mask is not None:
        for cell_type in unique_cell_types:
            try:
                color = colormap[cell_type]["color"]

            except KeyError:
                print("\nERROR: please provide a color for cell type:", cell_type)
                return
                
            cell_type_outline = np.zeros(segmentation_mask.shape)
            cell_ids = assignments.loc[assignments["FINAL_CELL_TYPE"] == cell_type, "CELL_IDENTIFIER"]
            
            for cell_id in cell_ids:
                cell_mask = (segmentation_mask == cell_id)

                eroded = binary_erosion(cell_mask)
                outline = cell_mask ^ eroded
                cell_type_outline[outline] = 1
            
            image_layer = viewer.add_image(cell_type_outline, name = cell_type, blending = "additive",
                                               colormap = napari.utils.colormaps.Colormap(["transparent", color]))
            image_layer.visible = False
        
    else:
        for cell_type in unique_cell_types:
            centroids = assignments.loc[(assignments["FINAL_CELL_TYPE"] == cell_type), ["Y", "X"]]
            
            try:
                color = colormap[cell_type]["color"]

            except KeyError:
                print("\nERROR: please provide a color for cell type:", cell_type)
                return

            image_layer = viewer.add_points(centroids, name = cell_type, border_color = "transparent", 
                                                face_color = color, size = 15)
            image_layer.visible = False
    
    napari.run()

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for celesta assignment visualization through napari")

    parser.add_argument("-a", "--assignments_path", help = "file path that points to the underlying location of the assignments.csv file")
    parser.add_argument("-c", "--colormap_path", help = "file path that points to the underlying location of the colormap.json file")
    parser.add_argument("-i", "--image_path", help = "file path that points to the underlying location of the proteomic TIF or TIFF file")

    parser.add_argument("-e", "--histology_path", help = "(optional) file path that points to the underlying location of the histology TIF or TIFF file")
    parser.add_argument("-m", "--mask_path", help = "(optional) file path that points to the underlying location of the segmentation mask, stored as a .npy file")

    parser.add_argument("--use_mesmer", action = "store_true", help = "toggle only if using a mesmer-generated segmentation mask")
    parser.add_argument("--use_cellpose", action = "store_true", help = "toggle only if using a cellpose-generated segmentation mask")

    return parser.parse_args()
    
def main():
    arguments = parse_arguments()

    assignments_path = arguments.assignments_path
    colormap_path = arguments.colormap_path
    image_path = arguments.image_path

    histology_path = arguments.histology_path
    mask_path = arguments.mask_path

    print("\nvisualizing sample", os.path.basename(os.path.splitext(image_path)[0]))

    assignments = pd.read_csv(assignments_path)

    with open(colormap_path, "r") as file:
        cell_type_info = json.load(file)
    
    spatial_images = {"proteomic": imread(image_path), 
                      "histology": imread(histology_path)} if histology_path \
                                                               else {"proteomic": imread(image_path)}
    
    segmentation_mask = None

    if arguments.use_mesmer:
        segmentation_mask = np.load(mask_path, allow_pickle = True)

    if arguments.use_cellpose:
        segmentation_mask = np.load(mask_path, allow_pickle = True).item()["outlines"]
    
    create_assignment_overlays(spatial_images, assignments, cell_type_info, segmentation_mask)

main()