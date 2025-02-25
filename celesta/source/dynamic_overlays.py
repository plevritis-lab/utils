import argparse
import json
import napari
import numpy as np
import os
import pandas as pd
from tifffile import imread

def create_assignment_overlays(spatial_images, assignments, colormap):
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
    for cell_type in unique_cell_types:
        centroids = assignments.loc[(assignments["FINAL_CELL_TYPE"] == cell_type), ["Y", "X"]]
        
        try:
            color = colormap[cell_type]["color"]

        except KeyError:
            print("\nERROR: please provide a color for cell type:", cell_type)
            return

        image_layer = viewer.add_points(centroids,
                                        border_color = "transparent",
                                        face_color = color,
                                        name = cell_type,
                                        size = 15)
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
    
    create_assignment_overlays(spatial_images, assignments, cell_type_info)

main()