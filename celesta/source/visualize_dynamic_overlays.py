import argparse
import json
import napari
import numpy as np
import os
import pandas as pd
from scipy.ndimage import binary_erosion
import seaborn as sns
from tifffile import imread

from widgets.cell_proportions import show_cell_type_proportions

def stratify_marker_probabilities(assignments):
    """stratifies cell marker probabilities into discrete categories (in place) for visualization

    args:
        assignments (dataframe): dataframe containing cell assignments and marker probabilities; \
                                 marker probability columns must be suffixed with "_PROBABILITY"
    """

    def categorize_expression(marker_probability):
        """categorizes marker probabilities into discrete bins

        args:
            marker_probability (float): probability value between 0 and 1
        """

        if marker_probability > 0.9:
            return ">0.9"
        elif marker_probability > 0.8:
            return ">0.8"
        elif marker_probability > 0.7:
            return ">0.7"
        elif marker_probability > 0.5:
            return ">0.5"
        else:
            return "<=0.5"
    
    probability_columns = assignments.columns[assignments.columns.str.endswith("_PROBABILITY")]
    assignments[probability_columns] = assignments[probability_columns].apply(lambda x: \
                                                                                  x.map(categorize_expression))
    
    return assignments

def create_cell_type_outline(assignments, segmentation_mask, cell_type):
    """creates an outline for a specific cell type using a segmentation mask

    args:
        assignments (dataframe): dataframe containing cell assignments and locations
        segmentation_mask (array): mask array where each cell has a unique identifier of shape (y, x)
        cell_type (string): name of the cell type to outline
    """

    cell_type_outline = np.zeros(segmentation_mask.shape)
    cell_ids = assignments.loc[assignments["FINAL_CELL_TYPE"] == cell_type, "CELL_IDENTIFIER"]
    
    for cell_id in cell_ids:
        cell_mask = (segmentation_mask == cell_id)

        eroded = binary_erosion(cell_mask)
        outline = cell_mask ^ eroded
        cell_type_outline[outline] = 1
    
    return cell_type_outline

def add_spatial_images(viewer, spatial_images, panel, markers_to_colors):
    """adds spatial images to the napari viewer with user-provided marker coloring

    args:
        viewer (viewer): napari viewer instance
        spatial_images (dictionary): dictionary mapping modality names to image arrays of shape (c, y, x)
        panel (list): list of marker names
        markers_to_colors (dictionary): mapping of marker names to hexadecimal color values
    """

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

def add_cell_type_visualization(viewer, assignments, colormap, segmentation_mask):
    """adds cell type visualization overlays to the napari viewer

    args:
        viewer (viewer): napari viewer instance
        assignments (dataframe): dataframe containing cell assignments and locations
        colormap (dictionary): dictionary mapping cell types to colors and marker information
        segmentation_mask (array, optional): segmentation mask array for cell outlines of shape (y, x)
    """

    unique_cell_types = assignments["FINAL_CELL_TYPE"].unique()

    if segmentation_mask is not None:
        for cell_type in unique_cell_types:
            try:
                color = colormap[cell_type]["color"]

            except KeyError:
                print("\nERROR: please provide a color for cell type:", cell_type)
                return False

            cell_type_outline = create_cell_type_outline(assignments, segmentation_mask, cell_type)

            image_layer = viewer.add_image(cell_type_outline, name = cell_type, blending = "additive", 
                                               colormap = napari.utils.colormaps.Colormap(["transparent", color]))
            image_layer.visible = False

    else:
        for cell_type in unique_cell_types:
            try:
                color = colormap[cell_type]["color"]

            except KeyError:
                print("\nERROR: please provide a color for cell type:", cell_type)
                return False

            centroids = assignments.loc[(assignments["FINAL_CELL_TYPE"] == cell_type), ["Y", "X"]]

            image_layer = viewer.add_points(centroids, name = cell_type, border_color = "transparent", 
                                                face_color = color, size = 15)
            image_layer.visible = False

    return True

def add_probability_visualization(viewer, assignments, panel, image_shape, legend_spacing = (50, 200)):
    """adds marker probability visualization overlays to the napari viewer as points

    args:
        viewer (viewer): napari viewer instance
        assignments (dataframe): dataframe containing cell assignments and probabilities
        panel (list): list of marker names to visualize
        image_shape (tuple): shape of the proteomic image array of form (y, x)
        legend_spacing (tuple): spacing between legend items of shape (y, x); \
                                y controls the horizontal spacing between the point cloud and legend items; \
                                x controls the vertical spacing between legend items; \
                                defaults to (50, 200)
    """

    probability_colors = sns.color_palette("light:b", 4)
    probability_color_mapping = {
        ">0.5": probability_colors[0],
        ">0.7": probability_colors[1],
        ">0.8": probability_colors[2],
        ">0.9": probability_colors[3]
    }

    centroids = assignments[["Y", "X"]].values

    num_centroids = len(centroids)
    num_markers = len(panel)
    num_legend_items = len(probability_color_mapping)
    num_total_points = num_markers * (num_centroids + num_legend_items + 1)

    point_cloud_sequence = np.zeros((num_total_points, 3))
    color_sequence = np.zeros((num_total_points, 3))
    point_sizes = np.full(num_total_points, 15)
    text_labels = [""] * num_total_points

    for i, marker in enumerate(panel):
        start_index = i * num_centroids
        end_index = (i + 1) * num_centroids

        point_cloud_sequence[start_index : end_index, 0] = i
        point_cloud_sequence[start_index : end_index, 1:] = centroids
        
        marker_probabilities = assignments[f"{marker}_PROBABILITY"].values
        for discretization, rgb_color in probability_color_mapping.items():
            color_sequence[start_index : end_index][marker_probabilities == discretization] = rgb_color

        title_index = num_centroids * num_markers + i
    
        point_cloud_sequence[title_index] = [i, image_shape[0] / 2, image_shape[1] + legend_spacing[1]]
        text_labels[title_index] = marker

        legend_start_index = num_markers * (num_centroids + 1) + (i * num_legend_items)
        for j, (probability_label, color) in enumerate(probability_color_mapping.items()):
            legend_end_index = legend_start_index + j
            vertical_offset = legend_spacing[0] + (j * legend_spacing[0])

            point_cloud_sequence[legend_end_index] = [i, image_shape[0] / 2 + vertical_offset, image_shape[1] + legend_spacing[1]]
            color_sequence[legend_end_index] = color
            point_sizes[legend_end_index] = 20
            text_labels[legend_end_index] = probability_label

    image_layer = viewer.add_points(point_cloud_sequence, name = "probabilities", border_color = "transparent",
                                        face_color = color_sequence, size = point_sizes, symbol = "x", text = {
                                            "string": text_labels, 
                                            "color": "white", 
                                            "size": 15,
                                            "translation": np.array([0, 3, 50])
                                        })
    image_layer.visible = False

def create_assignment_overlays(spatial_images, assignments, colormap, segmentation_mask = None):
    """creates a napari viewer with various visual overlays for cell assignments

    args:
        spatial_images (dictionary): dictionary mapping modality names ("proteomic", "histology") to image arrays of shape (c, y, x)
        assignments (dataframe): dataframe containing cell assignments, locations, and marker probabilities
        colormap (dictionary): dictionary mapping cell types to colors and marker information
        segmentation_mask (array): segmentation mask array of shape (y, x) where each cell has unique identifier; \
                                   defaults to None
    """

    viewer = napari.Viewer()

    exclude_columns = ["CELL_IDENTIFIER", "MAJOR_AXIS_LENGTH", "MINOR_AXIS_LENGTH", "X", "Y", 
                           "SIZE", "ECCENTRICITY", "ORIENTATION", "CELL_TYPE_NUMBER", "FINAL_CELL_TYPE"]
    panel = assignments.columns[~assignments.columns.str.contains("_PROBABILITY$") & 
                                    ~assignments.columns.isin(exclude_columns)].tolist()
    
    markers_to_colors = {info["marker"]: info["color"] for info in colormap.values() if info.get("marker")}
    
    add_spatial_images(viewer, spatial_images, panel, markers_to_colors)
    add_cell_type_visualization(viewer, assignments, colormap, segmentation_mask)
    add_probability_visualization(viewer, stratify_marker_probabilities(assignments), panel, 
                                      spatial_images["proteomic"][0].shape)
    
    show_cell_type_proportions(viewer, assignments)

    napari.run()

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for dynamic celesta assignment visualization through napari")

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

    assignments = pd.read_csv(assignments_path).dropna()

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