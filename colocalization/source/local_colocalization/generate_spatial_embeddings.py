import argparse
import glob
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from scipy.spatial.distance import cdist

BASE_SIZE = 20
POINT_SIZE = 8

# TODO - convert the LCLQ calculation in terms of MICRONS instead of PIXELS; convert pixels to microns for more accurate results

def compute_local_colocation_quotient(assignments, A, B, bandwidth, log_transform):
    """calculates a local indicator of the colocation quotient (LCLQ) between two cell types in a single image, which can serve as a per-cell feature embedding

    args:
        assignments (dataframe): generated cell type assignments from CELESTA
        A (str): first cell type of interest (e.g., epithelial)
        B (str): second cell type of interest (e.g., granulocyte)
        bandwidth (int): bandwidth parameter used to assign geospatial weights in the gaussian kernel density function; \
                         a smaller bandwidth typically produces smaller and spiky clusters whereas a larger bandwidth results in smoother and larger clusters
        log_transform (bool): whether to log-transform the colocalization statistics
    """

    A_cells = assignments[assignments["FINAL_CELL_TYPE"] == A]
    B_cells = assignments[assignments["FINAL_CELL_TYPE"] == B]

    N = assignments.shape[0]
    N_A = A_cells.shape[0]
    N_B = B_cells.shape[0]

    # TODO - add in other weighting mechanisms introduced by the paper (box kernel density, topological distance, etc.)

    distances = cdist(A_cells[["X", "Y"]].values, assignments[["X", "Y"]].values) # N_A x N
    gaussian_weights = np.exp(-0.5 * (distances ** 2 / bandwidth ** 2)) # N_A x N
    gaussian_weights[distances > bandwidth] = 0

    mask = np.zeros(N)
    mask[B_cells.index] = 1

    numerator = gaussian_weights.dot(mask) - gaussian_weights[np.arange(N_A), A_cells.index] * mask[A_cells.index]
    denominator = gaussian_weights.sum(axis = 1) - gaussian_weights[np.arange(N_A), A_cells.index]

    # TODO - handle case where there are no cells in bandwidth (divide by zero error in denominator, assign LCLQ value of 0? right now does nan, but make sure to suppress warning in console)

    N_AB = numerator / denominator

    P = N_B / (N - 1) if A != B else (N_A - 1) / (N - 1)

    colocalization_embedding = N_AB / P if P != 0 else np.full(N_A, np.nan)
    if log_transform:
        colocalization_embedding = np.log1p(colocalization_embedding)

    return colocalization_embedding

def generate_spatial_cell_embedding_matrix(assignments, save_path, sample_name, cell_types, log_transform, bandwidth):
    """constructs and saves a spatial embedding matrix for a single image, representing each cell as a function of their LCLQ value with respect to all other cell types

    args:
        assignments (dataframe): generated cell type assignments from CELESTA
        save_path (str): directory where LCLQ embeddings will be stored
        sample_name (str): name of the sample being processed
        cell_types (list): list of cell types to compute LCLQ embeddings for; \
                           if entries are not present in the assignments dataframe, its LCLQ value will be set to NA
        log_transform (bool): whether to log-transform the colocalization statistics
        bandwidth (int, optional): bandwidth parameter used to assign geospatial weights in the gaussian kernel density function; \
                                   a smaller bandwidth typically produces smaller and spiky clusters whereas a larger bandwidth results in smoother and larger clusters
    """

    unique_cell_types = [t for t in (["unknown"] + cell_types) if t in assignments["FINAL_CELL_TYPE"].unique()]

    composite_embedding_matrix = [None] * len(unique_cell_types)

    for i, A in enumerate(unique_cell_types):
        A_cells = assignments[assignments["FINAL_CELL_TYPE"] == A]

        A_embeddings = {B : compute_local_colocation_quotient(assignments, A, B, bandwidth, log_transform) for B in (["unknown"] + cell_types)}
        A_embedding_matrix = pd.DataFrame({
            "FINAL_CELL_TYPE" : A,
            "X" : A_cells["X"].values, 
            "Y" : A_cells["Y"].values,
            **A_embeddings
        })

        composite_embedding_matrix[i] = A_embedding_matrix
    
    transformation = "transformed" if log_transform else "untransformed"
            
    composite_embedding_matrix = pd.concat(composite_embedding_matrix, ignore_index = True)
    composite_embedding_matrix.to_csv(os.path.join(save_path, f"{sample_name}_local_colocalization_matrix_{transformation}_{bandwidth}_pixel_bandwidth.csv"), index = False, na_rep = "nan")

    return composite_embedding_matrix

def visualize_bandwidth_parameter(centroids, bandwidth, save_path):
    """plots and saves a visual representation of the bandwidth parameter (using Euclidean distances) on a randomly chosen point from user-provided centroids

    args:
        centroids (array): array of centroids of shape N x 2 where N is the total number of data points
        bandwidth (int): bandwidth parameter used to control the radius of the spatial neighborhood around a point
        save_path (str): sample directory where a visual representation of the bandwidth parameter will be stored
    """

    matplotlib.use("agg")

    fig, ax = plt.subplots(figsize = (10, 8))

    selected_cell = centroids[np.random.choice(centroids.shape[0])]
    ax.scatter(centroids[:, 0], centroids[:, 1], color = "white")
    ax.scatter(selected_cell[0], selected_cell[1], color = "red")    

    spatial_neighborhood = plt.Circle((selected_cell[0], selected_cell[1]), bandwidth, color = "red", 
                                           linestyle = "--", linewidth = 2, fill = False)
    plt.gca().add_patch(spatial_neighborhood)

    ax.set_facecolor("black")

    fig.savefig(os.path.join(save_path, f"{bandwidth}_bandwidth_visualization.pdf"))
    plt.close(fig)

def plot_image_pairwise_colocalizations(cell_embedding_matrix, A, B, save_path, cell_type_colors, bins = 4):
    """plots and saves the euclidean-based LCLQ measurements between two cell types for a single image, with A's measurements with respect to B graded into <= the number of user-specified number of categories

    args:
        cell_embedding_matrix (dataframe): generated LCLQ measurements for every cell in the image
        A (str): first cell type of interest (e.g., epithelial)
        B (str): second cell type of interest (e.g., granulocyte)
        save_path (str): sample directory where pairwise LCLQ plots will be stored under the reference cell type
        cell_type_colors (dictionary): dictionary mapping cell types to their colors
        bins (int, optional): used to determine the maximum number of LCLQ categories to display; \
                              defaults to 4
    """

    matplotlib.use("agg")

    x_range = cell_embedding_matrix["X"].max() - cell_embedding_matrix["X"].min()
    y_range = cell_embedding_matrix["Y"].max() - cell_embedding_matrix["Y"].min()

    aspect_ratio = x_range / y_range
    
    if aspect_ratio > 1:
        width = BASE_SIZE
        height = BASE_SIZE / aspect_ratio
    else:
        width = BASE_SIZE * aspect_ratio
        height = BASE_SIZE

    fig, (ax, legend_ax) = plt.subplots(
        nrows = 1, ncols = 2, 
        figsize = (width + 1.5, height), 
        gridspec_kw = {"width_ratios": [width, 0.2]},
        facecolor = "black"
    )

    ax.invert_yaxis()
    ax.set_axis_off()
    legend_ax.set_facecolor("black")

    A_cells = cell_embedding_matrix[cell_embedding_matrix["FINAL_CELL_TYPE"] == A]
    B_cells = cell_embedding_matrix[cell_embedding_matrix["FINAL_CELL_TYPE"] == B]

    if A_cells[B].isna().any():
        return
    
    cuts = pd.cut(A_cells[B], bins = bins)

    try:
        cuts = cuts.cat.rename_categories(lambda x: pd.Interval(abs(max(round(x.left, 2), 0)), 
                                                                round(x.right, 2), 
                                                                closed = x.closed))
    except:
        pass
    
    if A != B:
        ax.scatter(B_cells["X"], B_cells["Y"], s = 20, label = B, color = cell_type_colors[B]["color"])

    for category in sorted(cuts.unique()):
        A_category_identifiers = cuts[cuts == category].index
        A_category_cells = A_cells.loc[A_category_identifiers]

        bin_index = sorted(cuts.unique()).index(category)
        base_color = np.array(matplotlib.colors.to_rgb(cell_type_colors[A]["color"]))

        t = bin_index / max(len(sorted(cuts.unique())) - 1, 1)
        t = 0.9 * t + 0.1

        shade = tuple((1 - t) * np.ones(3) + t * base_color)

        ax.scatter(A_category_cells["X"], A_category_cells["Y"], s = 20, label = f"{A} {category}", color = shade)

    handles, labels = ax.get_legend_handles_labels()
    legend = legend_ax.legend(
        handles, labels,
        loc = "center",
        frameon = True
    )

    legend.get_frame().set_facecolor("black")

    for text in legend.get_texts():
        text.set_color("white")

    save_path = os.path.join(save_path, B)
    os.makedirs(save_path, exist_ok = True)

    fig.savefig(os.path.join(save_path, f"{B}_{A}.png"), dpi = 300, bbox_inches = "tight", pad_inches = 0.2)

    plt.close(fig)

def driver(assignments_directory, signature_matrix, save_path, filter, log_transform, bandwidth, cell_type_colors):
    """driver function that generates a cell embedding matrix (and associated plots) for each image in a user-specified directory

    args:
        assignments_directory (str): directory path that houses cell type assignment spreadsheets from CELESTA
        signature_matrix (str): path to a signature matrix file that contains the cell types of interest
        save_path (str): parent directory path where output (LCLQ embeddings, cluster visualizations, etc.) organized by sample subdirectories will be stored
        filter (str): comma-separated list of sample names to process, or 'all' to process everything
        log_transform (bool): whether to log-transform the colocalization statistics
        bandwidth (int): bandwidth parameter used to control the radius of the spatial neighborhood around a point
        cell_type_colors (dictionary): dictionary mapping cell types to their colors, used for plotting
    """

    colocalizations_directory = os.path.join(save_path, "colocalizations", "local_colocalizations")

    sample_assignments = glob.glob(os.path.join(assignments_directory, "*_assignments.csv"))
    if filter != "all":
        sample_assignments = [f for f in sample_assignments if os.path.basename(f) in 
                                 [f"{name}_assignments.csv" for name in filter.split(",")]]
        
    signature_matrix = pd.read_csv(signature_matrix)

    for sample_path in sample_assignments:
        sample_name = re.match(r"(.+?)_assignments", os.path.basename(sample_path)).group(1)

        print("computing LCLQ embeddings for sample", sample_name)
        
        assignments = pd.read_csv(sample_path).dropna().reset_index(drop = True)

        sample_visualizations_directory = os.path.join(save_path, "visualizations", "cell_colocalizations", sample_name)
        os.makedirs(sample_visualizations_directory, exist_ok = True)

        cell_embedding_matrix = generate_spatial_cell_embedding_matrix(assignments, colocalizations_directory, sample_name, signature_matrix["CELL_TYPE"].tolist(), log_transform, bandwidth)

        unique_cell_types = cell_embedding_matrix["FINAL_CELL_TYPE"].unique()
        for B in unique_cell_types:
            for A in unique_cell_types:
                plot_image_pairwise_colocalizations(cell_embedding_matrix, A, B, sample_visualizations_directory, cell_type_colors)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for generating LCLQ embeddings")

    parser.add_argument("-b", "--bandwidth", type = int, default = 100, help = "bandwidth radius of circular neighborhood to consider; \
                                                                                      defaults to 100")
    parser.add_argument("-c", "--colormap_path", help = "file path that points to the underlying location of the colormap.json file")
    parser.add_argument("-d", "--data_directory", help = "path to a data directory of cell assignment .csv files")
    parser.add_argument("-f", "--filter", default = "all", help = "comma-separated list of sample names to process, or 'all' to process everything; \
                                                                   defaults to 'all'")
    parser.add_argument("-m", "--signature_matrix", help = "path to a signature matrix file")
    parser.add_argument("-s", "--save_path", help = "path to save LCLQ embeddings and visualizations under colocalizations/local_colocalizations and visualizations/cell_colocalizations, respectively")
    
    parser.add_argument("--log_transform", action = "store_true", help = "whether to log-transform the colocalization statistics")

    return parser.parse_args()
    
def main():
    arguments = parse_arguments()

    bandwidth = arguments.bandwidth
    colormap_path = arguments.colormap_path
    data_directory = arguments.data_directory
    filter = arguments.filter
    signature_matrix = arguments.signature_matrix
    save_path = arguments.save_path
    log_transform = arguments.log_transform

    with open(colormap_path, "r") as file:
        cell_type_colors = json.load(file)

    os.makedirs(save_path, exist_ok = True)
    os.makedirs(os.path.join(save_path, "colocalizations", "local_colocalizations"), exist_ok = True)
    os.makedirs(os.path.join(save_path, "visualizations", "cell_colocalizations"), exist_ok = True)
    
    driver(data_directory, signature_matrix, save_path, filter, log_transform, bandwidth, cell_type_colors)

main()