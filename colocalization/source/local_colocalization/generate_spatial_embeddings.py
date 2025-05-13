import glob
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
from scipy.spatial.distance import cdist

def __compute_local_colocation_quotient(assignments, A, B, bandwidth):
    """calculates a local indicator of the colocation quotient (LCLQ) between two cell types in a single image, which can serve as a per-cell feature embedding

    args:
        assignments (dataframe): generated cell type assignments from CELESTA
        A (str): first cell type of interest (e.g., epithelial)
        B (str): second cell type of interest (e.g., granulocyte)
        bandwidth (int): bandwidth parameter used to assign geospatial weights in the gaussian kernel density function; \
                         a smaller bandwidth typically produces smaller and spiky clusters whereas a larger bandwidth results in smoother and larger clusters
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

    colocalization_embedding = N_AB / P if P != 0 else [np.nan]
    return colocalization_embedding

def generate_spatial_cell_embedding_matrix(assignments, save_path, bandwidth = 200):
    """constructs and saves a spatial embedding matrix for a single image, representing each cell as a function of their LCLQ value with respect to all other cell types

    args:
        assignments (dataframe): generated cell type assignments from CELESTA
        save_path (str): sample subdirectory where LCLQ embeddings will be stored
        bandwidth (int, optional): bandwidth parameter used to assign geospatial weights in the gaussian kernel density function; \
                                   a smaller bandwidth typically produces smaller and spiky clusters whereas a larger bandwidth results in smoother and larger clusters; \
                                   defaults to 200
    """

    unique_cell_types = assignments["FINAL_CELL_TYPE"].unique()

    composite_embedding_matrix = [None] * unique_cell_types.size

    for i, A in enumerate(unique_cell_types):
        A_cells = assignments[assignments["FINAL_CELL_TYPE"] == A]

        A_embeddings = {B : __compute_local_colocation_quotient(assignments, A, B, bandwidth) for B in unique_cell_types}
        A_embedding_matrix = pd.DataFrame({
            "FINAL_CELL_TYPE" : A,
            "X" : A_cells["X"].values, 
            "Y" : A_cells["Y"].values,
            **A_embeddings
        })

        composite_embedding_matrix[i] = A_embedding_matrix
            
    composite_embedding_matrix = pd.concat(composite_embedding_matrix, ignore_index = True)
    composite_embedding_matrix.to_csv(os.path.join(save_path, "spatial_embeddings.csv"), index = False, na_rep = "nan")

    return composite_embedding_matrix

def visualize_bandwidth_parameter(centroids, bandwidth, save_path):
    """plots and saves a visual representation of the bandwidth parameter (using Euclidean distances) on a randomly chosen point from user-provided centroids

    args:
        centroids (array): array of centroids of shape N x 2 where N is the total number of data points
        bandwidth (int): bandwidth parameter used to control the radius of the spatial neighborhood around a point
        save_path (str): sample directory where a visual representation of the bandwidth parameter will be stored
    """

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

def visualize_bandwidth_versus_cells(centroids, bandwidth_range, save_path):
    """plots and saves a line chart showing the average number of cells per image across a range of bandwidths

    args:
        centroids (array): array of centroids of shape N x 2 where N is the total number of data points
        bandwidth_range (array): range of bandwidth parameters to analyze
        save_path (str): sample directory where results from the bandwidth analysis will be stored
    """

    distances = cdist(centroids, centroids) # N x N
    average_cells = [0] * bandwidth_range.size

    # TODO - free up memory consumed by bandwidth_distances, write this function?

    for i, bandwidth in enumerate(bandwidth_range):
        bandwidth_distances = np.where(distances > bandwidth, 0, 1)
        average_cells[i] = bandwidth_distances.sum(axis = 1).mean()
    
    print(average_cells)

    # fig, ax = plt.subplots(figsize = (10, 8))

    # selected_cell = centroids[np.random.choice(centroids.shape[0])]
    # ax.scatter(centroids[:, 0], centroids[:, 1], color = "white")
    # ax.scatter(selected_cell[0], selected_cell[1], color = "red")    

    # spatial_neighborhood = plt.Circle((selected_cell[0], selected_cell[1]), bandwidth, color = "red", 
    #                                        linestyle = "--", linewidth = 2, fill = False)
    # plt.gca().add_patch(spatial_neighborhood)

    # ax.set_facecolor("black")

    # fig.savefig(os.path.join(save_path, f"{bandwidth}_bandwidth_visualization.pdf"))
    # plt.close(fig)

def plot_image_pairwise_colocalizations(cell_embedding_matrix, A, B, save_path, bins = 4):
    """plots and saves the euclidean-based LCLQ measurements between two cell types for a single image, with A's measurements with respect to B graded into <= the number of user-specified number of categories

    args:
        cell_embedding_matrix (dataframe): generated LCLQ measurements for every cell in the image
        A (str): first cell type of interest (e.g., epithelial)
        B (str): second cell type of interest (e.g., granulocyte)
        save_path (str): sample directory where pairwise LCLQ plots will be stored under their own subdirectory (pairwise_plots)
        bins (int, optional): used to determine the maximum number of LCLQ categories to display; \
                              defaults to 4
    """

    # TODO - clean up plots (fix formatting of numbers, consistent sizes across images despite legends, etc.)

    matplotlib.use("agg")

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
    
    fig, ax = plt.subplots(figsize = (10, 8))

    if A != B:
        if B == "PDPN+ cell":
            C = "Mesenchymal cell (podoplanin+)"
        else:
            C = B

        ax.scatter(B_cells["X"], B_cells["Y"], s = 20, label = C)

    for category in sorted(cuts.unique()):
        A_category_identifiers = cuts[cuts == category].index
        A_category_cells = A_cells.loc[A_category_identifiers]

        if A == "Epithelial cell (Cytokeratin+)":
            D = "Epithelial cell (cytokeratin+)"
        else:
            D = A

        ax.scatter(A_category_cells["X"], A_category_cells["Y"], s = 20, label = f"{D} {category}")

    ax.set_facecolor("black")

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.legend(loc = "center", bbox_to_anchor = (0.5, -0.1))

    fig.tight_layout()

    save_path = os.path.join(save_path, "pairwise_plots")
    os.makedirs(save_path, exist_ok = True)

    fig.savefig(os.path.join(save_path, f"{A}_{B}.pdf"))

    plt.close(fig)

def driver(assignments_directory, save_path):
    """driver function that generates a cell embedding matrix (and associated plots) for each image in a user-specified directory

    args:
        assignments_directory (str): directory path that houses cell type assignment spreadsheets from CELESTA
        save_path (str): parent directory path where output (LCLQ embeddings, cluster visualizations, etc.) organized by sample subdirectories will be stored
    """

    for sample_path in glob.glob(os.path.join(assignments_directory, "*")):
        sample_name = re.match(r"(.+?)_assignments", os.path.basename(sample_path)).group(1)

        print("computing LCLQ embeddings for sample", sample_name)
        
        assignments = pd.read_csv(sample_path).dropna().reset_index(drop = True)
        assignments.columns = assignments.columns.str.upper()
        assignments.columns = assignments.columns.str.replace(" ", "_")

        sample_directory = os.path.join(save_path, sample_name)
        os.makedirs(sample_directory, exist_ok = True)

        cell_embedding_matrix = generate_spatial_cell_embedding_matrix(assignments, sample_directory)

        unique_cell_types = cell_embedding_matrix["FINAL_CELL_TYPE"].unique()
        for A in unique_cell_types:
            for B in unique_cell_types:
                plot_image_pairwise_colocalizations(cell_embedding_matrix, A, B, sample_directory)

if __name__ == "__main__":
    # conditions = ["whole_slide_images/involved_lymph_node",
    #               "whole_slide_images/primary_node_negative",
    #               "whole_slide_images/primary_node_positive"]
    
    # num_processes = len(conditions)

    # processes = [None] * num_processes
    # for process_identifier in range(num_processes):
    #     process = multiprocessing.Process(target = driver, args = (os.path.join("../input", conditions[process_identifier]),
    #                                                                os.path.join("../output", conditions[process_identifier])))
        
    #     processes[process_identifier] = process
    #     processes[process_identifier].start()

    # for process_identifier in range(num_processes):
    #     processes[process_identifier].join()

    driver("/Users/rohit/Downloads/example", "/Users/rohit/Downloads/example_output")

# TODO - convert the LCLQ calculation in terms of MICRONS instead of PIXELS; convert pixels to microns for more accurate results