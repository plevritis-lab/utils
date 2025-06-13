import argparse
import glob
import gudhi
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import re
from sklearn.preprocessing import MinMaxScaler

def plot_persistence_diagram_density(persistence, A, save_path):
    """plots a persistence diagram and its density for a given cell type A
    
    args:
        persistence (list): list of persistence pairs, where each pair is a dictionary with keys dimension, birth, and death
        A (str): cell type of interest (e.g., epithelial)
        save_path (str): sample directory where persistence visualizations will be stored
    """

    matplotlib.use("agg")

    persistence = [(item["dimension"], (item["birth"], item["death"])) for item in persistence]

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5))

    gudhi.plot_persistence_diagram(persistence, axes = axes[0], max_intervals = 0)
    gudhi.plot_persistence_density(persistence, axes = axes[1], max_intervals = 0)

    fig.savefig(os.path.join(save_path, f"{A}_persistence_diagram_density.png"), dpi = 300)
    plt.close(fig)

def featurize_cell_type_topology(cells):
    """featurizes the topology of a given cell type using persistent homology

    args:
        cells (dataframe): dataframe containing cell coordinates X and Y
    """

    centroids = MinMaxScaler().fit_transform(cells[["X", "Y"]])

    alpha_complex = gudhi.AlphaComplex(points = centroids)
    simplex_tree = alpha_complex.create_simplex_tree()
    persistence = simplex_tree.persistence()

    return persistence

def compute_topological_features(assignments, save_path, sample_name, cell_types, count_criteria):
    """computes persistent homology features for a given sample across cell types individually
    
    args:
        assignments (dataframe): generated cell type assignments from CELESTA
        save_path (str): directory where persistent homology features will be saved
        sample_name (str): name of the sample being processed
        cell_types (list): list of cell types to compute persistent homology features for; \
                           if entries are not present in the assignments dataframe, its corresponding values will be set to NA
        count_criteria (int): minimum number of points required to compute persistent homology for a cell type
    """

    unique_cell_types = ["unknown"] + cell_types

    composite_persistence_matrix = {}

    for A in unique_cell_types:
        A_cells = assignments[assignments["FINAL_CELL_TYPE"] == A]

        if len(A_cells) < count_criteria:
            print("skipping cell type", A, "in sample", sample_name, "due to insufficient points")

            persistence_data = [
                {
                    "dimension": None,
                    "birth": None,
                    "death": None
                }
            ]

        else:
            persistence = featurize_cell_type_topology(A_cells)

            persistence_data = [
                {
                    "dimension": dimension,
                    "birth": birth,
                    "death": death
                }

                for dimension, (birth, death) in persistence
            ]

        composite_persistence_matrix[A] = persistence_data

    output_path = os.path.join(save_path, f"{sample_name}_persistence.json")
    pd.DataFrame([{"cell_type": A, "persistence": composite_persistence_matrix[A]} 
                     for A in composite_persistence_matrix]).to_json(output_path, orient = "records", indent = 2)
    
    return composite_persistence_matrix

def driver(assignments_directory, signature_matrix, save_path, filter, count_criteria):
    """driver function that computes persistent homology features for all samples in a given data directory

    args:
        assignments_directory (str): path to a data directory of cell assignment .csv files
        signature_matrix (str): path to a signature matrix file that contains the cell types of interest
        save_path (str): path to save persistent homology features and diagrams under topology and visualizations/persistence_diagrams, respectively
        filter (str): comma-separated list of sample names to process, or 'all' to process everything
        count_criteria (int): minimum number of points required to compute persistent homology for a cell type
    """

    topology_directory = os.path.join(save_path, "topology")

    sample_assignments = glob.glob(os.path.join(assignments_directory, "*_assignments.csv"))
    if filter != "all":
        sample_assignments = [f for f in sample_assignments if os.path.basename(f) in 
                                 [f"{name}_assignments.csv" for name in filter.split(",")]]
        
    signature_matrix = pd.read_csv(signature_matrix)
        
    for sample_path in sample_assignments:
        sample_name = re.match(r"(.+?)_assignments", os.path.basename(sample_path)).group(1)

        print("computing persistent homology features for sample", sample_name)
        
        assignments = pd.read_csv(sample_path).dropna().reset_index(drop = True)

        sample_visualizations_directory = os.path.join(save_path, "visualizations", "persistence_diagrams", sample_name)
        os.makedirs(sample_visualizations_directory, exist_ok = True)

        persistent_homology_features = compute_topological_features(assignments, topology_directory, sample_name, signature_matrix["CELL_TYPE"].tolist(), count_criteria)
        
        unique_cell_types = assignments["FINAL_CELL_TYPE"].unique()
        for A in unique_cell_types:
            if len(assignments[assignments["FINAL_CELL_TYPE"] == A]) >= count_criteria:
                plot_persistence_diagram_density(persistent_homology_features[A], A, sample_visualizations_directory)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for applying persistent homology to 2D spatial point patterns")
    
    parser.add_argument("-c", "--count_critera", type = int, default = 3, help = "minimum number of points required to compute persistent homology for a cell type; \
                                                                                  defaults to 3")
    parser.add_argument("-d", "--data_directory", help = "path to a data directory of cell assignment .csv files")
    parser.add_argument("-f", "--filter", default = "all", help = "comma-separated list of sample names to process, or 'all' to process everything; \
                                                                   defaults to 'all'")
    parser.add_argument("-m", "--signature_matrix", help = "path to a signature matrix file")
    parser.add_argument("-s", "--save_path", help = "path to save persistent homology features and diagrams under topology and visualizations/persistence_diagrams, respectively")

    return parser.parse_args()

def main():
    arguments = parse_arguments()

    count_criteria = arguments.count_critera
    data_directory = arguments.data_directory
    filter = arguments.filter
    signature_matrix = arguments.signature_matrix
    save_path = arguments.save_path

    os.makedirs(save_path, exist_ok = True)
    os.makedirs(os.path.join(save_path, "topology"), exist_ok = True)
    os.makedirs(os.path.join(save_path, "visualizations", "persistence_diagrams"), exist_ok = True)

    driver(data_directory, signature_matrix, save_path, filter, count_criteria)

main()