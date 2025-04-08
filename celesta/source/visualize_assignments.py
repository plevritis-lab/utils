import argparse
import json
import os
import pandas as pd
import matplotlib.pyplot as plt

BASE_SIZE = 20
POINT_SIZE = 8

def visualize_assignments(assignments, cell_type_info, display_cells, save_path):
    """visualizes cell type assignments statically

    args:
        assignments (dataframe): dataframe containing cell positions ["X", "Y"] and cell type assignments ["FINAL_CELL_TYPE"]
        cell_type_info (dictionary): dictionary mapping cell types to their display properties, including colors
        display_cells (dictionary): list of cell type names to display color in the visualization
        save_path (str): directory path where the visualization images will be saved
    """

    x_range = assignments["X"].max() - assignments["X"].min()
    y_range = assignments["Y"].max() - assignments["Y"].min()

    aspect_ratio = x_range / y_range
    
    if aspect_ratio > 1:
        width = BASE_SIZE
        height = BASE_SIZE / aspect_ratio
    else:
        width = BASE_SIZE * aspect_ratio
        height = BASE_SIZE

    fig, ax = plt.subplots(figsize = (width, height), facecolor = "black")
    
    ax.invert_yaxis()
    ax.set_axis_off()
    
    for cell_type in display_cells:
        subset = assignments[assignments["FINAL_CELL_TYPE"] == cell_type]
        
        x = subset["X"]
        y = subset["Y"]

        color = cell_type_info[cell_type]["color"]
        
        ax.scatter(x, y, c = color, label = cell_type, s = POINT_SIZE)
    
    legend_fig, legend_ax = plt.subplots(figsize = (5, len(display_cells) * 0.4), facecolor = "black")

    legend_ax.set_axis_off()
    
    handles, labels = ax.get_legend_handles_labels()
    legend_ax.legend(
        handles = handles, 
        labels = labels,
        loc = "center", 
        facecolor = "black", 
        edgecolor = "white", 
        labelcolor = "white",
        markerscale = 3,
        prop = {"size": 12, 
                "weight": "bold"}
    )

    legend_fig.tight_layout()
    fig.tight_layout()
    
    legend_fig.savefig(os.path.join(save_path, "assignments_legend.png"), dpi = 300)
    fig.savefig(os.path.join(save_path, "assignments.png"), dpi = 300)

def visualize_cell_proportions(assignments, cell_type_info, save_path):
    """visualizes cell proportions as a stacked bar chart
    
    args:
        assignments (dataframe): dataframe containing cell type assignments ["FINAL_CELL_TYPE"]
        cell_type_info (dictionary): dictionary mapping cell types to their display properties, including colors
        save_path (str): directory path where the visualization images will be saved
    """
    
    fig, ax = plt.subplots(figsize = (4, 8))

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 1)

    ax.set_xticks([])
    ax.set_yticks([0.25, 0.5, 0.75])
    ax.set_yticklabels(["25%", "50%", "75%"], fontsize = 8, 
                         weight = "bold", va = "center")
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis = "y", length = 0)

    proportions = assignments["FINAL_CELL_TYPE"].value_counts() / assignments.shape[0]
    proportions = [proportions.get(cell_type, 0) for cell_type in cell_type_info.keys()]
    
    colors = [cell_type_info[cell_type]["color"] for cell_type in cell_type_info.keys()]

    ax.bar([0], [1], color = "white", edgecolor = "black", linewidth = 1, width = 0.5)
    
    bottom = 0
    for proportion, color, cell_type in zip(proportions, colors, cell_type_info.keys()):
        ax.bar([0], [proportion], color = color, edgecolor = "black", linewidth = 1, 
                   width = 0.5, bottom = bottom, label = cell_type)
        bottom += proportion
    
    legend_fig, legend_ax = plt.subplots(figsize = (5, len(cell_type_info) * 0.4))

    legend_ax.set_axis_off()

    handles, labels = ax.get_legend_handles_labels()
    legend_ax.legend(
        handles = handles,
        labels = labels,
        loc = "center",
        facecolor = "white",
        edgecolor = "black",
        labelcolor = "black",
        markerscale = 1.5,
        prop = {"size": 10, 
                "weight": "bold"}
    )

    legend_fig.tight_layout()
    fig.tight_layout()

    legend_fig.savefig(os.path.join(save_path, "proportions_legend.png"), dpi = 300)    
    fig.savefig(os.path.join(save_path, "proportions.png"), dpi = 300)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for static celesta assignment visualization")

    parser.add_argument("-a", "--assignments_path", help = "file path that points to the underlying location of the assignments.csv file")
    parser.add_argument("-c", "--colormap_path", help = "file path that points to the underlying location of the colormap.json file")
    parser.add_argument("-d", "--display_cells", default = "all", help = "list of cell types to display in the visualization; \
                                                                          this argument should be a comma separated list with no spaces or 'all'; \
                                                                          defaults to 'all'")
    parser.add_argument("-s", "--save_path", help = "directory path that points to the underlying location where output will be written")

    return parser.parse_args()
    
def main():
    arguments = parse_arguments()

    assignments_path = arguments.assignments_path
    colormap_path = arguments.colormap_path
    display_cells = arguments.display_cells
    save_path = arguments.save_path

    with open(colormap_path, "r") as file:
        cell_type_info = json.load(file)

    display_cells = display_cells.split(",") if display_cells != "all" else list(cell_type_info.keys())
    assignments = pd.read_csv(assignments_path).dropna()

    os.makedirs(save_path, exist_ok = True)
    
    visualize_assignments(assignments, cell_type_info, display_cells, save_path)
    visualize_cell_proportions(assignments, cell_type_info, save_path)

main()