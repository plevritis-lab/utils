from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from tifffile import imread
from skimage.transform import resize

SPEED = 15

def visualize_movable_circle(proteomic_path, histology_path, segmentation_path, assignments_path, colocalization_path, panel_path, markers, output_directory, sample_name, radius = 300):
    """visualizes a movable circle in histology and proteomic images with keyboard control; \
       press arrow keys to move circle, enter to select region of interest
    
    args:
        proteomic_path (str): file path to the proteomic image of shape (c, y, x)
        histology_path (str): file path to the histology image of shape (c, y, x)
        segmentation_path (str): file path to the segmentation image of shape (c, y, x)
        assignments_path (str): file path to the assignments file
        colocalization_path (str): file path to the colocalization file
        panel_path (str): file path to the protein panel
        markers (dictionary): dictionary mapping protein markers to colors ('red', 'green', or 'yellow')
        output_directory (str): directory to save the output image
        sample_name (str): name of the sample to be used in the output file name
        radius (int): radius of the circle in pixels; \
                      defaults to 300
    """

    # TODO - adjust scripts so that every output image is always written in (c, y, x) format and as a tiff

    proteomic_image = imread(proteomic_path)
    histology_image = imread(histology_path)
    segmentation_image = imread(segmentation_path).transpose(2, 0, 1)

    assignments = pd.read_csv(assignments_path)
    assignments = assignments[assignments["FINAL_CELL_TYPE"].isin(markers.keys())]
    colocalizations = pd.read_csv(colocalization_path)
    colocalizations = colocalizations[colocalizations["SAMPLE"] == sample_name]
    colocalizations = colocalizations[colocalizations["FINAL_CELL_TYPE"].isin(markers.keys())]
    
    with open(panel_path, "r") as panel:
        protein_panel = [m.rstrip() for m in panel]
    
    proteomic_image = (((proteomic_image - proteomic_image.min()) / (proteomic_image.max() - proteomic_image.min())) * 255).astype(np.uint8)
    histology_image = (((histology_image - histology_image.min()) / (histology_image.max() - histology_image.min())) * 255).astype(np.uint8)
    segmentation_image = (((segmentation_image - segmentation_image.min()) / (segmentation_image.max() - segmentation_image.min())) * 255).astype(np.uint8)
    
    height, width = histology_image.shape[1:3]
    center = [width // 2, height // 2]
    
    dpi = plt.rcParams["figure.dpi"]
    fig, ax = plt.subplots(figsize = (width / dpi, height / dpi), dpi = dpi)
    
    circle = Circle(center, radius, fill = False, color = "red")

    ax.add_patch(circle)
    ax.imshow(histology_image.transpose(1, 2, 0))
    
    def on_key(event):
        nonlocal center

        if event.key == "left":
            center[0] -= SPEED
        elif event.key == "right":
            center[0] += SPEED
        elif event.key == "up":
            center[1] -= SPEED
        elif event.key == "down":
            center[1] += SPEED
        elif event.key == "enter":
            plt.close()
            return
        
        circle.center = center
        fig.canvas.draw_idle()
    
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    
    y, x = np.ogrid[:height, :width]
    mask = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2) <= radius
    
    mask_histology = np.repeat(mask[np.newaxis, :, :], histology_image.shape[0], axis = 0)
    mask_proteomics = np.repeat(mask[np.newaxis, :, :], proteomic_image.shape[0], axis = 0)
    mask_segmentation = np.repeat(mask[np.newaxis, :, :], segmentation_image.shape[0], axis = 0)
    
    masked_histology = (histology_image * mask_histology) / 255.0
    masked_proteomics = (proteomic_image * mask_proteomics) / 255.0
    masked_segmentation = (segmentation_image * mask_segmentation) / 255.0
    
    y_index, x_index = np.where(mask)
    y_min, y_max = y_index.min(), y_index.max() + 1
    x_min, x_max = x_index.min(), x_index.max() + 1
    
    cropped_histology = masked_histology[:, y_min : y_max, x_min : x_max]
    cropped_proteomics = masked_proteomics[:, y_min : y_max, x_min : x_max]
    cropped_segmentation = masked_segmentation[:, y_min : y_max, x_min : x_max]

    assignments_in_circle = []
    for _, row in assignments.iterrows():
        y, x = int(row["Y"]), int(row["X"])
        if mask[y, x]:
            assignments_in_circle.append(row)

    assignments_in_circle = pd.DataFrame(assignments_in_circle)

    colocalizations_in_circle = []
    for _, row in colocalizations.iterrows():
        y, x = int(row["Y"]), int(row["X"])
        if mask[y, x]:
            colocalizations_in_circle.append(row)

    colocalizations_in_circle = pd.DataFrame(colocalizations_in_circle)

    channels = []
    for cell_type, marker_info in markers.items():
        color = marker_info["color"]
        channel = cropped_proteomics[protein_panel.index(marker_info["name"])]

        if color == "red":
            channels.append([channel, 
                             np.zeros_like(channel), 
                             np.zeros_like(channel)])
        elif color == "green":
            channels.append([np.zeros_like(channel), 
                             channel, 
                             np.zeros_like(channel)])
        elif color == "yellow":
            channels.append([channel / 2, 
                             channel / 2, 
                             np.zeros_like(channel)])
        elif color == "blue":
            channels.append([np.zeros_like(channel),
                             np.zeros_like(channel),
                             channel])
        elif color == "orange":
            channels.append([channel,
                             channel / 2,
                             np.zeros_like(channel)])
    
    cropped_proteomics = np.max(channels, axis = 0)

    dot_plot = np.zeros_like(cropped_histology)
    if not assignments_in_circle.empty:
        for _, row in assignments_in_circle.iterrows():
            y, x = int(row["Y"]) - y_min, int(row["X"]) - x_min
            cell_type = row["FINAL_CELL_TYPE"]
            if 0 <= y < dot_plot.shape[1] and 0 <= x < dot_plot.shape[2]:
                color = markers[cell_type]["color"]
                # Create a sharper, larger circular mask for the dot
                yy, xx = np.ogrid[:dot_plot.shape[1], :dot_plot.shape[2]]
                mask = (yy - y) ** 2 + (xx - x) ** 2 <= 4 ** 2  # radius 4 for larger, sharper dots
                if color == "red":
                    dot_plot[0][mask] = 1.0
                    dot_plot[1][mask] = 0.0
                    dot_plot[2][mask] = 0.0
                elif color == "green":
                    dot_plot[0][mask] = 0.0
                    dot_plot[1][mask] = 1.0
                    dot_plot[2][mask] = 0.0
                elif color == "yellow":
                    dot_plot[0][mask] = 1.0
                    dot_plot[1][mask] = 1.0
                    dot_plot[2][mask] = 0.0
                elif color == "blue":
                    dot_plot[0][mask] = 0.0
                    dot_plot[1][mask] = 0.0
                    dot_plot[2][mask] = 1.0
                elif color == "orange":
                    dot_plot[0][mask] = 1.0
                    dot_plot[1][mask] = 0.5
                    dot_plot[2][mask] = 0.0

    # dot_plot_colocalization = np.zeros_like(cropped_histology)
    # if not colocalizations_in_circle.empty:
    #     fig_coloc, ax_coloc = plt.subplots(figsize=(cropped_histology.shape[2] / dpi, cropped_histology.shape[1] / dpi), dpi=dpi)
    #     ax_coloc.imshow(np.zeros_like(cropped_histology).transpose(1, 2, 0))  # blank background

    #     for _, row in colocalizations_in_circle.iterrows():
    #         y, x = int(row["Y"]) - y_min, int(row["X"]) - x_min
    #         cell_type = row["FINAL_CELL_TYPE"]
    #         cluster_label = str(row["CLUSTER"])
    #         if 0 <= y < dot_plot_colocalization.shape[1] and 0 <= x < dot_plot_colocalization.shape[2]:
    #             color = markers[cell_type]["color"]
    #             yy, xx = np.ogrid[:dot_plot_colocalization.shape[1], :dot_plot_colocalization.shape[2]]
    #             mask = (yy - y) ** 2 + (xx - x) ** 2 <= 4 ** 2
    #             if color == "red":
    #                 dot_plot_colocalization[0][mask] = 1.0
    #                 dot_plot_colocalization[1][mask] = 0.0
    #                 dot_plot_colocalization[2][mask] = 0.0
    #                 text_color = "red"
    #             elif color == "green":
    #                 dot_plot_colocalization[0][mask] = 0.0
    #                 dot_plot_colocalization[1][mask] = 1.0
    #                 dot_plot_colocalization[2][mask] = 0.0
    #                 text_color = "green"
    #             elif color == "yellow":
    #                 dot_plot_colocalization[0][mask] = 1.0
    #                 dot_plot_colocalization[1][mask] = 1.0
    #                 dot_plot_colocalization[2][mask] = 0.0
    #                 text_color = "gold"
    #             elif color == "blue":
    #                 dot_plot_colocalization[0][mask] = 0.0
    #                 dot_plot_colocalization[1][mask] = 0.0
    #                 dot_plot_colocalization[2][mask] = 1.0
    #                 text_color = "blue"
    #             elif color == "orange":
    #                 dot_plot_colocalization[0][mask] = 1.0
    #                 dot_plot_colocalization[1][mask] = 0.5
    #                 dot_plot_colocalization[2][mask] = 0.0
    #                 text_color = "orange"
    #             # Draw text label at centroid
    #             ax_coloc.text(x, y, cluster_label, color=text_color, fontsize=8, ha='center', va='center',
    #                           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'))

    #     ax_coloc.axis("off")
    #     fig_coloc.subplots_adjust(left=0, right=1, top=1, bottom=0)
    #     fig_coloc.canvas.draw()
    #     # Convert figure to numpy array
    #     label_img = np.frombuffer(fig_coloc.canvas.tostring_rgb(), dtype=np.uint8)
    #     width, height = fig_coloc.canvas.get_width_height()
    #     label_img = label_img.reshape((height, width, 3))
    #     label_img = label_img.astype(np.float32) / 255.0
    #     # Resize to match dot_plot_colocalization shape if needed
    #     if label_img.shape[0] != dot_plot_colocalization.shape[1] or label_img.shape[1] != dot_plot_colocalization.shape[2]:
    #         label_img = resize(label_img, (dot_plot_colocalization.shape[1], dot_plot_colocalization.shape[2], 3), preserve_range=True)
    #     # Overlay text labels on dot plot
    #     dot_plot_colocalization = np.clip(dot_plot_colocalization + label_img.transpose(2, 0, 1), 0, 1)
    #     plt.close(fig_coloc)

    spacing = np.zeros((3, cropped_histology.shape[1], 20), dtype=cropped_histology.dtype)
    stacked_image = np.concatenate([
        cropped_histology, spacing, 
        cropped_proteomics, spacing, 
        cropped_segmentation, spacing, 
        dot_plot
    ], axis=2)

    sample_name = os.path.splitext(os.path.basename(proteomic_path))[0]

    fig = plt.figure(figsize=(16, 8))
    plt.axis("off")
    plt.imshow(stacked_image.transpose(1, 2, 0))
    plt.savefig(os.path.join(output_directory, f"{sample_name}_masked_regions.pdf"), bbox_inches="tight")
    plt.close(fig)

proteomic_path = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/521_S1_reg024/data/521_S1_reg024.tif"
histology_path = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/histology/521_S1_reg024.tif"
segmentation_path = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/521_S1_reg024/full/combined_algorithm_overlays/image_1_cellpose_PDGFRB + CD3 + CYTOKERATIN_mesmer_PDGFRB + CD3 + CYTOKERATIN_display_PDGFRB + CD3 + CYTOKERATIN.tif"
assignments_path = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/assignments/521_S1_reg024_assignments.csv"
colocalization_path = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/local_spatial_summary/microarray_clusters.csv"

panel_path = "/Users/rohit/Desktop/plevritis_analysis/data/channel_names.txt"

markers = {
    "T cell (CD3+)" : {
        "color": "blue",
        "name": "CD3"
    },

    "endothelial cell (CD31+)" : {
        "color": "red",
        "name": "CD31"
    },

    "epithelial cell (cytokeratin+)" : {
        "color": "yellow",
        "name": "CYTOKERATIN"
    },

    "cytotoxic T cell (CD8+)" : {
        "color": "orange",
        "name": "CD8"
    },
}

output_directory = "/Users/rohit/Downloads/figures"

visualize_movable_circle(proteomic_path, histology_path, segmentation_path, assignments_path, colocalization_path, panel_path, markers, output_directory, "521_S1_reg024")