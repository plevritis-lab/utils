import numpy as np
from tifffile import imread
from matplotlib.patches import Circle
import matplotlib.pyplot as plt

def visualize_movable_circle(proteomic_path, histology_path, panel_path, markers, radius = 300):
    """visualizes a movable circle in histology and proteomic images with keyboard control; \
       press arrow keys to move circle, enter to select region of interest
    
    args:
        proteomic_path (str): file path to the proteomic image of shape (c, y, x)
        histology_path (str): file path to the histology image of shape (c, y, x)
        panel_path (str): file path to the protein panel
        markers (dictionary): dictionary mapping protein markers to colors ('red', 'green', or 'yellow')
        radius (int): radius of the circle in pixels; \
                      defaults to 300
    """

    proteomic_image = imread(proteomic_path)
    histology_image = imread(histology_path)
    
    with open(panel_path, "r") as panel:
        protein_panel = [m.rstrip() for m in panel]
    
    proteomic_image = (((proteomic_image - proteomic_image.min()) / (proteomic_image.max() - proteomic_image.min())) * 255).astype(np.uint8)
    histology_image = (((histology_image - histology_image.min()) / (histology_image.max() - histology_image.min())) * 255).astype(np.uint8)
    
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
            center[0] -= 5
        elif event.key == "right":
            center[0] += 5
        elif event.key == "up":
            center[1] -= 5
        elif event.key == "down":
            center[1] += 5
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
    
    masked_histology = (histology_image * mask_histology) / 255.0
    masked_proteomics = (proteomic_image * mask_proteomics) / 255.0
    
    y_index, x_index = np.where(mask)
    y_min, y_max = y_index.min(), y_index.max() + 1
    x_min, x_max = x_index.min(), x_index.max() + 1
    
    cropped_histology = masked_histology[:, y_min : y_max, x_min : x_max]
    cropped_proteomics = masked_proteomics[:, y_min : y_max, x_min : x_max]
    
    marker_indices = [protein_panel.index(marker) for marker in markers.keys()]
    
    channels = []
    for index, color in zip(marker_indices, markers.values()):
        channel = cropped_proteomics[index]
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
    
    cropped_proteomics = np.max(channels, axis = 0)
    
    spacing = np.zeros((3, 20, cropped_histology.shape[2]), dtype = cropped_histology.dtype)
    stacked_image = np.concatenate([cropped_histology, spacing, cropped_proteomics], axis = 1)
    
    fig = plt.figure(figsize = (12, 8))
    
    plt.axis("off")
    plt.imshow(stacked_image.transpose(1, 2, 0))
    plt.savefig("masked_regions.pdf", bbox_inches = "tight")
    plt.close(fig)