import numpy as np
from tifffile import imread
from matplotlib.patches import Circle
import matplotlib.pyplot as plt

def visualize_circular_regions(proteomic_path, histology_path, panel_path, markers):
    proteomic_image = imread(proteomic_path)
    histology_image = imread(histology_path)

    with open(panel_path, "r") as panel:
        protein_panel = [m.rstrip() for m in panel]

    proteomic_image = (((proteomic_image - proteomic_image.min()) / (proteomic_image.max() - proteomic_image.min())) * 255).astype(np.uint8)
    histology_image = (((histology_image - histology_image.min()) / (histology_image.max() - histology_image.min())) * 255).astype(np.uint8)

    class CircleDrawer:
        def __init__(self, ax):
            self.ax = ax
            self.start_point = None
            self.circle = None
            self.center = None
            self.radius = None
            
            self.press = self.ax.figure.canvas.mpl_connect("button_press_event", self.on_press)
            self.move = self.ax.figure.canvas.mpl_connect("motion_notify_event", self.on_move)
            self.release = self.ax.figure.canvas.mpl_connect("button_release_event", self.on_release)

        def on_press(self, event):
            if event.inaxes != self.ax:
                return
            
            self.start_point = (event.xdata, event.ydata)
            self.center = self.start_point
            if self.circle is not None:
                self.circle.remove()

            self.circle = Circle(self.center, 0, fill = False, color = "red")
            self.ax.add_patch(self.circle)

        def on_move(self, event):
            if self.start_point is None or event.inaxes != self.ax:
                return
            
            self.radius = np.sqrt((event.xdata - self.center[0]) ** 2 + 
                                  (event.ydata - self.center[1]) ** 2)
            self.circle.set_radius(self.radius)
            
            self.ax.figure.canvas.draw_idle()

        def on_release(self, event):
            if self.start_point is None or event.inaxes != self.ax:
                return
            
            self.process_selection()
            self.start_point = None

        def process_selection(self):
            if self.center is None or self.radius is None:
                return

            c, y, x = histology_image.shape
            y, x = np.ogrid[:y, :x]
            distance_from_center = np.sqrt((x - self.center[0]) ** 2 + (y - self.center[1]) ** 2)
            mask = distance_from_center <= self.radius

            mask_histology = np.repeat(mask[np.newaxis, :, :], c, axis = 0)
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

    dpi = plt.rcParams["figure.dpi"]
    height, width = histology_image.shape[1:3]
    figure_size = (width / dpi, height / dpi)

    fig, ax = plt.subplots(figsize = figure_size, dpi = dpi)
    ax.imshow(histology_image.transpose(1, 2, 0))

    ax.set_xticks([])
    ax.set_yticks([])

    def on_circle_drawn(event):
        if event.button == 1:
            plt.close(fig)

    fig.canvas.mpl_connect("button_release_event", on_circle_drawn)

    CircleDrawer(ax)
    plt.show()