from magicgui import magicgui
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas # type: ignore
from PyQt5.QtWidgets import QVBoxLayout, QWidget

class CellProportionWidget(QWidget):
    """widget to display the proportion of each cell type in the dataset"""
    
    def __init__(self, assignments):
        """initialize the cell proportion widget

        args:
            assignments (dataframe): dataframe containing cell type assignments with a "FINAL_CELL_TYPE" column
        """

        super().__init__()

        self.assignments = assignments
        self.initialize_interface()

    def initialize_interface(self):
        """initialize the user interface components for the cell proportion display widget"""
        
        layout = QVBoxLayout()

        self.canvas = FigureCanvas(plt.figure(figsize = (6, 8)))
        
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.update_plot()

    def update_plot(self):
        """creates a bar chart showing the percentage of each cell type in the dataset"""

        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111)
        
        cell_counts = self.assignments["FINAL_CELL_TYPE"].value_counts()
        cell_proportions = (cell_counts / cell_counts.sum()) * 100
        
        ax.bar(cell_proportions.index, cell_proportions.values, color = sns.color_palette("husl", len(cell_proportions)))
        ax.set_xlabel("cell type")
        ax.set_ylabel("percentage (%)")
        ax.tick_params(axis = "x", rotation = 90)

        plt.tight_layout()
        
        self.canvas.draw()
    
@magicgui()
def show_cell_type_proportions(viewer, assignments):
    """display a widget showing the proportions of each cell type in the dataset

    args:
        viewer: the napari viewer instance to add the widget to
        assignments: dataframe containing cell type assignments with a "FINAL_CELL_TYPE" column
    """

    viewer.window.add_dock_widget(CellProportionWidget(assignments), name = "cell type proportions")