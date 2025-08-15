# celesta

the celesta module will visualize overlays of assignments, marker probabilities, histology, and proteomic images **in real-time**, allowing the user to toggle back-and-forth between these different modes

source files:

1. `source/apply_celesta.R`: batch applies celesta to a data directory of CSV files
2. `source/dynamic_overlays.py`: visualizes composite overlays dynamically using `napari`

---

## detailed module overview and directory structure

### directory structure

- `source/`: contains all main processing and visualization scripts for celesta.
  - `apply_celesta.r`: batch applies celesta to a directory of csv files, assigning cell types using a signature matrix and thresholds.
  - `visualize_assignments.py`: generates static visualizations of cell type assignments and cell proportions from assignment csvs and colormap jsons.
  - `visualize_dynamic_overlays.py`: provides an interactive napari-based viewer for overlays of assignments, marker probabilities, histology, and proteomic images, with dynamic toggling and stratification.
  - `generate_thresholds.py`: generates or updates threshold spreadsheets for spatial proteomics data, based on a signature matrix and image directory structure.
  - `generate_quality_control.py`: produces a quality control spreadsheet for spatial proteomics data, combining clinical data and protein panel information.
  - `colormaps/`: contains json colormap templates and project-specific colormaps for visualizations.
  - `widgets/`: contains reusable python widgets for visualization, such as `cell_proportions.py` for interactive cell proportion displays.
- `bash/`: contains shell scripts for batch processing.
  - `celesta_processing_template.sh`: template for batch celesta processing.
- `requirements.txt`: lists python dependencies for running celesta scripts.

### usage notes
- the celesta module is designed for both batch and interactive processing of spatial proteomics data, supporting both R and python workflows.
- visualization scripts support both static (matplotlib) and interactive (napari) outputs, with customizable colormaps and overlays.
- batch scripts in `bash/` are organized by project and sample type for reproducible, large-scale processing.