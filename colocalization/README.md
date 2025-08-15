# colocalization

this module provides tools for analyzing spatial colocalization at both global and local levels in spatial proteomics data.

## directory structure

- `source/`: main analysis scripts and submodules.
  - `global_colocalization/`: r and python scripts for global colocalization analysis, including summary generation, visualization, and density correlation analysis.
  - `local_colocalization/`: python scripts for identifying spatial clusters and generating spatial embeddings.
- `bash/`: shell scripts for batch processing.
  - `local_colocalization_template.sh`, `global_colocalization_template.sh`: templates for batch colocalization analysis.
  - `plevritis/`: project-specific batch scripts for different sample types (e.g., `local_colocalization_edges_negative_plevritis.sh`).

## usage notes
- the module supports both r and python workflows for comprehensive colocalization analysis.
- batch scripts are organized by project/sample type for reproducibility.

## directory structure changes
- the module is organized to separate source code (`source/`), batch scripts (`bash/`), and project-specific workflows (`plevritis/`).