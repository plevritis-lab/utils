# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EPOCHS is a modular CODEX spatial proteomics pipeline. It converts raw multiplexed images into:
- segmentation masks and per-cell measurements
- CELESTA cell type assignments
- spatial colocalization metrics
- topology features (persistent homology)

## Development Setup

Install dependencies:

```bash
uv sync
```

For HPC environments:

```bash
uv sync --extra hpc
```

## Pipeline Summary

```
Raw QPTIFF/TIF images
    -> 1_segmentation   (cell masks + per-cell measurements)
    -> 2_celesta        (cell type assignments + QC plots)
    -> 3_colocalization (global + local spatial co-occurrence)
    -> 4_topology       (persistent homology features)
```

Modules can run independently if expected inputs exist.

## Critical Inputs Across Modules

- images: `.qptiff` / `.tif` / `.tiff`
- `channel_names.txt` (needed when metadata is missing / plain TIF workflows)
- `signature_matrix.csv` with `CELL_TYPE`
- `colormap.json` with matching cell type labels and `"unknown"`

## Module Contracts

- `1_segmentation` -> outputs `*_cell_measurements.csv`
- `2_celesta` -> outputs `*_assignments.csv` with `FINAL_CELL_TYPE`
- `3_colocalization` -> outputs local / global colocalization matrices
- `4_topology` -> outputs `*_persistence.json` + persistence diagrams

## Running Scripts

Segmentation pipeline (config-driven, runs reformat + segment + quantify):

```bash
uv run python3 1_segmentation/src/run_pipeline.py --config path/to/config.yaml
uv run python3 1_segmentation/src/run_pipeline.py --config path/to/config.yaml --force
```

See `1_segmentation/conf/config_template.yaml` for the config format.

CELESTA pipeline (config-driven, runs thresholds + assignment + visualization):

```bash
uv run python3 2_celesta/src/run_pipeline.py --config path/to/config.yaml
uv run python3 2_celesta/src/run_pipeline.py --config path/to/config.yaml --filter sample1,sample2
```

See `2_celesta/conf/config_template.yaml` for the config format.

Other Python scripts:

```bash
uv run python3 3_colocalization/src/local_colocalization/generate_spatial_embeddings.py --help
uv run python3 4_topology/src/apply_persistent_homology.py --help
```

R scripts:

```bash
Rscript 3_colocalization/src/global_colocalization/generate_condition_summaries.R --help
```

### Module Layout

- **1_segmentation**: cell detection + per-cell quantification
- **2_celesta**: CELESTA assignment + static/interactive visualization
- **3_colocalization**: local LCLQ (Python) + global CLQ (R/metadisco)
- **4_topology**: persistent homology with GUDHI

### Key Conventions

- Image arrays use shape `(channels, height, width)`
- Cell data CSVs use columns: `CELL_IDENTIFIER`, `X`, `Y`, `FINAL_CELL_TYPE`
- Marker intensity columns are uppercase (for example `E_CADHERIN`)
- Colormaps map cell type names to hex colors
- Python CLIs use `argparse` and support `--help`
