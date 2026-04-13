# EPOCHS

EPOCHS is a modular pipeline for CODEX spatial proteomics. It turns raw multiplexed images into:
- cell masks and per-cell measurements
- cell type assignments
- spatial colocalization statistics
- topology features

## Quick Start

Install dependencies:

```bash
uv sync
```

For HPC environments:

```bash
uv sync --extra hpc
```

## Pipeline At A Glance

```
Raw QPTIFF/TIF images
    -> 1_segmentation   (cell masks + per-cell measurements)
    -> 2_celesta        (cell type assignments + QC plots)
    -> 3_colocalization (global + local spatial co-occurrence)
    -> 4_topology       (persistent homology features)
```

Modules are usually run in this order, but each can run independently if its required inputs already exist.

## Required Inputs

| File | Format | Required For | Notes |
|---|---|---|---|
| Proteomic images | `.qptiff`, `.tif`, `.tiff` | Segmentation | One multi-channel image per sample |
| Channel names | `.txt` | Segmentation (sometimes) | One marker per line, in channel order; needed for plain TIF or missing QPTIFF metadata |
| Signature matrix | `.csv` | CELESTA | Must include `CELL_TYPE` |
| Colormap | `.json` | CELESTA + plots | Keys should match CELESTA output labels; include `"unknown"` |

## How To Run

Start with each module README:
1. [**1_segmentation**](1_segmentation/README.md)
2. [**2_celesta**](2_celesta/README.md)
3. [**3_colocalization**](3_colocalization/README.md)
4. [**4_topology**](4_topology/README.md)

## Core Outputs

After a full run, look for these key outputs:

- `quantifications/` -> per-cell marker measurements (`*_cell_measurements.csv`)
- `assignments/` -> CELESTA cell type assignments (`*_assignments.csv`)
- `colocalizations/` -> local + global co-occurrence matrices
- `topology/` -> persistence features (`*_persistence.json`)
- `visualizations/` -> assignment plots, legends, colocalization maps, persistence diagrams

## Minimal Directory Map

```
your_project/
├── data/{condition}/
│   ├── {sample}/data/{sample}.qptiff
│   ├── quantifications/{method}/*_cell_measurements.csv
│   ├── assignments/*_assignments.csv
│   ├── colocalizations/{local_colocalizations,global_colocalizations}/
│   ├── topology/*_persistence.json
│   └── visualizations/
└── celesta/
    ├── signature_matrices/{condition}/
    └── thresholds/{condition}/
```

## Conventions

- Image array shape: `(channels, height, width)`
- Expected assignment columns: `CELL_IDENTIFIER`, `X`, `Y`, `FINAL_CELL_TYPE`
- Marker columns are uppercase in CSVs (example: `E_CADHERIN`)
- Python scripts use `argparse` and support `--help`
