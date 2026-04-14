# Segmentation

This module detects cells from multiplexed CODEX images and quantifies per-cell marker expression.

- **Step 1:** input reformatting (organizes flat images into expected directory layout)
- **Step 2:** segmentation mask generation (cellpose)
- **Step 3:** per-cell measurement extraction into CSV files

## Prerequisites

- Standard install: `uv sync`
- HPC/Linux extras: `uv sync --extra hpc`
- Segmentation is compute-heavy; use HPC for full-resolution runs.

## Required Inputs

| File | Format | Notes |
|---|---|---|
| Proteomic image | `.qptiff`, `.tif`, `.tiff` | One multi-channel image per sample |
| `channel_names.txt` | text | One marker per line; required for plain TIF or missing QPTIFF metadata |

Place all sample images as a flat list inside the image directory:

```
/path/to/images/
├── sample_001.qptiff
├── sample_002.qptiff
└── channel_names.txt
```

The pipeline's reformat step automatically reorganizes them into the expected layout:

```
/path/to/images/
├── sample_001/
│   └── data/
│       └── sample_001.qptiff
├── sample_002/
│   └── data/
│       └── sample_002.qptiff
└── channel_names.txt
```

## Config File

Create a YAML config from the template at `1_segmentation/conf/config_template.yaml`:

```yaml
image_directory: /path/to/images
panel_path: /path/to/images/channel_names.txt
nuclear_channel: DAPI
segment_channel: PanCK, CD45
segmentation_method: cellpose
```

## Run the Pipeline

The unified pipeline runs all three steps (reformat, segment, quantify) in sequence:

```bash
uv run python3 1_segmentation/src/run_pipeline.py --config path/to/config.yaml
```

To re-run all samples (ignoring previously completed outputs):

```bash
uv run python3 1_segmentation/src/run_pipeline.py --config path/to/config.yaml --force
```

The pipeline logs progress to `segmentation_pipeline.log` inside the image directory and skips samples that already have a mask and quantification CSV.

## Key Outputs

- Masks: `{sample}/full/{method}/{channel}/image_1_seg.npy`
- Quantifications: `quantifications/{method}/{sample}_cell_measurements.csv`

Critical columns in the output CSV:
- `CELL_IDENTIFIER`, `X`, `Y`
- Morphology features (`SIZE`, axis lengths, etc.)
- One intensity column per marker (uppercase marker names)

## Running Individual Steps

The pipeline wraps these scripts, which can also be run standalone:

```bash
uv run python3 1_segmentation/src/reformat_input.py --data_directory /path/to/images
uv run python3 1_segmentation/src/apply_segmentation.py --help
uv run python3 1_segmentation/src/quantify_expression.py --help
```

## Help

```bash
uv run python3 1_segmentation/src/run_pipeline.py --help
```
