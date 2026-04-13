# Segmentation

This module detects cells from multiplexed CODEX images and quantifies per-cell marker expression.

- **Step 1:** segmentation mask generation (cellpose or mesmer)
- **Step 2:** per-cell measurement extraction into CSV files

## Prerequisites

- Standard install: `uv sync`
- HPC/Linux extras (required for mesmer): `uv sync --extra hpc`
- Segmentation is compute-heavy; use HPC for full-resolution runs.

## Required Inputs

| File | Format | Notes |
|---|---|---|
| Proteomic image | `.qptiff`, `.tif`, `.tiff` | One multi-channel image per sample |
| `channel_names.txt` | text | One marker per line; required for plain TIF or missing QPTIFF metadata |

Expected sample layout:

```
/path/to/images/
└── sample_001/
    └── data/
        └── sample_001.qptiff
```

If your images are flat, reorganize first:

```bash
uv run python 1_segmentation/src/reformat_input.py --data_directory /path/to/images
```

## Run Segmentation

```bash
uv run python 1_segmentation/src/apply_segmentation.py \
    --image_path /path/to/sample_001/data/sample_001.qptiff \
    --panel_path /path/to/channel_names.txt \
    --nuclear_channel DAPI \
    --segment_channel PanCK,CD45 \
    --apply_cellpose
```

Use one of:
- `--apply_cellpose`
- `--apply_mesmer` (requires `DEEPCELL_ACCESS_TOKEN`)

Useful options:
- `--overlay_masks` to compare cellpose vs mesmer
- `--debug` to run on random high-coverage patches

Key output files:
- mask: `{sample}/full/{method}/{channel}/image_1_seg.npy`
- quick image: `{sample}/full/{method}/{channel}/image_1.tif`

## Run Quantification

```bash
uv run python 1_segmentation/src/quantify_expression.py \
    --image_path /path/to/sample_001/data/sample_001.qptiff \
    --mask_path /path/to/sample_001/full/cellpose/PanCK/image_1_seg.npy \
    --panel_path /path/to/channel_names.txt \
    --save_path /path/to/images \
    --apply_cellpose
```

Output:

```
/path/to/images/quantifications/{method}/{sample}_cell_measurements.csv
```

Critical columns in output CSV:
- `CELL_IDENTIFIER`, `X`, `Y`
- morphology features (`SIZE`, axis lengths, etc.)
- one intensity column per marker (uppercase marker names)

## Batch Processing

Use templates in `1_segmentation/bash/`:
- `segmentation_template.sh`
- `quantify_expressions_template.sh`

Workflow:
1. Copy template.
2. Replace each `<TODO>`.
3. Run script.

## Help

```bash
uv run python 1_segmentation/src/apply_segmentation.py --help
uv run python 1_segmentation/src/quantify_expression.py --help
```
