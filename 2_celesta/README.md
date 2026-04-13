# CELESTA

This module assigns cell types from segmentation quantification CSVs and provides interactive/static visualization.

## Required Inputs

| File | Format | Notes |
|---|---|---|
| Cell measurements | `*_cell_measurements.csv` | Output from `1_segmentation` |
| Signature matrix | CSV | Must include `CELL_TYPE` |
| Colormap | JSON | Keys should match final cell type labels; include `"unknown"` |
| `channel_names.txt` | text | Needed for QC generation |
| `clinical_data.csv` | CSV (optional) | For QC sheet; needs `TMA`, `TMA_PART`, `CORE_IMAGE_ID` |

Templates:
- colormap: `2_celesta/src/colormaps/colormap_template.json`
- batch script: `2_celesta/bash/celesta_processing_template.sh`

## Typical Run Order

1. Generate thresholds
2. Run CELESTA assignment (R)
3. Review overlays/plots
4. Tune thresholds and rerun

## 1) Generate Thresholds

```bash
uv run python 2_celesta/src/generate_thresholds.py \
    --image_directory /path/to/images \
    --signature_matrix /path/to/signature_matrix.csv \
    --save_path /path/to/thresholds
```

Output: `{sample}_thresholds.csv` per sample.

## 2) Run CELESTA Assignments (R)

```bash
Rscript 2_celesta/src/apply_celesta.R \
    --data_directory /path/to/quantifications/cellpose \
    --save_path /path/to/assignments \
    --signature_matrix /path/to/signature_matrix.csv \
    --thresholds_directory /path/to/thresholds
```

Output: `{sample}_assignments.csv` per sample.

Critical output columns:
- `CELL_IDENTIFIER`, `X`, `Y`
- `{MARKER}_PROBABILITY`
- `FINAL_CELL_TYPE`

## 3) Visualize Assignments

Interactive (napari):

```bash
uv run python 2_celesta/src/visualize_dynamic_overlays.py \
    --assignments_path /path/to/sample_001_assignments.csv \
    --colormap_path /path/to/colormap.json \
    --image_path /path/to/sample_001/data/sample_001.qptiff
```

Static plots:

```bash
uv run python 2_celesta/src/visualize_assignments.py \
    --assignments_path /path/to/sample_001_assignments.csv \
    --colormap_path /path/to/colormap.json \
    --display_cells all \
    --save_path /path/to/visualizations
```

Key outputs:
- `visualizations/cell_plots/*_assignments.png`
- `visualizations/cell_proportions/*_proportions.png`

## Optional: QC Spreadsheet

```bash
uv run python 2_celesta/src/generate_quality_control.py \
    --clinical_data /path/to/clinical_data.csv \
    --panel_path /path/to/channel_names.txt \
    --save_path /path/to/output
```

Output: `{save_path}_quality_control.csv`

## Batch Processing

Use `2_celesta/bash/celesta_processing_template.sh`:
1. Copy template.
2. Replace all `<TODO>` values.
3. Run once per condition.

## Help

```bash
uv run python 2_celesta/src/generate_thresholds.py --help
uv run python 2_celesta/src/visualize_assignments.py --help
Rscript 2_celesta/src/apply_celesta.R --help
```
