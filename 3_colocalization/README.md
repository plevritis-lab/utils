# Colocalization

This module measures spatial co-occurrence between cell types.

- **Local (Python):** per-cell Local CLQ (LCLQ) embeddings
- **Global (R):** per-sample CLQ matrices for condition-level comparison

## Required Inputs

| File | Format | Notes |
|---|---|---|
| Assignment CSVs | `*_assignments.csv` | Must include `X`, `Y`, `FINAL_CELL_TYPE` |
| Signature matrix | CSV | Must include `CELL_TYPE` |
| Colormap | JSON | Needed for local visualization |

## Local Colocalization (Python)

Compute LCLQ embeddings:

```bash
uv run python 3_colocalization/src/local_colocalization/generate_spatial_embeddings.py \
    --data_directory /path/to/assignments \
    --colormap_path /path/to/colormap.json \
    --signature_matrix /path/to/signature_matrix.csv \
    --save_path /path/to/output \
    --bandwidth 100
```

Key output:

```
{save_path}/colocalizations/local_colocalizations/*_local_colocalization_matrix_*.csv
{save_path}/visualizations/cell_colocalizations/
```

Interpretation:
- LCLQ > 1: colocalization (enrichment)
- LCLQ < 1: dispersion (depletion)

Optional clustering (`identify_spatial_clusters.py`) can be run programmatically on saved local matrices.

## Global Colocalization (R)

Compute CLQ summaries per sample:

```bash
Rscript 3_colocalization/src/global_colocalization/generate_condition_summaries.R \
    --data_directory /path/to/assignments \
    --signature_matrix /path/to/signature_matrix.csv \
    --save_path /path/to/colocalizations \
    --bandwidth 100 \
    --n_neighbors 15 \
    --weight_scheme linear
```

Key output:

```
{save_path}/counts/{condition}_count_matrix.csv
{save_path}/global_colocalizations/{condition}_global_colocalization_matrix_*.csv
```

Weighting options: `constant`, `linear`, `squared_exponential`.

## Additional Analyses

- Differential global analysis: `src/global_colocalization/apply_global_colocation_analysis.R`
- Distribution/heatmap scripts:
  - `plot_colocalization_histograms.R`
  - `visualize_colocalizations.R`
  - `epithelial_proportions.R`
- Density-correlation analysis on raw markers:
  - `run_density_correlation_analysis.py`

## Batch Processing

Use templates in `3_colocalization/bash/`:
- `local_colocalization_template.sh`
- `global_colocalization_template.sh`

Workflow:
1. Copy template.
2. Fill `<TODO>` values.
3. Run once per condition.

## Help

```bash
uv run python 3_colocalization/src/local_colocalization/generate_spatial_embeddings.py --help
Rscript 3_colocalization/src/global_colocalization/generate_condition_summaries.R --help
```
