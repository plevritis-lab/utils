# Topology

This module computes persistent homology features from cell coordinates to summarize spatial structure by cell type.

## Prerequisite

`gudhi` is Linux-only in this workflow:

```bash
uv sync --extra hpc
```

## Required Inputs

| File | Format | Notes |
|---|---|---|
| Assignment CSVs | `*_assignments.csv` | Must include `X`, `Y`, `FINAL_CELL_TYPE` |
| Signature matrix | CSV | Must include `CELL_TYPE` |

## Run

```bash
uv run python 4_topology/src/apply_persistent_homology.py \
    --data_directory /path/to/assignments \
    --signature_matrix /path/to/signature_matrix.csv \
    --save_path /path/to/output \
    --count_critera 3
```

Key option:
- `--count_critera`: minimum cells required to compute persistence for a cell type

## Outputs

```
{save_path}/topology/{sample}_persistence.json
{save_path}/visualizations/persistence_diagrams/{sample}/{cell_type}_persistence_diagram_density.png
```

Interpretation:
- **H0** features describe connected components (clusters)
- **H1** features describe loop-like structures
- Long-lived features (death - birth) indicate stronger spatial structure

Cell types below the count threshold are recorded as null/NA entries.

## Batch Processing

Use `4_topology/bash/topology_template.sh`:
1. Copy template.
2. Fill `<TODO>` values.
3. Run once per condition.

## Help

```bash
uv run python 4_topology/src/apply_persistent_homology.py --help
```
