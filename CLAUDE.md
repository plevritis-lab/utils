# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EPOCHS is a collection of utility scripts for processing CODEX (multiplexed immunohistochemistry) spatial proteomics data. It is organized as independent modules, each with `source/` (Python/R scripts) and `bash/` (batch processing templates) subdirectories. Some modules also use R.

## Development Setup

```bash
uv sync                        # install all dependencies from lockfile
uv sync --extra hpc            # also install gudhi (Linux HPC only)
uv run pre-commit install      # set up pre-commit hooks
```

`deepcell` (Mesmer) requires numpy<2 and TensorFlow 2.8 -- it must be installed manually in a separate HPC environment (`pip install deepcell`). Cellpose works on Apple M-chips; Mesmer does not.

## Common Commands

```bash
uv run pytest                  # run tests
uv run ruff check .            # lint
uv run ruff format .           # format
uv add <package>               # add a runtime dependency
uv add --dev <package>         # add a dev dependency
```

## Running Scripts

Scripts are invoked via `uv run` and use `argparse` for CLI arguments:

```bash
uv run python segmentation/source/apply_segmentation.py --help
uv run python colocalization/source/local_colocalization/generate_spatial_embeddings.py --help
```

Batch processing uses bash templates (`bash/*_template.sh`) with `<TODO>` placeholders that are filled in for each project. Project-specific scripts live in `bash/<project_name>/` subdirectories. Segmentation is computationally expensive and should be run on an HPC cluster.

## Architecture

Each module is independent with clear data contracts between them:

```
Raw QPTIFF Images
  -> segmentation (cellpose/mesmer) -> cell masks (.npy) + measurements (.csv)
  -> celesta (R) -> cell type assignments (.csv)
  -> colocalization (Python + R) / topology (Python + GUDHI) / differential_testing (Python)
  -> visualization (matplotlib/napari)
```

### Module Layout

- **segmentation**: cellpose/mesmer wrappers for cell detection and per-cell measurement quantification
- **celesta**: R-based cell type assignment + napari interactive visualization of overlays
- **colocalization**: global (R, metadisco package) and local (Python, LCLQ) spatial co-occurrence analysis
- **topology**: persistent homology feature extraction via GUDHI alpha complexes
- **differential_testing**: statistical comparisons (PERMANOVA, Mann-Whitney U, Bray-Curtis) across conditions
- **visualization**: ROI selection and publication figure generation
- **pathology**: histology image post-processing

### Key Conventions

- Image arrays use shape `(channels, height, width)`
- Cell data CSVs use columns: `CELL_IDENTIFIER`, `X`, `Y`, `FINAL_CELL_TYPE`
- Colormaps are stored as JSON files mapping cell types to hex colors
- Python scripts use `argparse`; bash templates loop over sample directories
- Imports reference modules as packages: `from segmentation.source.utils import extract_proteomic_panel`
