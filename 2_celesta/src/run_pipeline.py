import argparse
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generate_threshold_report import generate_threshold_report
from generate_thresholds import process_thresholds
from visualize_assignments import visualize_assignments, visualize_cell_proportions

logger = logging.getLogger(__name__)


def generate_colormap(
    signature_matrix_path: Path,
    save_path: Path,
) -> dict:
    """auto-generates a colormap from signature matrix cell types using a perceptually distinct palette"""

    signature_matrix = pd.read_csv(signature_matrix_path)
    cell_types = signature_matrix["CELL_TYPE"].tolist()

    palette = plt.colormaps.get_cmap("tab20").resampled(max(len(cell_types), 20))
    colormap = {}

    for i, cell_type in enumerate(cell_types):
        rgb = palette(i % 20)[:3]
        hex_color = "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        )
        colormap[cell_type] = {"color": hex_color}

    colormap["unknown"] = {"color": "#7F7F7F"}

    with open(save_path, "w") as f:
        json.dump(colormap, f, indent=4)

    return colormap


def run_pipeline(
    config_path: Path,
    sample_filter: list[str] | None = None,
) -> None:
    """orchestrates the full celesta pipeline: thresholds, assignment, visualization"""

    with open(config_path) as f:
        config = yaml.safe_load(f)

    image_directory = Path(config["image_directory"])
    signature_matrix = Path(config["signature_matrix"])
    panel_path = Path(config["panel_path"])
    segmentation_method = config["segmentation_method"]
    colormap_path = Path(config["colormap"]) if "colormap" in config else None

    with panel_path.open() as f:
        panel = [marker.rstrip() for marker in f]

    log_path = image_directory / "celesta_pipeline.log"
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.handlers.clear()
    root_logger.addHandler(logging.StreamHandler(sys.stdout))
    root_logger.addHandler(logging.FileHandler(log_path, mode="a"))
    for handler in root_logger.handlers:
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s — %(message)s")
        )

    logger.info("")
    logger.info("=" * 80)
    logger.info("CELESTA PIPELINE RUN — %s", datetime.now().isoformat())
    logger.info("=" * 80)
    logger.info("config: %s", config)

    if sample_filter:
        logger.info("filter: %s", sample_filter)

    thresholds_directory = image_directory / "thresholds"
    assignments_directory = image_directory / "assignments"

    logger.info("[thresholds] generating thresholds")
    process_thresholds(
        image_directory=image_directory,
        signature_matrix_path=signature_matrix,
        save_path=thresholds_directory,
        sample_filter=sample_filter,
    )
    logger.info("[thresholds] done")

    quantifications_directory = (
        image_directory / "quantifications" / segmentation_method
    )
    celesta_script = Path(__file__).resolve().parent / "apply_celesta.R"

    # fmt: off
    celesta_command = [
        "Rscript", str(celesta_script), "--data_directory", str(quantifications_directory), "--save_path", str(assignments_directory),
            "--signature_matrix", str(signature_matrix), "--thresholds_directory", str(thresholds_directory)
    ]
    # fmt: on

    if sample_filter:
        celesta_command.extend(["--filter", ",".join(sample_filter)])

    logger.info("[celesta] running: %s", " ".join(celesta_command))
    result = subprocess.run(
        celesta_command,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        logger.info("[celesta] stdout:\n%s", result.stdout)

    if result.returncode != 0:
        logger.error("[celesta] failed with exit code %d", result.returncode)
        if result.stderr:
            logger.error("[celesta] stderr:\n%s", result.stderr)
        raise RuntimeError(
            f"Rscript apply_celesta.R failed with exit code {result.returncode}"
        )

    logger.info("[celesta] done")

    if colormap_path and colormap_path.is_file():
        with open(colormap_path) as f:
            cell_type_info = json.load(f)
        logger.info("[colormap] loaded from %s", colormap_path)
    else:
        generated_path = image_directory / "colormap.json"
        cell_type_info = generate_colormap(signature_matrix, generated_path)
        logger.info("[colormap] auto-generated at %s", generated_path)

    cell_plots_directory = assignments_directory / "cell_plots"
    cell_proportions_directory = assignments_directory / "cell_proportions"
    cell_plots_directory.mkdir(parents=True, exist_ok=True)
    cell_proportions_directory.mkdir(parents=True, exist_ok=True)

    display_cells = list(cell_type_info.keys())

    assignment_files = sorted(assignments_directory.glob("*_assignments.csv"))

    if sample_filter:
        assignment_files = [
            f
            for f in assignment_files
            if f.stem.replace("_assignments", "") in sample_filter
        ]

    for assignment_file in assignment_files:
        sample_name = assignment_file.stem.replace("_assignments", "")
        logger.info("[visualize] %s", sample_name)

        assignments = pd.read_csv(assignment_file).dropna()

        visualize_assignments(
            assignments,
            cell_type_info,
            display_cells,
            str(cell_plots_directory),
            sample_name,
        )
        visualize_cell_proportions(
            assignments,
            cell_type_info,
            str(cell_proportions_directory),
            sample_name,
        )

    IMAGE_EXTENSIONS = {".qptiff", ".tif", ".tiff"}
    signature_matrix_df = pd.read_csv(signature_matrix)

    for assignment_file in assignment_files:
        sample_name = assignment_file.stem.replace("_assignments", "")
        sample_data_directory = image_directory / sample_name / "data"

        if not sample_data_directory.is_dir():
            logger.warning("[report] no data directory for %s, skipping", sample_name)
            continue

        image_files = [
            f for f in sample_data_directory.iterdir() if f.suffix in IMAGE_EXTENSIONS
        ]

        if not image_files:
            logger.warning("[report] no image found for %s, skipping", sample_name)
            continue

        threshold_path = thresholds_directory / f"{sample_name}_thresholds.csv"
        if not threshold_path.is_file():
            logger.warning("[report] no thresholds for %s, skipping", sample_name)
            continue

        logger.info("[report] %s", sample_name)
        report_path = generate_threshold_report(
            assignments_path=assignment_file,
            image_path=image_files[0],
            panel=panel,
            signature_matrix=signature_matrix_df,
            thresholds=pd.read_csv(threshold_path),
            colormap=cell_type_info,
            sample_name=sample_name,
            save_path=assignments_directory,
        )
        logger.info("[report] saved to %s", report_path)

    logger.info("pipeline complete — processed %d samples", len(assignment_files))


def main() -> None:
    # fmt: off
    parser = argparse.ArgumentParser(description="unified celesta pipeline: thresholds, assignment, visualization")
    parser.add_argument("--config", required=True, help="path to YAML config file")
    parser.add_argument("--filter", default=None, help="comma-separated sample names to process")
    # fmt: on
    arguments = parser.parse_args()

    sample_filter = (
        [s.strip() for s in arguments.filter.split(",")] if arguments.filter else None
    )

    run_pipeline(config_path=Path(arguments.config), sample_filter=sample_filter)


if __name__ == "__main__":
    main()
