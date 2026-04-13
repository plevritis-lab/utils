import argparse
import logging
import sys
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml

from apply_segmentation import segment_sample
from quantify_expression import quantify_sample
from reformat_input import prepare_input
from utils import extract_proteomic_panel

logger = logging.getLogger(__name__)

SKIP_DIRECTORIES = {"histology", "quantifications", "assignments"}
IMAGE_EXTENSIONS = {".qptiff", ".tif", ".tiff"}
REQUIRED_CONFIG_KEYS = {
    "image_directory",
    "nuclear_channel",
    "segment_channel",
    "segmentation_method",
}


def load_config(config_path: Path) -> dict:
    """loads and validates the pipeline configuration from a YAML file"""

    with open(config_path) as f:
        config = yaml.safe_load(f)

    missing_keys = REQUIRED_CONFIG_KEYS - set(config.keys())
    if missing_keys:
        raise ValueError(
            f"missing required config keys: {', '.join(sorted(missing_keys))}"
        )

    if config["segmentation_method"] not in ("cellpose", "mesmer"):
        raise ValueError(
            f"segmentation_method must be 'cellpose' or 'mesmer', "
            f"got '{config['segmentation_method']}'"
        )

    config["image_directory"] = Path(config["image_directory"])
    if not config["image_directory"].is_dir():
        raise FileNotFoundError(
            f"image_directory does not exist: {config['image_directory']}"
        )

    if "panel_path" in config:
        config["panel_path"] = Path(config["panel_path"])
        if not config["panel_path"].is_file():
            raise FileNotFoundError(
                f"panel_path does not exist: {config['panel_path']}"
            )

    config["segment_channel"] = [
        s.strip() for s in config["segment_channel"].split(",")
    ]

    return config


def detect_image_extension(image_directory: Path) -> str:
    """auto-detects the image file extension used in sample subdirectories"""

    for sample_directory in sorted(image_directory.iterdir()):
        if not sample_directory.is_dir() or sample_directory.name in SKIP_DIRECTORIES:
            continue

        data_directory = sample_directory / "data"
        if not data_directory.is_dir():
            continue

        for file in data_directory.iterdir():
            if file.suffix in IMAGE_EXTENSIONS:
                return file.suffix

    raise FileNotFoundError(
        f"no image files with extensions {IMAGE_EXTENSIONS} found in {image_directory}"
    )


def discover_samples(
    image_directory: Path,
    image_extension: str,
) -> list[Path]:
    """discovers sample image paths within the expected directory structure"""

    samples = []

    for sample_directory in sorted(image_directory.iterdir()):
        if not sample_directory.is_dir() or sample_directory.name in SKIP_DIRECTORIES:
            continue

        sample_name = sample_directory.name
        image_path = sample_directory / "data" / f"{sample_name}{image_extension}"

        if image_path.is_file():
            samples.append(image_path)

    return samples


def sample_is_complete(
    image_path: Path,
    segmentation_method: str,
    segment_channel: list[str],
) -> bool:
    """checks whether segmentation and quantification outputs already exist for a sample"""

    base_directory = image_path.parent.parent
    sample_name = image_path.stem

    segment_channel_name = (
        " + ".join(segment_channel) if len(segment_channel) > 1 else segment_channel[0]
    )

    mask_path = (
        base_directory
        / "full"
        / segmentation_method
        / segment_channel_name
        / "image_1_seg.npy"
    )

    csv_directory = (
        image_path.parent.parent.parent / "quantifications" / segmentation_method
    )
    csv_path = csv_directory / f"{sample_name}_cell_measurements.csv"

    return mask_path.is_file() and csv_path.is_file()


def extract_panel(
    image_path: Path,
    panel_path: Path | None,
) -> list[str]:
    """extracts the protein panel, falling back to panel_path if metadata extraction fails"""

    try:
        return extract_proteomic_panel(str(image_path))
    except (AttributeError, KeyError, ET.ParseError):
        if panel_path is not None:
            logger.info("metadata extraction failed, using panel_path: %s", panel_path)
            return extract_proteomic_panel(str(image_path), str(panel_path))

        raise RuntimeError(
            "unable to extract marker metadata from image file; "
            "add 'panel_path' to your config.yaml pointing to a channel_names.txt file"
        )


def setup_logging(image_directory: Path) -> None:
    """configures logging to both console and a log file in the image directory"""

    log_path = image_directory / "pipeline.log"

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    root_logger.handlers.clear()

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    root_logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_path, mode="a")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s — %(message)s")
    )
    root_logger.addHandler(file_handler)


def run_pipeline(config_path: Path, force: bool = False) -> None:
    """orchestrates the full segmentation pipeline: reformat, segment, quantify"""

    config = load_config(config_path)
    image_directory = config["image_directory"]

    setup_logging(image_directory)

    logger.info("=" * 60)
    logger.info("segmentation pipeline started at %s", datetime.now().isoformat())
    logger.info("config: %s", config)
    logger.info("force: %s", force)
    logger.info("=" * 60)

    logger.info("reformatting input directory if needed")
    prepare_input(image_directory)

    image_extension = detect_image_extension(image_directory)
    logger.info("detected image extension: %s", image_extension)

    samples = discover_samples(image_directory, image_extension)
    logger.info("discovered %d samples", len(samples))

    if not samples:
        logger.warning("no samples found in %s", image_directory)
        return

    panel_path = config.get("panel_path")
    panel = None

    succeeded = []
    failed = []
    skipped = []

    for image_path in samples:
        sample_name = image_path.stem

        if not force and sample_is_complete(
            image_path, config["segmentation_method"], config["segment_channel"]
        ):
            logger.info("[skip] %s — outputs already exist", sample_name)
            skipped.append(sample_name)
            continue

        logger.info("[start] %s", sample_name)

        try:
            if panel is None:
                panel = extract_panel(image_path, panel_path)

            mask_path = segment_sample(
                image_path=image_path,
                nuclear_channel=config["nuclear_channel"],
                segment_channel=config["segment_channel"],
                segmentation_method=config["segmentation_method"],
                panel=panel,
            )

            quantify_sample(
                image_path=image_path,
                mask_path=mask_path,
                panel=panel,
                segmentation_method=config["segmentation_method"],
                save_directory=image_directory,
            )

            logger.info("[done] %s", sample_name)
            succeeded.append(sample_name)

        except Exception:
            logger.error("[fail] %s\n%s", sample_name, traceback.format_exc())
            failed.append(sample_name)

    logger.info("=" * 60)
    logger.info("pipeline complete")
    logger.info("  succeeded: %d  %s", len(succeeded), succeeded)
    logger.info("  skipped:   %d  %s", len(skipped), skipped)
    logger.info("  failed:    %d  %s", len(failed), failed)
    logger.info("=" * 60)

    if failed:
        logger.error(
            "the following samples failed — check pipeline.log for tracebacks: %s",
            failed,
        )


def parse_arguments() -> argparse.Namespace:
    """parses command line arguments for the segmentation pipeline"""

    parser = argparse.ArgumentParser(
        description="unified segmentation pipeline: reformat, segment, quantify"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="path to YAML configuration file",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-run all samples even if outputs already exist",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    run_pipeline(
        config_path=Path(arguments.config),
        force=arguments.force,
    )


if __name__ == "__main__":
    main()
