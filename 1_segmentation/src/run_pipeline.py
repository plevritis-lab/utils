import argparse
import logging
import sys
import traceback
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml

from apply_segmentation import segment_sample
from quantify_expression import quantify_sample
from reformat_input import prepare_input

logger = logging.getLogger(__name__)

IMAGE_EXTENSIONS = {".qptiff", ".tif", ".tiff"}


def run_pipeline(
    config_path: Path,
    force: bool = False,
) -> None:
    """orchestrates the full segmentation pipeline: reformat, segment, quantify"""

    with open(config_path) as f:
        config = yaml.safe_load(f)

    image_directory = Path(config["image_directory"])
    panel_path = Path(config["panel_path"])
    segment_channel = [s.strip() for s in config["segment_channel"].split(",")]
    segmentation_method = config.get("segmentation_method", "cellpose")

    with panel_path.open() as f:
        panel = [marker.rstrip() for marker in f]

    log_path = image_directory / "segmentation_pipeline.log"
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
    logger.info("SEGMENTATION PIPELINE RUN — %s", datetime.now().isoformat())
    logger.info("=" * 80)
    logger.info("config: %s", config)

    prepare_input(image_directory)

    segment_channel_name = (
        " + ".join(segment_channel) if len(segment_channel) > 1 else segment_channel[0]
    )

    succeeded, failed, skipped = [], [], []

    for sample_directory in sorted(image_directory.iterdir()):
        if not sample_directory.is_dir():
            continue

        sample_name = sample_directory.name
        image_files = (
            [
                f
                for f in (sample_directory / "data").iterdir()
                if f.suffix in IMAGE_EXTENSIONS
            ]
            if (sample_directory / "data").is_dir()
            else []
        )

        if not image_files:
            continue

        image_path = image_files[0]

        mask_path = (
            sample_directory
            / "full"
            / segmentation_method
            / segment_channel_name
            / "image_1_seg.npy"
        )
        csv_path = (
            image_directory
            / "quantifications"
            / segmentation_method
            / f"{sample_name}_cell_measurements.csv"
        )

        if not force and mask_path.is_file() and csv_path.is_file():
            logger.info("[skip] %s", sample_name)
            skipped.append(sample_name)
            continue

        logger.info("[start] %s", sample_name)

        try:
            generated_mask = segment_sample(
                image_path=image_path,
                nuclear_channel=config["nuclear_channel"],
                segment_channel=segment_channel,
                segmentation_method=segmentation_method,
                panel=panel,
            )
            quantify_sample(
                image_path=image_path,
                mask_path=generated_mask,
                panel=panel,
                segmentation_method=segmentation_method,
                save_directory=image_directory,
            )
            logger.info("[done] %s", sample_name)
            succeeded.append(sample_name)
        except Exception:
            logger.error("[fail] %s\n%s", sample_name, traceback.format_exc())
            failed.append(sample_name)

    logger.info(
        "pipeline complete — succeeded: %d, skipped: %d, failed: %d",
        len(succeeded),
        len(skipped),
        len(failed),
    )
    if failed:
        logger.error("failed samples: %s", failed)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="unified segmentation pipeline: reformat, segment, quantify"
    )
    parser.add_argument("--config", required=True, help="path to YAML config file")
    parser.add_argument("--force", action="store_true", help="re-run all samples")
    arguments = parser.parse_args()
    run_pipeline(config_path=Path(arguments.config), force=arguments.force)


if __name__ == "__main__":
    main()
