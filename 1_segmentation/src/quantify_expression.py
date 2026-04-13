from pathlib import Path

import numpy as np
import pandas as pd
from skimage.measure import regionprops
from tifffile import imread


def quantify_expression(
    image: np.ndarray,
    segmentation_mask: np.ndarray,
    panel: list[str],
    save_path: str,
) -> None:
    """quantifies expression of every marker in the panel per cell"""

    cell_data = {}

    for channel_number in range(image.shape[0]):
        channel_image = image[channel_number, :, :]

        properties = regionprops(segmentation_mask, channel_image)

        for cell_property in properties:
            cell_identifier = cell_property.label

            if cell_identifier not in cell_data:
                cell_data[cell_identifier] = {"CELL_IDENTIFIER": cell_identifier}

            cell_data[cell_identifier]["MAJOR_AXIS_LENGTH"] = (
                cell_property.axis_major_length
            )
            cell_data[cell_identifier]["MINOR_AXIS_LENGTH"] = (
                cell_property.axis_minor_length
            )

            cell_centroid_y, cell_centroid_x = cell_property.centroid
            cell_data[cell_identifier]["X"] = round(cell_centroid_x)
            cell_data[cell_identifier]["Y"] = round(cell_centroid_y)
            cell_data[cell_identifier]["SIZE"] = round(cell_property.num_pixels)
            cell_data[cell_identifier]["ECCENTRICITY"] = cell_property.eccentricity
            cell_data[cell_identifier]["ORIENTATION"] = cell_property.orientation

            cell_data[cell_identifier][panel[channel_number].upper()] = (
                cell_property.intensity_mean
            )

    cell_expressions = pd.DataFrame.from_dict(cell_data, orient="index")
    cell_expressions.to_csv(f"{save_path}_cell_measurements.csv", index=False)


def quantify_sample(
    image_path: Path,
    mask_path: Path,
    panel: list[str],
    segmentation_method: str,
    save_directory: Path,
) -> Path:
    """loads image and mask, quantifies expression, and returns path to output CSV"""

    image = imread(str(image_path))

    if segmentation_method == "cellpose":
        segmentation_mask = np.load(str(mask_path), allow_pickle=True).item()[
            "outlines"
        ]
    else:
        segmentation_mask = np.load(str(mask_path), allow_pickle=True)

    output_directory = save_directory / "quantifications" / segmentation_method
    output_directory.mkdir(parents=True, exist_ok=True)

    sample_name = image_path.stem
    save_path = str(output_directory / sample_name)

    quantify_expression(image, segmentation_mask, panel, save_path)

    return Path(f"{save_path}_cell_measurements.csv")
