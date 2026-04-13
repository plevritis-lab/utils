import logging
import os
import warnings
from pathlib import Path

import numpy as np
from skimage.segmentation import find_boundaries
from tifffile import imread, imwrite

warnings.filterwarnings(
    "ignore", category=FutureWarning, message=".*torch.load.*weights_only=False.*"
)

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.getLogger("tensorflow").setLevel(logging.ERROR)

logger = logging.getLogger(__name__)


def construct_pseudochannel(
    image: np.ndarray,
    segment_channel: list[int],
) -> np.ndarray:
    """constructs a pseudochannel that merges intensities from a list of markers"""

    pseudochannel = np.sum(image[segment_channel, :, :], axis=0)
    pseudochannel = np.clip(pseudochannel, 0, np.iinfo(image.dtype).max).astype(
        image.dtype
    )
    pseudochannel = pseudochannel[np.newaxis, :, :]

    return np.concatenate((image, pseudochannel), axis=0)


def randomly_sample_patches(
    image: np.ndarray,
    nuclear_channel: int,
    segment_channel: int,
    patch_size: tuple[int, int] = (1000, 1000),
    num_patches: int = 5,
    intensity_threshold: int = 60,
    patch_threshold: float = 0.4,
) -> list[np.ndarray]:
    """selects image regions with enough marker coverage for performance testing"""

    y, x = image.shape[1], image.shape[2]

    upper_y = y - patch_size[0]
    upper_x = x - patch_size[1]

    patches = []

    nuclear_intensity_criterion = np.percentile(
        image[nuclear_channel, :, :], intensity_threshold
    )
    segment_intensity_criterion = np.percentile(
        image[segment_channel, :, :], intensity_threshold
    )

    while len(patches) < num_patches:
        sampled_y = np.random.choice(upper_y)
        sampled_x = np.random.choice(upper_x)

        patch = image[
            :,
            sampled_y : sampled_y + patch_size[0],
            sampled_x : sampled_x + patch_size[1],
        ]

        nuclear_channel_intensities = patch[nuclear_channel, :, :]
        segment_channel_intensities = patch[segment_channel, :, :]

        nuclear_criterion_count = np.count_nonzero(
            nuclear_channel_intensities >= nuclear_intensity_criterion
        )
        segment_criterion_count = np.count_nonzero(
            segment_channel_intensities >= segment_intensity_criterion
        )

        nuclear_criteria = nuclear_criterion_count / (patch_size[0] * patch_size[1])
        segment_criteria = segment_criterion_count / (patch_size[0] * patch_size[1])

        if nuclear_criteria > patch_threshold and segment_criteria > patch_threshold:
            patches.append(patch)

    return patches


def compress_channels(
    images: list[np.ndarray],
    nuclear_channel: int,
    segment_channel: tuple[str, int],
    save_path: str,
) -> list[np.ndarray]:
    """compresses multichannel image into RGB (R=empty, G=segment, B=nuclear)"""

    compressed_images = [None] * len(images)

    segment_channel_name, segment_channel_index = segment_channel

    save_path = os.path.join(save_path, "original")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok=True)

    for i, image in enumerate(images):
        compressed_image = np.zeros(
            (3, image.shape[1], image.shape[2]), dtype=image.dtype
        )
        compressed_image[1, :, :] = image[segment_channel_index, :, :]
        compressed_image[2, :, :] = image[nuclear_channel, :, :]

        imwrite(
            os.path.join(save_path, f"image_{i + 1}.tif"),
            np.transpose(compressed_image, (1, 2, 0)),
        )
        compressed_images[i] = compressed_image

    return compressed_images


def apply_cellpose(
    compressed_images: list[np.ndarray],
    segment_channel_name: str,
    save_path: str,
) -> None:
    """applies cellpose-SAM model to segment compressed RGB images"""

    from cellpose import io, models, utils

    model = models.CellposeModel(pretrained_model="cpsam")

    save_path = os.path.join(save_path, "cellpose")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok=True)

    for i, compressed_image in enumerate(compressed_images):
        masks, flows, _ = model.eval(compressed_image, channels=[2, 3])
        compressed_image[0, :, :] = (
            utils.masks_to_outlines(masks) * np.iinfo(compressed_image.dtype).max
        )

        imwrite(
            os.path.join(save_path, f"image_{i + 1}.tif"),
            np.transpose(compressed_image, (1, 2, 0)),
        )
        io.masks_flows_to_seg(
            compressed_image,
            masks,
            flows,
            os.path.join(save_path, f"image_{i + 1}"),
            channels=[2, 3],
        )


def apply_mesmer(
    compressed_images: list[np.ndarray],
    segment_channel_name: str,
    save_path: str,
) -> None:
    """applies deepcell mesmer model to segment compressed RGB images"""

    from deepcell.applications import Mesmer  # type: ignore

    model = Mesmer()

    save_path = os.path.join(save_path, "mesmer")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok=True)

    for i, compressed_image in enumerate(compressed_images):
        reshaped_image = np.transpose(compressed_image, (2, 1, 0)).copy()

        reshaped_image = reshaped_image[:, :, 1:]
        reshaped_image = reshaped_image[..., ::-1]

        reshaped_image = np.expand_dims(reshaped_image, axis=0)

        segmentation_predictions = model.predict(
            reshaped_image, batch_size=1, image_mpp=0.377, compartment="whole-cell"
        )[0, :, :, 0]
        segmentation_predictions = np.transpose(segmentation_predictions, (1, 0))

        reshaped_image = reshaped_image[0, :, :, :]
        reshaped_image = np.stack(
            (
                np.zeros(
                    (reshaped_image.shape[0], reshaped_image.shape[1]),
                    dtype=reshaped_image.dtype,
                ),
                reshaped_image[:, :, 1],
                reshaped_image[:, :, 0],
            ),
            axis=-1,
        )

        reshaped_image = np.transpose(reshaped_image, (1, 0, 2))

        located_boundaries = find_boundaries(
            segmentation_predictions, connectivity=1, mode="inner"
        )
        reshaped_image[:, :, 0][located_boundaries > 0] = np.iinfo(
            reshaped_image.dtype
        ).max

        imwrite(os.path.join(save_path, f"image_{i + 1}.tif"), reshaped_image)
        np.save(
            os.path.join(save_path, f"image_{i + 1}_seg.npy"),
            segmentation_predictions,
        )


def overlay_masks(
    compressed_images: list[np.ndarray],
    cellpose_segment_channel_name: str,
    mesmer_segment_channel_name: str,
    display_segment_channel_name: str,
    load_path: str,
    save_path: str,
) -> None:
    """overlays cellpose and mesmer masks onto the same image for comparison"""

    save_path = os.path.join(save_path, "combined_algorithm_overlays")
    os.makedirs(save_path, exist_ok=True)

    for i, compressed_image in enumerate(compressed_images):
        overlaid_image = compressed_image.copy()

        try:
            cellpose_segmentation = np.load(
                os.path.join(
                    load_path,
                    f"cellpose/{cellpose_segment_channel_name}/image_{i + 1}_seg.npy",
                ),
                allow_pickle=True,
            ).item()["outlines"]
            mesmer_segmentation = np.load(
                os.path.join(
                    load_path,
                    f"mesmer/{mesmer_segment_channel_name}/image_{i + 1}_seg.npy",
                ),
                allow_pickle=True,
            )

        except FileNotFoundError:
            logger.error(
                "please double-check that your provided segmentation masks are in "
                "their specified locations\n"
                "currently using cellpose segmentation mask path: %s\n"
                "currently using mesmer segmentation mask path: %s",
                os.path.join(
                    load_path,
                    f"cellpose/{cellpose_segment_channel_name}/image_{i + 1}_seg.npy",
                ),
                os.path.join(
                    load_path,
                    f"mesmer/{mesmer_segment_channel_name}/image_{i + 1}_seg.npy",
                ),
            )
            return

        mesmer_boundaries = find_boundaries(
            mesmer_segmentation, connectivity=1, mode="inner"
        )

        yellow_mask = cellpose_segmentation > 0
        red_mask = mesmer_boundaries > 0

        overlaid_image[0, :, :][yellow_mask] = np.iinfo(overlaid_image.dtype).max
        overlaid_image[1, :, :][yellow_mask] = np.iinfo(overlaid_image.dtype).max
        overlaid_image[0, :, :][red_mask] = np.iinfo(overlaid_image.dtype).max
        overlaid_image[1, :, :][red_mask & yellow_mask] = int(
            np.iinfo(overlaid_image.dtype).max * 2 / 3
        )

        imwrite(
            os.path.join(
                save_path,
                f"image_{i + 1}_cellpose_{cellpose_segment_channel_name}"
                f"_mesmer_{mesmer_segment_channel_name}"
                f"_display_{display_segment_channel_name}.tif",
            ),
            np.transpose(overlaid_image, (1, 2, 0)),
        )


def segment_sample(
    image_path: Path,
    nuclear_channel: str,
    segment_channel: list[str],
    segmentation_method: str,
    panel: list[str],
) -> Path:
    """segments a single sample image and returns the path to the generated mask"""

    np.random.seed(42)

    nuclear_channel_index = panel.index(nuclear_channel)
    segment_channel_index = [panel.index(s) for s in segment_channel]

    image = imread(str(image_path))

    if len(segment_channel) > 1:
        image = construct_pseudochannel(image, segment_channel_index)
        segment_channel_name = " + ".join(segment_channel)
        segment_channel_final_index = image.shape[0] - 1
    else:
        segment_channel_name = segment_channel[0]
        segment_channel_final_index = segment_channel_index[0]

    base_directory = image_path.parent.parent
    save_path = str(base_directory / "full")

    compressed_images = compress_channels(
        [image],
        nuclear_channel_index,
        (segment_channel_name, segment_channel_final_index),
        save_path,
    )

    if segmentation_method == "cellpose":
        apply_cellpose(compressed_images, segment_channel_name, save_path)
    elif segmentation_method == "mesmer":
        apply_mesmer(compressed_images, segment_channel_name, save_path)

    mask_path = (
        base_directory
        / "full"
        / segmentation_method
        / segment_channel_name
        / "image_1_seg.npy"
    )

    return mask_path
