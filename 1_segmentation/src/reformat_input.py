import shutil
from pathlib import Path

import numpy as np
from tifffile import TiffFile, imwrite

IMAGE_EXTENSIONS = (".tiff", ".tif", ".qptiff")


def prepare_input(data_directory: Path) -> None:
    """reorganizes flat image files into sample/data/sample.ext subdirectories"""

    image_files = [f for f in data_directory.iterdir() if f.suffix in IMAGE_EXTENSIONS]

    if not image_files:
        return

    for image_file in image_files:
        sample_name = image_file.stem
        sample_directory = data_directory / sample_name / "data"
        sample_directory.mkdir(parents=True, exist_ok=True)
        shutil.move(str(image_file), str(sample_directory / image_file.name))


def reformat_image(image_path: Path) -> None:
    """reformats an image by transposing its dimensions to (c, y, x), overwriting the original"""

    with TiffFile(str(image_path)) as tif:
        image = tif.asarray()
        axes = tif.series[0].axes

    assert all(a in axes for a in ("C", "X", "Y"))

    if "Z" in axes:
        image = np.moveaxis(
            image,
            [axes.index("Z"), axes.index("C"), axes.index("Y"), axes.index("X")],
            [0, 1, 2, 3],
        )
        image = image.reshape(-1, image.shape[2], image.shape[3])
    else:
        image = np.moveaxis(
            image,
            [axes.index("C"), axes.index("Y"), axes.index("X")],
            [0, 1, 2],
        )

    imwrite(str(image_path), image)
