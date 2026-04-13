import numpy as np
from pathlib import Path
from tifffile import imwrite

from reformat_input import prepare_input


def _write_dummy_tif(path: Path, shape: tuple = (3, 64, 64)) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    imwrite(str(path), np.zeros(shape, dtype=np.uint16))


def test_prepare_input_moves_flat_images(tmp_path):
    _write_dummy_tif(tmp_path / "sample_001.tif")
    _write_dummy_tif(tmp_path / "sample_002.tiff")

    prepare_input(tmp_path)

    assert (tmp_path / "sample_001" / "data" / "sample_001.tif").is_file()
    assert (tmp_path / "sample_002" / "data" / "sample_002.tiff").is_file()
    assert not (tmp_path / "sample_001.tif").exists()
    assert not (tmp_path / "sample_002.tiff").exists()


def test_prepare_input_is_noop_when_already_structured(tmp_path):
    sample_dir = tmp_path / "sample_001" / "data"
    _write_dummy_tif(sample_dir / "sample_001.tif")

    prepare_input(tmp_path)

    assert (sample_dir / "sample_001.tif").is_file()
    assert len(list(tmp_path.iterdir())) == 1
