import numpy as np
import pytest
import yaml
from pathlib import Path
from tifffile import imwrite

from run_pipeline import (
    detect_image_extension,
    discover_samples,
    load_config,
    sample_is_complete,
)


def _create_sample_directory(base: Path, name: str, extension: str = ".tif") -> Path:
    data_dir = base / name / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    image_path = data_dir / f"{name}{extension}"
    imwrite(str(image_path), np.zeros((3, 64, 64), dtype=np.uint16))
    return image_path


def _write_config(path: Path, overrides: dict | None = None) -> Path:
    config = {
        "image_directory": str(path),
        "nuclear_channel": "DAPI",
        "segment_channel": "E_CADHERIN",
        "segmentation_method": "cellpose",
    }
    if overrides:
        config.update(overrides)

    config_path = path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


def test_load_config_valid(tmp_path):
    _create_sample_directory(tmp_path, "sample_001")
    config_path = _write_config(tmp_path)

    config = load_config(config_path)

    assert config["image_directory"] == tmp_path
    assert config["nuclear_channel"] == "DAPI"
    assert config["segment_channel"] == ["E_CADHERIN"]
    assert config["segmentation_method"] == "cellpose"


def test_load_config_comma_separated_channels(tmp_path):
    _create_sample_directory(tmp_path, "sample_001")
    config_path = _write_config(
        tmp_path, {"segment_channel": "PDGFRB, CD3, CYTOKERATIN"}
    )

    config = load_config(config_path)

    assert config["segment_channel"] == ["PDGFRB", "CD3", "CYTOKERATIN"]


def test_load_config_missing_key(tmp_path):
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump({"image_directory": str(tmp_path)}, f)

    with pytest.raises(ValueError, match="missing required config keys"):
        load_config(config_path)


def test_load_config_invalid_method(tmp_path):
    _create_sample_directory(tmp_path, "sample_001")
    config_path = _write_config(tmp_path, {"segmentation_method": "invalid"})

    with pytest.raises(ValueError, match="segmentation_method must be"):
        load_config(config_path)


def test_detect_image_extension(tmp_path):
    _create_sample_directory(tmp_path, "sample_001", ".tif")

    assert detect_image_extension(tmp_path) == ".tif"


def test_detect_image_extension_qptiff(tmp_path):
    _create_sample_directory(tmp_path, "sample_001", ".qptiff")

    assert detect_image_extension(tmp_path) == ".qptiff"


def test_detect_image_extension_no_images(tmp_path):
    (tmp_path / "empty_dir").mkdir()

    with pytest.raises(FileNotFoundError, match="no image files"):
        detect_image_extension(tmp_path)


def test_discover_samples(tmp_path):
    _create_sample_directory(tmp_path, "sample_001")
    _create_sample_directory(tmp_path, "sample_002")
    (tmp_path / "histology").mkdir()
    (tmp_path / "quantifications").mkdir()

    samples = discover_samples(tmp_path, ".tif")

    assert len(samples) == 2
    names = [s.stem for s in samples]
    assert "sample_001" in names
    assert "sample_002" in names


def test_sample_is_complete_false_when_no_outputs(tmp_path):
    image_path = _create_sample_directory(tmp_path, "sample_001")

    assert not sample_is_complete(image_path, "cellpose", ["E_CADHERIN"])


def test_sample_is_complete_true_when_outputs_exist(tmp_path):
    image_path = _create_sample_directory(tmp_path, "sample_001")

    mask_dir = tmp_path / "sample_001" / "full" / "cellpose" / "E_CADHERIN"
    mask_dir.mkdir(parents=True)
    np.save(str(mask_dir / "image_1_seg.npy"), np.zeros((64, 64)))

    csv_dir = tmp_path / "quantifications" / "cellpose"
    csv_dir.mkdir(parents=True)
    (csv_dir / "sample_001_cell_measurements.csv").write_text("header\n")

    assert sample_is_complete(image_path, "cellpose", ["E_CADHERIN"])
