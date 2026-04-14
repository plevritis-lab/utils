import numpy as np
import yaml
from pathlib import Path
from tifffile import imwrite


def _create_sample_directory(base: Path, name: str, extension: str = ".tif") -> Path:
    data_dir = base / name / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    image_path = data_dir / f"{name}{extension}"
    imwrite(str(image_path), np.zeros((3, 64, 64), dtype=np.uint16))
    return image_path


def _write_panel(path: Path) -> Path:
    panel_path = path / "channel_names.txt"
    panel_path.write_text("DAPI\nE_CADHERIN\nCD3\n")
    return panel_path


def _write_config(path: Path, overrides: dict | None = None) -> Path:
    panel_path = _write_panel(path)
    config = {
        "image_directory": str(path),
        "panel_path": str(panel_path),
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


def test_load_config_and_discover_samples(tmp_path):
    _create_sample_directory(tmp_path, "sample_001")
    _create_sample_directory(tmp_path, "sample_002")
    (tmp_path / "histology").mkdir()
    config_path = _write_config(tmp_path)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    assert config["nuclear_channel"] == "DAPI"
    assert config["segmentation_method"] == "cellpose"

    image_directory = Path(config["image_directory"])
    sample_dirs = [
        d
        for d in sorted(image_directory.iterdir())
        if d.is_dir() and d.name not in {"histology", "quantifications", "assignments"}
    ]
    assert len(sample_dirs) == 2


def test_skip_when_outputs_exist(tmp_path, segmentation_pipeline):
    _create_sample_directory(tmp_path, "sample_001")

    mask_dir = tmp_path / "sample_001" / "full" / "cellpose" / "E_CADHERIN"
    mask_dir.mkdir(parents=True)
    np.save(str(mask_dir / "image_1_seg.npy"), np.zeros((64, 64)))

    csv_dir = tmp_path / "quantifications" / "cellpose"
    csv_dir.mkdir(parents=True)
    (csv_dir / "sample_001_cell_measurements.csv").write_text("header\n")

    config_path = _write_config(tmp_path)

    segmentation_pipeline.run_pipeline(config_path, force=False)

    assert (tmp_path / "segmentation_pipeline.log").is_file()
    log_content = (tmp_path / "segmentation_pipeline.log").read_text()
    assert "[skip] sample_001" in log_content
