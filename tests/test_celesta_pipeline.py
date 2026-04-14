from unittest.mock import patch

import pandas as pd
import yaml


def test_colormap_auto_generation(tmp_path, celesta_pipeline):
    signature_matrix_path = tmp_path / "signature_matrix.csv"
    pd.DataFrame(
        {
            "CELL_TYPE": ["epithelial", "immune", "stromal"],
            "Lineage_level": ["1_0_1", "1_0_2", "1_0_3"],
        }
    ).to_csv(signature_matrix_path, index=False)

    colormap_path = tmp_path / "colormap.json"
    colormap = celesta_pipeline.generate_colormap(signature_matrix_path, colormap_path)

    assert colormap_path.exists()
    assert "epithelial" in colormap
    assert "immune" in colormap
    assert "stromal" in colormap
    assert "unknown" in colormap
    assert colormap["unknown"]["color"] == "#7F7F7F"

    for cell_type in ["epithelial", "immune", "stromal"]:
        assert "color" in colormap[cell_type]
        assert colormap[cell_type]["color"].startswith("#")


def test_config_loading_and_filter_passthrough(tmp_path, celesta_pipeline):
    image_directory = tmp_path / "images"
    (image_directory / "sample_a").mkdir(parents=True)
    (image_directory / "quantifications" / "cellpose").mkdir(parents=True)

    signature_matrix_path = tmp_path / "signature_matrix.csv"
    pd.DataFrame({"CELL_TYPE": ["epithelial"], "Lineage_level": ["1_0_1"]}).to_csv(
        signature_matrix_path, index=False
    )

    config_path = tmp_path / "config.yaml"
    config = {
        "image_directory": str(image_directory),
        "signature_matrix": str(signature_matrix_path),
        "segmentation_method": "cellpose",
        "segment_channel": "CYTOKERATIN",
    }
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    captured_commands = []

    def mock_subprocess_run(command, **kwargs):
        captured_commands.append(command)

        class MockResult:
            returncode = 0
            stdout = ""
            stderr = ""

        return MockResult()

    with patch.object(
        celesta_pipeline.subprocess, "run", side_effect=mock_subprocess_run
    ):
        celesta_pipeline.run_pipeline(config_path, sample_filter=["sample_a"])

    assert len(captured_commands) == 1
    command = captured_commands[0]
    assert "Rscript" in command[0]
    assert "--filter" in command
    filter_index = command.index("--filter")
    assert command[filter_index + 1] == "sample_a"
