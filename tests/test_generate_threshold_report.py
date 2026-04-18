import numpy as np
import pandas as pd
from tifffile import imwrite

from generate_threshold_report import generate_threshold_report


def _create_test_data(tmp_path):
    assignments = pd.DataFrame(
        {
            "CELL_IDENTIFIER": [1, 2, 3, 4, 5],
            "X": [100.0, 200.0, 300.0, 400.0, 500.0],
            "Y": [100.0, 200.0, 300.0, 400.0, 500.0],
            "CD3_PROBABILITY": [0.95, 0.85, 0.3, 0.1, 0.6],
            "CYTOKERATIN_PROBABILITY": [0.1, 0.05, 0.92, 0.88, 0.2],
            "CELL_TYPE_NUMBER": [1, 1, 2, 2, 0],
            "FINAL_CELL_TYPE": [
                "T cell",
                "T cell",
                "epithelial",
                "epithelial",
                "unknown",
            ],
        }
    )
    assignments_path = tmp_path / "sample_a_assignments.csv"
    assignments.to_csv(assignments_path, index=False)

    image = np.random.randint(0, 255, (3, 64, 64), dtype=np.uint16)
    image_path = tmp_path / "sample_a.tif"
    imwrite(str(image_path), image)

    panel = ["DRAQ5", "CD3", "CYTOKERATIN"]

    signature_matrix = pd.DataFrame(
        {
            "CELL_TYPE": ["T cell", "epithelial"],
            "Lineage_level": ["1_0_1", "1_0_2"],
            "CD3": [1, 0],
            "CYTOKERATIN": [0, 1],
        }
    )

    thresholds = pd.DataFrame(
        {
            "CELL_TYPE": ["T cell", "epithelial"],
            "ANCHOR": [0.8, 0.7],
            "ITERATION": [0.5, 0.5],
        }
    )

    colormap = {
        "T cell": {"color": "#FF0000"},
        "epithelial": {"color": "#00FF00"},
        "unknown": {"color": "#7F7F7F"},
    }

    return assignments_path, image_path, panel, signature_matrix, thresholds, colormap


def test_report_generates_html_file(tmp_path):
    (
        assignments_path,
        image_path,
        panel,
        signature_matrix,
        thresholds,
        colormap,
    ) = _create_test_data(tmp_path)

    output_path = generate_threshold_report(
        assignments_path=assignments_path,
        image_path=image_path,
        panel=panel,
        signature_matrix=signature_matrix,
        thresholds=thresholds,
        colormap=colormap,
        sample_name="sample_a",
        save_path=tmp_path,
    )

    assert output_path.exists()
    assert output_path.suffix == ".html"
    assert output_path.stat().st_size > 0


def test_report_contains_expected_content(tmp_path):
    (
        assignments_path,
        image_path,
        panel,
        signature_matrix,
        thresholds,
        colormap,
    ) = _create_test_data(tmp_path)

    output_path = generate_threshold_report(
        assignments_path=assignments_path,
        image_path=image_path,
        panel=panel,
        signature_matrix=signature_matrix,
        thresholds=thresholds,
        colormap=colormap,
        sample_name="sample_a",
        save_path=tmp_path,
    )

    html = output_path.read_text()

    assert "sample_a" in html
    assert "T cell" in html
    assert "epithelial" in html
    assert "ANCHOR: 0.8" in html
    assert "ANCHOR: 0.7" in html
    assert "expected: positive" in html
    assert "expected: negative" in html
    assert "data:image/png;base64," in html
