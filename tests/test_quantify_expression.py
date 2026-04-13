import numpy as np
import pandas as pd
from pathlib import Path

from quantify_expression import quantify_expression


def test_quantify_expression_produces_expected_columns(tmp_path):
    image = np.random.randint(0, 255, (3, 64, 64), dtype=np.uint16)

    mask = np.zeros((64, 64), dtype=np.int32)
    mask[10:20, 10:20] = 1
    mask[30:40, 30:40] = 2

    panel = ["DAPI", "E_CADHERIN", "CD3"]
    save_path = str(tmp_path / "test_sample")

    quantify_expression(image, mask, panel, save_path)

    csv_path = Path(f"{save_path}_cell_measurements.csv")
    assert csv_path.is_file()

    df = pd.read_csv(csv_path)

    expected_columns = {
        "CELL_IDENTIFIER",
        "X",
        "Y",
        "SIZE",
        "MAJOR_AXIS_LENGTH",
        "MINOR_AXIS_LENGTH",
        "ECCENTRICITY",
        "ORIENTATION",
        "DAPI",
        "E_CADHERIN",
        "CD3",
    }
    assert set(df.columns) == expected_columns
    assert len(df) == 2


def test_quantify_expression_cell_positions_are_reasonable(tmp_path):
    image = np.random.randint(0, 255, (2, 100, 100), dtype=np.uint16)

    mask = np.zeros((100, 100), dtype=np.int32)
    mask[40:60, 40:60] = 1

    panel = ["MARKER_A", "MARKER_B"]
    save_path = str(tmp_path / "position_test")

    quantify_expression(image, mask, panel, save_path)

    df = pd.read_csv(f"{save_path}_cell_measurements.csv")

    assert 40 <= df.iloc[0]["X"] <= 60
    assert 40 <= df.iloc[0]["Y"] <= 60
