import pandas as pd

from generate_thresholds import process_thresholds


def test_seeds_default_thresholds_for_new_samples(tmp_path):
    image_directory = tmp_path / "images"
    (image_directory / "sample_a").mkdir(parents=True)
    (image_directory / "sample_b").mkdir(parents=True)

    signature_matrix_path = tmp_path / "signature_matrix.csv"
    pd.DataFrame(
        {"CELL_TYPE": ["epithelial", "immune"], "Lineage_level": ["1_0_1", "1_0_2"]}
    ).to_csv(signature_matrix_path, index=False)

    save_path = tmp_path / "thresholds"

    process_thresholds(image_directory, signature_matrix_path, save_path)

    for sample in ["sample_a", "sample_b"]:
        threshold_file = save_path / f"{sample}_thresholds.csv"
        assert threshold_file.exists()

        thresholds = pd.read_csv(threshold_file)
        assert list(thresholds.columns) == ["CELL_TYPE", "ANCHOR", "ITERATION"]
        assert list(thresholds["CELL_TYPE"]) == ["epithelial", "immune"]
        assert all(thresholds["ANCHOR"] == 0.7)
        assert all(thresholds["ITERATION"] == 0.5)


def test_preserves_existing_thresholds_on_rerun(tmp_path):
    image_directory = tmp_path / "images"
    (image_directory / "sample_a").mkdir(parents=True)

    signature_matrix_path = tmp_path / "signature_matrix.csv"
    pd.DataFrame(
        {"CELL_TYPE": ["epithelial", "immune"], "Lineage_level": ["1_0_1", "1_0_2"]}
    ).to_csv(signature_matrix_path, index=False)

    save_path = tmp_path / "thresholds"
    save_path.mkdir(parents=True)

    pd.DataFrame(
        {
            "CELL_TYPE": ["epithelial", "immune"],
            "ANCHOR": [0.9, 0.6],
            "ITERATION": [0.8, 0.3],
        }
    ).to_csv(save_path / "sample_a_thresholds.csv", index=False)

    process_thresholds(image_directory, signature_matrix_path, save_path)

    thresholds = pd.read_csv(save_path / "sample_a_thresholds.csv")
    assert list(thresholds["ANCHOR"]) == [0.9, 0.6]
    assert list(thresholds["ITERATION"]) == [0.8, 0.3]


def test_filter_limits_threshold_generation(tmp_path):
    image_directory = tmp_path / "images"
    (image_directory / "sample_a").mkdir(parents=True)
    (image_directory / "sample_b").mkdir(parents=True)

    signature_matrix_path = tmp_path / "signature_matrix.csv"
    pd.DataFrame({"CELL_TYPE": ["epithelial"], "Lineage_level": ["1_0_1"]}).to_csv(
        signature_matrix_path, index=False
    )

    save_path = tmp_path / "thresholds"

    process_thresholds(
        image_directory,
        signature_matrix_path,
        save_path,
        sample_filter=["sample_a"],
    )

    assert (save_path / "sample_a_thresholds.csv").exists()
    assert not (save_path / "sample_b_thresholds.csv").exists()
