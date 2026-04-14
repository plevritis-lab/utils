from pathlib import Path

import pandas as pd

EXCLUDED_DIRECTORIES = {
    "assignments",
    "quantifications",
    "thresholds",
}


def process_thresholds(
    image_directory: Path,
    signature_matrix_path: Path,
    save_path: Path,
    sample_filter: list[str] | None = None,
) -> None:
    """generates or updates per-sample threshold CSVs from a signature matrix"""

    signature_matrix = pd.read_csv(signature_matrix_path)
    save_path.mkdir(parents=True, exist_ok=True)

    sample_directories = [
        d.name
        for d in sorted(image_directory.iterdir())
        if d.is_dir() and d.name not in EXCLUDED_DIRECTORIES
    ]

    if sample_filter is not None:
        sample_directories = [s for s in sample_directories if s in sample_filter]

    existing_threshold_files = {
        f.stem.replace("_thresholds", ""): f
        for f in save_path.iterdir()
        if f.name.endswith("_thresholds.csv")
    }

    for sample_name in sample_directories:
        threshold_path = save_path / f"{sample_name}_thresholds.csv"

        if sample_name in existing_threshold_files:
            thresholds = pd.read_csv(existing_threshold_files[sample_name])

            rows = []
            for cell_type in signature_matrix["CELL_TYPE"]:
                if cell_type in thresholds["CELL_TYPE"].values:
                    row = thresholds[thresholds["CELL_TYPE"] == cell_type].iloc[0]
                    rows.append(
                        {
                            "CELL_TYPE": cell_type,
                            "ANCHOR": row["ANCHOR"],
                            "ITERATION": row["ITERATION"],
                        }
                    )
                else:
                    rows.append(
                        {"CELL_TYPE": cell_type, "ANCHOR": 0.7, "ITERATION": 0.5}
                    )

            updated_thresholds = pd.DataFrame(rows)
            updated_thresholds.to_csv(threshold_path, index=False)
        else:
            thresholds = pd.DataFrame(
                {
                    "CELL_TYPE": signature_matrix["CELL_TYPE"],
                    "ANCHOR": [0.7] * len(signature_matrix),
                    "ITERATION": [0.5] * len(signature_matrix),
                }
            )
            thresholds.to_csv(threshold_path, index=False)
