import base64
import io
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tifffile import imread

matplotlib.use("Agg")

PROBABILITY_BINS = ["<=0.5", ">0.5", ">0.7", ">0.8", ">0.9"]
PROBABILITY_COLORS = sns.color_palette("light:b", 5)
POINT_SIZE = 6


def _figure_to_base64(fig: plt.Figure) -> str:
    """converts a matplotlib figure to a base64-encoded PNG string"""

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=150, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    buffer.seek(0)
    return base64.b64encode(buffer.read()).decode("utf-8")


def _bin_probabilities(values: pd.Series) -> list[str]:
    """bins probability values into discrete categories"""

    bins = []
    for value in values:
        if value > 0.9:
            bins.append(">0.9")
        elif value > 0.8:
            bins.append(">0.8")
        elif value > 0.7:
            bins.append(">0.7")
        elif value > 0.5:
            bins.append(">0.5")
        else:
            bins.append("<=0.5")
    return bins


def _render_overlay(
    raw_channel: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    probability_bins: list[str],
) -> str:
    """renders probability points overlaid on a raw fluorescence channel"""

    fig, ax = plt.subplots(figsize=(8, 8), facecolor="black")
    ax.set_facecolor("black")

    ax.imshow(raw_channel, cmap="gray", aspect="equal")

    bin_to_color = dict(zip(PROBABILITY_BINS, PROBABILITY_COLORS))

    for probability_bin in PROBABILITY_BINS:
        mask = [b == probability_bin for b in probability_bins]
        if not any(mask):
            continue
        ax.scatter(
            x[mask],
            y[mask],
            c=[bin_to_color[probability_bin]],
            s=POINT_SIZE,
            alpha=0.7,
            label=probability_bin,
        )

    ax.set_xlim(0, raw_channel.shape[1])
    ax.set_ylim(raw_channel.shape[0], 0)
    ax.set_axis_off()

    legend = ax.legend(
        loc="upper right",
        facecolor="black",
        edgecolor="white",
        labelcolor="white",
        fontsize=8,
        markerscale=2,
    )
    legend.set_alpha(0.8)

    return _figure_to_base64(fig)


def _render_raw_channel(raw_channel: np.ndarray) -> str:
    """renders a raw fluorescence channel as grayscale"""

    fig, ax = plt.subplots(figsize=(8, 8), facecolor="black")
    ax.set_facecolor("black")
    ax.imshow(raw_channel, cmap="gray", aspect="equal")
    ax.set_axis_off()
    return _figure_to_base64(fig)


def _render_assignment_scatter(
    assignments: pd.DataFrame,
    colormap: dict,
) -> str:
    """renders a spatial scatter plot of cell type assignments"""

    fig, ax = plt.subplots(figsize=(8, 8), facecolor="black")
    ax.set_facecolor("black")
    ax.invert_yaxis()
    ax.set_axis_off()

    for cell_type in colormap:
        subset = assignments[assignments["FINAL_CELL_TYPE"] == cell_type]
        if subset.empty:
            continue
        ax.scatter(
            subset["X"],
            subset["Y"],
            c=colormap[cell_type]["color"],
            label=cell_type,
            s=POINT_SIZE,
        )

    ax.legend(
        loc="upper right",
        facecolor="black",
        edgecolor="white",
        labelcolor="white",
        fontsize=7,
        markerscale=2,
    )

    return _figure_to_base64(fig)


def _render_proportions(
    assignments: pd.DataFrame,
    colormap: dict,
) -> str:
    """renders a cell proportions bar chart"""

    fig, ax = plt.subplots(figsize=(4, 6), facecolor="white")

    counts = assignments["FINAL_CELL_TYPE"].value_counts()
    proportions = (counts / counts.sum()) * 100

    cell_types = [ct for ct in colormap if ct in proportions.index]
    values = [proportions.get(ct, 0) for ct in cell_types]
    colors = [colormap[ct]["color"] for ct in cell_types]

    ax.barh(cell_types, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("percentage (%)")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)

    return _figure_to_base64(fig)


def _expression_label(value: float) -> str:
    """converts a signature matrix value to a human-readable label"""

    if pd.isna(value):
        return "N/A"
    elif value == 1:
        return "positive"
    elif value == 0:
        return "negative"
    return str(value)


def _expression_badge_class(value: float) -> str:
    """returns a CSS class name for the expression badge"""

    if pd.isna(value):
        return "badge-na"
    elif value == 1:
        return "badge-positive"
    elif value == 0:
        return "badge-negative"
    return "badge-na"


def generate_threshold_report(
    assignments_path: Path,
    image_path: Path,
    panel: list[str],
    signature_matrix: pd.DataFrame,
    thresholds: pd.DataFrame,
    colormap: dict,
    sample_name: str,
    save_path: Path,
) -> Path:
    """generates a self-contained HTML threshold tuning report for a single sample"""

    assignments = pd.read_csv(assignments_path).dropna()
    image_data = imread(str(image_path))

    x = assignments["X"].values
    y = assignments["Y"].values

    cell_types = signature_matrix["CELL_TYPE"].tolist()
    sig_markers = [
        col
        for col in signature_matrix.columns
        if col not in ("CELL_TYPE", "Lineage_level")
    ]

    assignment_scatter_b64 = _render_assignment_scatter(assignments, colormap)
    proportions_b64 = _render_proportions(assignments, colormap)

    cell_type_sections = []

    for cell_type in cell_types:
        row = signature_matrix[signature_matrix["CELL_TYPE"] == cell_type].iloc[0]
        threshold_row = thresholds[thresholds["CELL_TYPE"] == cell_type]

        anchor_val = (
            threshold_row["ANCHOR"].values[0] if not threshold_row.empty else "N/A"
        )
        iteration_val = (
            threshold_row["ITERATION"].values[0] if not threshold_row.empty else "N/A"
        )

        marker_sections = []

        for marker in sig_markers:
            expected_value = row[marker]
            label = _expression_label(expected_value)
            badge_class = _expression_badge_class(expected_value)

            probability_column = f"{marker}_PROBABILITY"

            if probability_column in assignments.columns and marker in panel:
                channel_index = panel.index(marker)
                raw_channel = image_data[channel_index]

                probabilities = assignments[probability_column].values
                bins = _bin_probabilities(pd.Series(probabilities))

                overlay_b64 = _render_overlay(raw_channel, x, y, bins)
                raw_b64 = _render_raw_channel(raw_channel)

                marker_sections.append(
                    f"""
                    <div class="marker-section">
                        <h3>{marker} <span class="badge {badge_class}">expected: {label}</span></h3>
                        <div class="image-pair">
                            <div class="image-container">
                                <p class="image-label">probability overlay</p>
                                <img src="data:image/png;base64,{overlay_b64}" />
                            </div>
                            <div class="image-container">
                                <p class="image-label">raw channel</p>
                                <img src="data:image/png;base64,{raw_b64}" />
                            </div>
                        </div>
                    </div>
                    """
                )
            else:
                marker_sections.append(
                    f"""
                    <div class="marker-section">
                        <h3>{marker} <span class="badge {badge_class}">expected: {label}</span></h3>
                        <p class="no-data">no probability data or channel not in panel</p>
                    </div>
                    """
                )

        cell_type_sections.append(
            {
                "name": cell_type,
                "anchor": anchor_val,
                "iteration": iteration_val,
                "markers_html": "\n".join(marker_sections),
            }
        )

    sidebar_links = "\n".join(
        f'<a href="#ct-{i}">{section["name"]}</a>'
        for i, section in enumerate(cell_type_sections)
    )

    body_sections = "\n".join(
        f"""
        <section id="ct-{i}" class="cell-type-section">
            <h2>{section["name"]}</h2>
            <div class="threshold-info">
                <span class="threshold">ANCHOR: {section["anchor"]}</span>
                <span class="threshold">INDEX: {section["iteration"]}</span>
            </div>
            {section["markers_html"]}
        </section>
        """
        for i, section in enumerate(cell_type_sections)
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>threshold report — {sample_name}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #1a1a2e; color: #e0e0e0; display: flex; }}

    .sidebar {{
        position: fixed; top: 0; left: 0; width: 240px; height: 100vh;
        background: #16213e; overflow-y: auto; padding: 20px 0;
        border-right: 1px solid #0f3460;
    }}
    .sidebar h1 {{ font-size: 14px; padding: 0 16px 12px; color: #a0c4ff; border-bottom: 1px solid #0f3460; margin-bottom: 8px; }}
    .sidebar a {{
        display: block; padding: 8px 16px; color: #c0c0c0; text-decoration: none;
        font-size: 12px; border-left: 3px solid transparent; transition: all 0.15s;
    }}
    .sidebar a:hover {{ background: #0f3460; color: #fff; border-left-color: #a0c4ff; }}

    .sidebar .nav-section {{ padding: 12px 16px 4px; font-size: 11px; color: #607090; text-transform: uppercase; letter-spacing: 1px; }}

    .main {{ margin-left: 240px; padding: 32px; width: calc(100% - 240px); }}

    .header {{ margin-bottom: 32px; }}
    .header h1 {{ font-size: 24px; color: #a0c4ff; }}
    .header p {{ color: #808080; font-size: 13px; margin-top: 4px; }}

    .overview {{ display: flex; gap: 24px; margin-bottom: 40px; flex-wrap: wrap; }}
    .overview .image-container {{ background: #16213e; border-radius: 8px; padding: 12px; }}
    .overview .image-container img {{ max-width: 500px; border-radius: 4px; }}

    .cell-type-section {{
        margin-bottom: 48px; padding-bottom: 32px;
        border-bottom: 1px solid #0f3460;
    }}
    .cell-type-section h2 {{ font-size: 20px; color: #a0c4ff; margin-bottom: 8px; }}

    .threshold-info {{ margin-bottom: 20px; }}
    .threshold {{
        display: inline-block; background: #0f3460; padding: 4px 12px;
        border-radius: 4px; font-size: 13px; margin-right: 8px; font-family: monospace;
    }}

    .marker-section {{ margin-bottom: 24px; }}
    .marker-section h3 {{ font-size: 15px; margin-bottom: 8px; color: #d0d0d0; }}

    .badge {{
        display: inline-block; padding: 2px 8px; border-radius: 3px;
        font-size: 11px; font-weight: normal; margin-left: 8px;
    }}
    .badge-positive {{ background: #1b4332; color: #95d5b2; }}
    .badge-negative {{ background: #3c1518; color: #e5989b; }}
    .badge-na {{ background: #343a40; color: #adb5bd; }}

    .image-pair {{ display: flex; gap: 16px; flex-wrap: wrap; }}
    .image-container {{ background: #16213e; border-radius: 8px; padding: 8px; }}
    .image-container img {{ max-width: 450px; width: 100%; border-radius: 4px; }}
    .image-label {{ font-size: 11px; color: #808080; margin-bottom: 4px; text-transform: uppercase; letter-spacing: 0.5px; }}

    .no-data {{ color: #606060; font-style: italic; font-size: 13px; padding: 8px 0; }}
</style>
</head>
<body>

<nav class="sidebar">
    <h1>{sample_name}</h1>
    <div class="nav-section">overview</div>
    <a href="#overview">assignment overview</a>
    <div class="nav-section">cell types</div>
    {sidebar_links}
</nav>

<div class="main">
    <div class="header">
        <h1>threshold tuning report</h1>
        <p>{sample_name}</p>
    </div>

    <section id="overview">
        <h2>assignment overview</h2>
        <div class="overview">
            <div class="image-container">
                <p class="image-label">spatial assignments</p>
                <img src="data:image/png;base64,{assignment_scatter_b64}" />
            </div>
            <div class="image-container">
                <p class="image-label">cell proportions</p>
                <img src="data:image/png;base64,{proportions_b64}" />
            </div>
        </div>
    </section>

    {body_sections}
</div>

</body>
</html>"""

    output_path = save_path / f"{sample_name}_threshold_report.html"
    output_path.write_text(html)

    return output_path
