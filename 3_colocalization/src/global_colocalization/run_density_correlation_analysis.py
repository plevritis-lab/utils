import argparse
import cv2
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from scipy.stats import pearsonr, ttest_ind
import seaborn as sns
from skimage.io import imread

# TODO: normalize values!


def identify_metastatic_lymph_nodes():
    clinical_data = pd.read_csv("../../2_celesta/input/data/clinical_data.csv")

    lymph_node_data = clinical_data[clinical_data["Site"].str.contains("LN", na=False)]
    lymph_node_positive = lymph_node_data[
        ~lymph_node_data["N"].str.contains("N0", na=False)
    ]

    microarrays = (
        lymph_node_positive["TMA#"].astype(str)
        + "_"
        + lymph_node_positive["TMA_part"].astype(str)
    ).tolist()

    filtered_cores = lymph_node_positive[lymph_node_positive["Core_image_ID"] != -99]
    cores = filtered_cores["Core_image_ID"].apply(lambda x: f"reg{int(x):03}").tolist()

    return list(zip(microarrays, cores))


def compute_stain_correlation(
    image_path,
    mode,
    microarray,
    core,
    condition,
    marker_one,
    marker_two,
    tile_size,
    marker_one_threshold,
    marker_two_threshold,
):
    image = imread(image_path)  # (30, 1440, 1920, 4)

    panel = []
    with open("../../2_celesta/input/data/channel_names.txt") as file:
        for line in file:
            panel.append(line.rstrip())

    assert image.shape[0] * image.shape[-1] == len(panel)

    marker_one_index = panel.index(marker_one)
    marker_two_index = panel.index(marker_two)

    marker_one_stain = image[marker_one_index // 4, :, :, marker_one_index % 4]
    marker_two_stain = image[marker_two_index // 4, :, :, marker_two_index % 4]

    binary_marker_one_stain = (
        marker_one_stain >= (marker_one_threshold * np.max(marker_one_stain))
    ).astype(np.uint8)
    binary_marker_two_stain = (
        marker_two_stain >= (marker_two_threshold * np.max(marker_two_stain))
    ).astype(np.uint8)

    rows, columns = binary_marker_one_stain.shape
    num_tiles_y = rows // tile_size + (rows % tile_size != 0)
    num_tiles_x = columns // tile_size + (columns % tile_size != 0)

    rgb_marker_one_stain = cv2.cvtColor(
        (binary_marker_one_stain * 255).astype(np.uint8), cv2.COLOR_GRAY2BGR
    )
    rgb_marker_two_stain = cv2.cvtColor(
        (binary_marker_two_stain * 255).astype(np.uint8), cv2.COLOR_GRAY2BGR
    )

    for i in range(1, num_tiles_x):
        x = i * tile_size
        cv2.line(rgb_marker_one_stain, (x, 0), (x, rows), (0, 0, 255), 1)
        cv2.line(rgb_marker_two_stain, (x, 0), (x, rows), (0, 0, 255), 1)

    for i in range(1, num_tiles_y):
        y = i * tile_size
        cv2.line(rgb_marker_one_stain, (0, y), (columns, y), (0, 0, 255), 1)
        cv2.line(rgb_marker_two_stain, (0, y), (columns, y), (0, 0, 255), 1)

    os.makedirs(f"../output/{mode}/{condition}/{microarray}_{core}", exist_ok=True)
    cv2.imwrite(
        f"../output/{mode}/{condition}/{microarray}_{core}/{marker_one}_mask_threshold_{marker_one_threshold}.png",
        rgb_marker_one_stain,
    )
    cv2.imwrite(
        f"../output/{mode}/{condition}/{microarray}_{core}/{marker_two}_mask_threshold_{marker_two_threshold}.png",
        rgb_marker_two_stain,
    )

    marker_one_density_map = np.zeros((num_tiles_y, num_tiles_x))
    marker_two_density_map = np.zeros((num_tiles_y, num_tiles_x))

    for i in range(0, rows, tile_size):
        for j in range(0, columns, tile_size):
            marker_one_tile = binary_marker_one_stain[
                i : min(i + tile_size, rows), j : min(j + tile_size, columns)
            ]
            marker_two_tile = binary_marker_two_stain[
                i : min(i + tile_size, rows), j : min(j + tile_size, columns)
            ]

            marker_one_proportion = np.sum(marker_one_tile) / (
                marker_one_tile.shape[0] * marker_one_tile.shape[1]
            )
            marker_two_proportion = np.sum(marker_two_tile) / (
                marker_two_tile.shape[0] * marker_two_tile.shape[1]
            )

            marker_one_density_map[i // tile_size, j // tile_size] = (
                marker_one_proportion
            )
            marker_two_density_map[i // tile_size, j // tile_size] = (
                marker_two_proportion
            )

    fig, (ax_one, ax_two) = plt.subplots(1, 2, figsize=(10, 6))

    im = ax_one.imshow(marker_one_density_map, cmap="hot", interpolation="nearest")
    fig.colorbar(im)
    ax_one.set_title(f"density map of {marker_one}")

    im = ax_two.imshow(marker_two_density_map, cmap="hot", interpolation="nearest")
    fig.colorbar(im)
    ax_two.set_title(f"density map of {marker_two}")

    fig.savefig(
        f"../output/{mode}/{condition}/{microarray}_{core}/{marker_one}_{marker_two}_heatmap.pdf"
    )
    plt.close(fig)

    flattened_marker_one_density_map = marker_one_density_map.flatten()
    flattened_marker_two_density_map = marker_two_density_map.flatten()

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(
        flattened_marker_one_density_map, flattened_marker_two_density_map, alpha=0.5
    )
    ax.set_xlabel(f"{marker_one} tile proportions")
    ax.set_ylabel(f"{marker_two} tile proportions")

    polynomial = np.poly1d(
        np.polyfit(
            flattened_marker_one_density_map, flattened_marker_two_density_map, 1
        )
    )
    line_of_best_fit = polynomial(flattened_marker_one_density_map)

    ax.plot(flattened_marker_one_density_map, line_of_best_fit, color="red")

    correlation_coefficient, _ = pearsonr(
        flattened_marker_one_density_map, flattened_marker_two_density_map
    )

    ax.text(
        x=0.05,
        y=0.95,
        s=f"pearson $r$: {correlation_coefficient:.2f}",
        fontsize=12,
        verticalalignment="top",
        horizontalalignment="left",
        transform=plt.gca().transAxes,
    )

    fig.savefig(
        f"../output/{mode}/{condition}/{microarray}_{core}/{marker_one}_{marker_two}_correlation_scatter.pdf"
    )
    plt.close(fig)

    return correlation_coefficient


def parse_args():
    parser = argparse.ArgumentParser(
        description="interface for selecting comparator markers"
    )

    parser.add_argument("-m", "--mode")
    parser.add_argument("-m1", "--marker_one")
    parser.add_argument("-m2", "--marker_two")
    parser.add_argument("-t", "--tile_size", type=int)
    parser.add_argument("-t1", "--threshold_one", type=float)
    parser.add_argument("-t2", "--threshold_two", type=float)

    return parser.parse_args()


def main():
    args = parse_args()

    mode = args.mode
    marker_one = args.marker_one
    marker_two = args.marker_two
    tile_size = args.tile_size
    marker_one_threshold = args.threshold_one
    marker_two_threshold = args.threshold_two

    negative_assignments_directory = f"../input/{mode}/negative"
    positive_assignments_directory = f"../input/{mode}/positive"
    # lymph_node_samples = identify_metastatic_lymph_nodes()

    negative_correlations = []
    positive_correlations = []
    # lymph_node_correlations = []

    for file in os.listdir(negative_assignments_directory):
        match = re.match(r"^(.*?)_(reg\d+)_.*$", file)

        if match:
            microarray = match.group(1)
            core = match.group(2)

            image_path = f"../../2_celesta/input/data/{microarray}/processed_imaging_data/{core}_X01_Y01_Z01.tif"

            negative_correlations.append(
                compute_stain_correlation(
                    image_path,
                    mode,
                    microarray,
                    core,
                    "node_negative_stain",
                    marker_one,
                    marker_two,
                    tile_size,
                    marker_one_threshold,
                    marker_two_threshold,
                )
            )

    for file in os.listdir(positive_assignments_directory):
        match = re.match(r"^(.*?)_(reg\d+)_.*$", file)

        if match:
            microarray = match.group(1)
            core = match.group(2)

            image_path = f"../../2_celesta/input/data/{microarray}/processed_imaging_data/{core}_X01_Y01_Z01.tif"

            positive_correlations.append(
                compute_stain_correlation(
                    image_path,
                    mode,
                    microarray,
                    core,
                    "node_positive_stain",
                    marker_one,
                    marker_two,
                    tile_size,
                    marker_one_threshold,
                    marker_two_threshold,
                )
            )

    # for microarray, core in lymph_node_samples:
    #     image_path = f"../../2_celesta/input/data/{microarray}/processed_imaging_data/{core}_X01_Y01_Z01.tif"

    #     lymph_node_correlations.append(compute_stain_correlation(image_path, mode, microarray, core, "lymph_node_stain", marker_one, marker_two,
    #                                                                  tile_size, marker_one_threshold, marker_two_threshold))

    _, p_value = ttest_ind(positive_correlations, negative_correlations)
    # _, p_value = ttest_ind(positive_correlations, lymph_node_correlations)
    # _, p_value = ttest_ind(negative_correlations, lymph_node_correlations)

    data = pd.DataFrame(
        {
            "density correlation": positive_correlations + negative_correlations,
            "category": ["N+"] * len(positive_correlations)
            + ["N0"] * len(negative_correlations),
        }
    )

    # data = pd.DataFrame({
    #     "density correlation": positive_correlations + lymph_node_correlations,
    #     "category": ["N+"] * len(positive_correlations) + ["LN"] * len(lymph_node_correlations)
    # })

    # data = pd.DataFrame({
    #     "density correlation": negative_correlations + lymph_node_correlations,
    #     "category": ["N0"] * len(negative_correlations) + ["LN"] * len(lymph_node_correlations)
    # })

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.boxplot(
        x="category",
        y="density correlation",
        data=data,
        color="lightgray",
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        x="category",
        y="density correlation",
        data=data,
        color="black",
        jitter=True,
        size=7,
        alpha=0.7,
        ax=ax,
    )

    y_max = ax.get_ylim()[1] * 0.95
    ax.text(
        0.5,
        y_max,
        f"$p$-value = {p_value:.4f}",
        ha="center",
        va="center",
        fontsize=12,
        color="red",
    )

    os.makedirs(f"../output/{mode}/density_correlation_analysis", exist_ok=True)
    fig.savefig(
        f"../output/{mode}/density_correlation_analysis/{marker_one}_{marker_two}_tile_size_{tile_size}_boxplot.pdf"
    )


if __name__ == "__main__":
    main()
