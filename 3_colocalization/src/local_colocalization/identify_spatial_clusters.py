import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
import umap


def concatenate_spatial_embeddings(condition):
    """concatenates LCLQ spatial embeddings across all images from a user-specified condition (only "complete" embeddings are retained)

    args:
        condition (str): parent directory (e.g., center/negative) that contains sample subdirectories where spatial_embeddings.csv is located
    """

    condition_spatial_embeddings = []

    for root, _, files in os.walk(condition):
        for file in files:
            if "local_colocalization_matrix" in file:
                path_parts = root.split(os.path.sep)
                regionalization = path_parts[-4]
                nodal_status = path_parts[-3]
                sample_name = re.match(
                    r"([A-Za-z0-9_]+)_local_colocalization_matrix", file
                ).group(1)

                cell_embedding_matrix = pd.read_csv(os.path.join(root, file))

                cell_embedding_matrix["REGION"] = regionalization
                cell_embedding_matrix["N"] = nodal_status
                cell_embedding_matrix["SAMPLE"] = sample_name

                condition_spatial_embeddings.append(cell_embedding_matrix)

    condition_spatial_embeddings = pd.concat(
        condition_spatial_embeddings, ignore_index=True
    )

    return condition_spatial_embeddings


def plot_data_distributions(A_cells, A, save_path):
    """plots and saves a series of images visualizing LCLQ distributions between all cell type pairs involving A

    args:
        A_cells (dataframe): generated LCLQ measurements for every A cell
        A (str): cell type of interest (e.g., epithelial)
        save_path (str): parent directory path where output (histograms) will be stored under their own subdirectory (distributions/A)
    """

    A_cells = A_cells.drop(columns=["REGION", "N", "X", "Y", "SAMPLE"])

    save_path = os.path.join(save_path, f"distributions/{A}")
    os.makedirs(save_path, exist_ok=True)

    for cell_type in A_cells.columns:
        fig, ax = plt.subplots()

        ax.hist(A_cells[cell_type])

        fig.savefig(os.path.join(save_path, f"histogram_{cell_type}.png"), dpi=300)
        plt.close(fig)


def reduce_data_dimensionality(A_cells, A, color_by, save_path, method="PCA"):
    """reduces the dimensionality of A's LCLQ spatial embeddings to 2 by performing PCA or UMAP, saving the ensuing visualization

    args:
        A_cells (dataframe): generated LCLQ measurements for every A cell
        A (str): cell type of interest (e.g., epithelial)
        color_by (str): condition by which points in the PCA or UMAP visualization should be colored by (either REGION, N, or NONE)
        save_path (str): parent directory path where output (PCA or UMAP visualizations) will be stored under their own subdirectory (PCA or UMAP)
        method (str, optional): method that conducts dimensionality reduction (either PCA or UMAP); \
                                defaults to PCA
    """

    regionalization = A_cells["REGION"].values
    nodal_status = A_cells["N"].values

    # print(regionalization)
    # print(nodal_status)

    A_cells = A_cells.drop(columns=["REGION", "N", "X", "Y", "SAMPLE"])
    A_cells = A_cells.fillna(0)

    if method == "PCA":
        reduced_representation = PCA(n_components=2).fit_transform(A_cells)

    elif method == "UMAP":
        reduced_representation = umap.UMAP(
            n_neighbors=min(A_cells.shape[0] - 1, 15)
        ).fit_transform(A_cells)

    fig, ax = plt.subplots()

    if color_by != "NONE":
        if color_by == "REGION":
            subset_one = reduced_representation[np.where(regionalization == "centers")]
            subset_two = reduced_representation[np.where(regionalization == "edges")]

            labels = ["centers", "edges"]

        elif color_by == "N":
            subset_one = reduced_representation[
                np.where(nodal_status == "node_positive")
            ]
            subset_two = reduced_representation[
                np.where(nodal_status == "node_negative")
            ]

            labels = ["node_positive", "node_negative"]

        ax.scatter(
            subset_one[:, 0], subset_one[:, 1], color="blue", alpha=0.5, label=labels[0]
        )
        ax.scatter(
            subset_two[:, 0], subset_two[:, 1], color="red", alpha=0.5, label=labels[1]
        )

        ax.legend()

    else:
        ax.scatter(reduced_representation[:, 0], reduced_representation[:, 1])

    save_path = os.path.join(save_path, method)
    os.makedirs(save_path, exist_ok=True)

    fig.savefig(os.path.join(save_path, f"{method}_{A}_{color_by}.png"), dpi=300)

    plt.close(fig)


def cluster_cell_type_embeddings(A_cells, A, save_path, method="k_means"):
    """clusters LCLQ spatial embeddings for a user-provided cell type (amenable to both per-image and per-condition clustering), saving results from an accompanying silhouette and elbow analysis to determine optimal k

    args:
        A_cells (dataframe): generated LCLQ measurements for every A cell; \
                             modified in place through the addition of a CLUSTER column
        A (str): cell type of interest (e.g., epithelial)
        save_path (str): parent directory path where output (silhouette plots, elbow curves) will be stored under their own subdirectory (k_means or DBSCAN)
        method (str, optional): method that carries out clustering (either k_means or DBSCAN); \
                                defaults to k_means
    """

    # TODO - do you try and scale the values to be on a consistent scale (like StandardScaler() with scikit learn or just keep the raw values which may have some intuitive meaning?)
    # TODO - improve quality of clusters (silhouette score is very low, maybe try a different clustering algorithm or the elbow method - do the latter?!)

    A_embeddings = A_cells.drop(columns=["REGION", "N", "X", "Y", "SAMPLE"])
    A_embeddings = A_embeddings.fillna(0)

    if method == "k_means":
        K = list(range(2, min(A_embeddings.shape[0], 11)))

        silhouette_scores = [None] * len(K)
        squared_errors = [None] * len(K)

        for i, k in enumerate(K):
            kmeans = KMeans(n_clusters=k, random_state=42).fit(A_embeddings)

            silhouette_scores[i] = silhouette_score(A_embeddings, kmeans.labels_)
            squared_errors[i] = kmeans.inertia_

        if silhouette_scores == []:
            A_cells["CLUSTER"] = 0

            return

        optimal_k = K[np.argmax(silhouette_scores)]

        fig, [ax_one, ax_two] = plt.subplots(1, 2, figsize=(12, 5))

        ax_one.plot(K, silhouette_scores, marker="o")
        ax_one.axvline(
            x=optimal_k, color="red", linestyle="--", label=f"optimal k = {optimal_k}"
        )

        ax_one.set_xlabel("k")
        ax_one.set_ylabel("silhouette score")

        ax_two.plot(K, squared_errors, marker="o")
        ax_two.axvline(
            x=optimal_k, color="red", linestyle="--", label=f"optimal k = {optimal_k}"
        )

        ax_two.set_xlabel("k")
        ax_two.set_ylabel("sum of squared errors")

        ax_one.legend()
        ax_two.legend()

        save_path = os.path.join(save_path, method)
        os.makedirs(save_path, exist_ok=True)

        fig.savefig(os.path.join(save_path, f"{method}_{A}.pdf"))

        plt.close(fig)

        cluster_labels = KMeans(n_clusters=optimal_k, random_state=42).fit_predict(
            A_embeddings
        )
        A_cells["CLUSTER"] = cluster_labels

    elif method == "DBSCAN":
        min_samples = 2 * A_embeddings.shape[1]

        nearest_neighbors = NearestNeighbors(n_neighbors=min_samples).fit(A_embeddings)
        neighbor_distances, _ = nearest_neighbors.kneighbors(A_embeddings)
        neighbor_distances = np.sort(neighbor_distances[:, -1], axis=0)

        n_points = len(neighbor_distances)
        all_points = np.arange(len(neighbor_distances))

        line_vector = np.array([n_points - 1, neighbor_distances[-1]]) - np.array(
            [0, neighbor_distances[0]]
        )
        line_vector_norm = line_vector / np.sqrt(np.sum(line_vector**2))

        vector_from_first = np.vstack([all_points, neighbor_distances]).T - np.array(
            [0, neighbor_distances[0]]
        )
        scalar_product = np.sum(
            vector_from_first * np.tile(line_vector_norm, (n_points, 1)), axis=1
        )
        vector_from_line = vector_from_first - np.outer(
            scalar_product, line_vector_norm
        )
        distance_to_line = np.sqrt(np.sum(vector_from_line**2, axis=1))

        optimal_index = np.argmax(distance_to_line)
        optimal_epsilon = neighbor_distances[optimal_index]

        fig, ax = plt.subplots()

        ax.plot(neighbor_distances)
        ax.axhline(
            y=optimal_epsilon,
            color="r",
            linestyle="--",
            label=f"optimal epsilon = {optimal_epsilon:.3f}",
        )
        ax.scatter(
            optimal_index, optimal_epsilon, color="red", marker="o", label="elbow point"
        )

        ax.set_xlabel("data points")
        ax.set_ylabel("distance to {}-th nearest neighbor".format(min_samples))
        ax.set_title("k-distance graph")

        ax.legend()

        save_path = os.path.join(save_path, method)
        os.makedirs(save_path, exist_ok=True)

        fig.savefig(os.path.join(save_path, f"{method}_{A}.pdf"))

        plt.close(fig)

        cluster_labels = DBSCAN(
            eps=optimal_epsilon, min_samples=min_samples
        ).fit_predict(A_embeddings)
        A_cells["CLUSTER"] = cluster_labels


def save_spatial_clusters(A_cells, A, save_path):
    """plots and saves the identified colocalization states, mapping them back to their original location in the image, as well as a LCLQ cluster summary spreadsheet

    args:
        A_cells (dataframe): generated LCLQ measurements and cluster labels for every A cell
        A (str): cell type of interest (e.g., epithelial)
        save_path (str): parent directory path where output (spatial cluster maps, summary spreadsheets) will be stored under their own subdirectory (cluster_maps, cluster_summaries)
    """

    cluster_maps = os.path.join(save_path, "cluster_maps")
    cluster_summaries = os.path.join(save_path, "cluster_summaries")

    os.makedirs(cluster_maps, exist_ok=True)
    os.makedirs(cluster_summaries, exist_ok=True)

    for sample_name, sample_A_cells in A_cells.groupby("SAMPLE"):
        fig, ax = plt.subplots(figsize=(10, 8))

        for cluster_label in sorted(sample_A_cells["CLUSTER"].unique()):
            cluster_points = sample_A_cells[sample_A_cells["CLUSTER"] == cluster_label]

            # TODO - add MEM-like enrichment profile for each cluster label (top three most enriched LCLQ profiles, top three lowest enriched LCLQ profiles)

            plt.scatter(
                cluster_points["X"],
                cluster_points["Y"],
                s=20,
                label=f"{A} cluster {cluster_label}",
            )

        ax.set_facecolor("black")

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.legend(loc="center", bbox_to_anchor=(0.5, -0.1))

        fig.tight_layout()

        sample_cluster_maps = os.path.join(cluster_maps, sample_name)
        os.makedirs(sample_cluster_maps, exist_ok=True)

        fig.savefig(os.path.join(sample_cluster_maps, f"{A}_clusters.pdf"))

        plt.close(fig)

    A_embeddings = A_cells.drop(columns=["REGION", "N", "X", "Y", "SAMPLE"])
    summary_statistics = A_embeddings.groupby("CLUSTER").mean()

    summary_statistics.to_csv(
        os.path.join(cluster_summaries, f"{A}_cluster_summary.csv")
    )


def cluster_spatial_embeddings(condition, log_normalization=True):
    """clusters LCLQ spatial embeddings for all cell types across a series of images from a user-specified condition

    args:
        condition (str): parent directory (e.g., center/negative) that contains sample subdirectories where spatial_embeddings.csv is located
        log_normalization (bool, optional): applies log + 1 transformation to LCLQ values to account for distributional skewness; \
                                            defaults to True (recommended to always enable this setting)
    """

    matplotlib.use("agg")

    cell_embedding_matrix = concatenate_spatial_embeddings(condition)

    condition = os.path.join(condition, "local_spatial_summary")
    os.makedirs(condition, exist_ok=True)

    if log_normalization:
        normalize_columns = cell_embedding_matrix.columns.difference(
            ["FINAL_CELL_TYPE", "REGION", "N", "X", "Y", "SAMPLE"]
        )
        cell_embedding_matrix[normalize_columns] = np.log1p(
            cell_embedding_matrix[normalize_columns]
        )

    cluster_dataframes = []

    unique_cell_types = cell_embedding_matrix["FINAL_CELL_TYPE"].unique()
    for A in unique_cell_types:
        A_cells = cell_embedding_matrix[
            cell_embedding_matrix["FINAL_CELL_TYPE"] == A
        ].drop(columns="FINAL_CELL_TYPE")

        # plot_data_distributions(A_cells, A, condition)
        # reduce_data_dimensionality(A_cells, A, "REGION", condition, method = "PCA")
        # reduce_data_dimensionality(A_cells, A, "REGION", condition, method = "UMAP")
        # reduce_data_dimensionality(A_cells, A, "N", condition, method = "PCA")
        # reduce_data_dimensionality(A_cells, A, "N", condition, method = "UMAP")
        cluster_cell_type_embeddings(A_cells, A, condition, method="k_means")
        # cluster_cell_type_embeddings(A_cells, A, condition, method = "DBSCAN")
        save_spatial_clusters(A_cells, A, condition)

        A_cells["FINAL_CELL_TYPE"] = A

        cluster_dataframes.append(A_cells)

    cluster_dataframes = pd.concat(cluster_dataframes, ignore_index=True)
    cluster_dataframes.to_csv(
        os.path.join(condition, "microarray_clusters.csv"), index=False
    )


if __name__ == "__main__":
    # conditions = ["../output/microarray_cores/center/negative",
    #               "../output/microarray_cores/center/positive",
    #               "../output/microarray_cores/edge/negative",
    #               "../output/microarray_cores/edge/positive",
    #               "../output/whole_slide_images/involved_lymph_nodes"]

    # for condition in conditions:
    #     for sub_directory in glob.glob(os.path.join(condition, "*")):
    #         print("computing LCLQ spatial clusters for sample", os.path.basename(sub_directory))

    #         cluster_spatial_embeddings(sub_directory)

    cluster_spatial_embeddings(
        "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors"
    )

"""
    # TODO - THINGS TO TRY MYSELF FOR SANITY CHECKS
    # implement the following ideas:
    # try non-negative matrix factorization and map to different regions
    # try downsampling to improve the clustering outcomes
    # try z score normalization after taking the log of the data
    # bring distributions into a reasonable range
    # normalize data within sample first since N varies across images so much
    # maybe try on the whole slide images to see if the clusters come out more cleanly (proof of concept nothing is wrong with underlying idea)
    # jake experimented around with scaling the data (z-scoring the data) and then clustered on that (recommended to log and then z-score)
    # maybe simulate conditions where there is real difference in the tma cores just to prove that the metric works
    # look at the cross k function for computing a spatial embedding for a cell?! - generates a vector of measurements and may be better; is this idea novel and can be improved on?
    # can you compare LCLQ values across images given that the denominator varies across so many images? maybe plot a distribution of the denominator or maybe just get rid of the denominator entirely to enable easier comparison
    # convert the LCLQ values into probabilities by applying softmax perhaps?

    # in summary:
    # 1. try on the whole slide image just to try as a proof of concept (can it do something not on TMA) or on simulated data of TMA cores
    # 2. try z score normalization
    # 3. try non negative matrix factorization (maybe not gonna get you farther than PCA)
    # 4. try something other than a LCLQ (graph neural network)
    # 5. try applying softmax to the LCLQ embedding to generate probability measurements - could be interesting for normalization??? remove denominator???
    # 6. other different metrics besides silhouette since they are naturally always very low
    # 7. look into different weighting functions (gaussian exponential decay, but what about others, plot curve as function of distance and bandwidth)
    # 8. try on simulated data with real colocalization potential
"""
