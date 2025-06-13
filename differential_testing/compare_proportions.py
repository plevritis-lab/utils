import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scipy.spatial.distance import braycurtis
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
import seaborn as sns
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

# == BASIC DATA LOADING AND PREPARATION ==

colormap = "/Users/rohit/Desktop/EPOCHS/celesta/source/colormaps/plevritis/colormap_primary_plevritis.json"

centers_negative = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/colocalizations/counts/assignments_count_matrix.csv"
centers_positive = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_positive/colocalizations/counts/assignments_count_matrix.csv"
edges_negative = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/edges/node_negative/colocalizations/counts/assignments_count_matrix.csv"
edges_positive = "/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/edges/node_positive/colocalizations/counts/assignments_count_matrix.csv"

centers_negative = pd.read_csv(centers_negative)
centers_positive = pd.read_csv(centers_positive)
centers = pd.merge(centers_negative, centers_positive, on = "FINAL_CELL_TYPE")

edges_negative = pd.read_csv(edges_negative)
edges_positive = pd.read_csv(edges_positive)
edges = pd.merge(edges_negative, edges_positive, on = "FINAL_CELL_TYPE")

with open(colormap, "r") as file:
    colormap = json.load(file)

def plot_global_proportions(centers, edges, colormap):
    color_mapping = {k: v["color"] for k, v in colormap.items()}
    
    centers_counts = centers.drop(columns = ["FINAL_CELL_TYPE"]).values.sum(axis = 1)
    edges_counts = edges.drop(columns = ["FINAL_CELL_TYPE"]).values.sum(axis = 1)
    centers_proportions = centers_counts / centers_counts.sum()
    edges_proportions = edges_counts / edges_counts.sum()

    fig, axes = plt.subplots(1, 2, figsize = (8, 8))

    for ax in axes:
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks([])

        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.tick_params(axis = "y", length = 0)

        ax.bar([0], [1], color = "white", edgecolor = "black", linewidth = 1, width = 0.5)
    
    axes[0].set_yticks([0.25, 0.5, 0.75])
    axes[0].set_yticklabels(["25%", "50%", "75%"], fontsize = 8, weight = "bold", va = "center")
    axes[0].spines[["top", "right", "left", "bottom"]].set_visible(False)
    axes[0].tick_params(axis = "y", length = 0)
    
    bottom = 0
    for proportion, color, cell_type in zip(centers_proportions, color_mapping.values(), color_mapping.keys()):
        axes[0].bar([0], [proportion], color = color, edgecolor = "black", linewidth = 1, 
                   width = 0.5, bottom = bottom, label = cell_type)
        bottom += proportion
    
    bottom = 0
    for proportion, color, cell_type in zip(edges_proportions, color_mapping.values(), color_mapping.keys()):
        axes[1].bar([0], [proportion], color = color, edgecolor = "black", linewidth = 1, 
                   width = 0.5, bottom = bottom, label = cell_type)
        bottom += proportion
    
    legend_fig, legend_ax = plt.subplots(figsize = (5, len(color_mapping) * 0.4))

    legend_ax.set_axis_off()

    handles, labels = axes[0].get_legend_handles_labels()
    legend_ax.legend(
        handles = handles,
        labels = labels,
        loc = "center",
        facecolor = "white",
        edgecolor = "black",
        labelcolor = "black",
        markerscale = 1.5,
        prop = {"size": 10, 
                "weight": "bold"}
    )

    legend_fig.tight_layout()
    fig.tight_layout()

    fig.savefig("/Users/rohit/Downloads/figures/proportions.png", dpi = 300)
    legend_fig.savefig("/Users/rohit/Downloads/figures/proportions_legend.png", dpi = 300)

    plt.close(legend_fig)
    plt.close(fig)

plot_global_proportions(centers, edges, colormap)

# == PERMANOVA TESTING ==

def run_permanova(centers, edges):
    centers = centers.set_index("FINAL_CELL_TYPE").T
    edges = edges.set_index("FINAL_CELL_TYPE").T

    combined = pd.concat([centers, edges])
    condition = ["centers"] * len(centers) + ["edges"] * len(edges)

    distances = np.array([
        braycurtis(u, v) for u in combined.values for v in combined.values
    ]).reshape(len(combined), len(combined))

    distance_matrix = DistanceMatrix(distances, ids = combined.index)

    results = permanova(distance_matrix, grouping = condition)

    return condition, distance_matrix, results

def plot_principal_coordinates_analysis(condition, distance_matrix):
    results = pcoa(distance_matrix)
    coords = results.samples.iloc[:, :2].copy()
    coords["condition"] = condition

    var_explained = results.proportion_explained * 100
    x_label = f"pcoa1 ({var_explained.iloc[0]:.2f}% variance)"
    y_label = f"pcoa2 ({var_explained.iloc[1]:.2f}% variance)"

    fig, ax = plt.subplots(figsize = (6, 5))
    sns.scatterplot(
        data = coords,
        x = coords.columns[0],
        y = coords.columns[1],
        hue = "condition",
        s = 100,
        edgecolor = "black",
        ax = ax
    )

    ax.set_title("pcoa (bray-curtis distance)")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    fig.tight_layout()
    
    fig.savefig("/Users/rohit/Downloads/figures/pcoa.png", dpi = 300)

    plt.close(fig)

condition, distance_matrix, results = run_permanova(centers, edges)
plot_principal_coordinates_analysis(condition, distance_matrix)

print(results)

# == AUXILIARIES ==

def test_cell_type_differences(centers, edges):
    centers = centers.set_index("FINAL_CELL_TYPE")
    edges = edges.set_index("FINAL_CELL_TYPE")

    celltypes = centers.index
    pvals = []

    for ct in celltypes:
        x = centers.loc[ct].values
        y = edges.loc[ct].values
        _, p = mannwhitneyu(x, y, alternative="two-sided")
        pvals.append(p)

    # Multiple testing correction
    corrected = multipletests(pvals, method="fdr_bh")
    
    results_df = pd.DataFrame({
        "cell_type": celltypes,
        "p_value": pvals,
        "adj_p_value": corrected[1],
        "significant": corrected[0],
    })

    return results_df.sort_values("adj_p_value")

# Barplot visualization of significant cell types
def plot_significant_barplot(centers, edges, results_df):
    significant = results_df[results_df["significant"]]
    significant_types = significant["cell_type"].tolist()

    centers = centers.set_index("FINAL_CELL_TYPE")
    edges = edges.set_index("FINAL_CELL_TYPE")

    mean_centers = centers.loc[significant_types].mean(axis=1)
    mean_edges = edges.loc[significant_types].mean(axis=1)

    df = pd.DataFrame({
        "cell_type": significant_types,
        "centers": mean_centers,
        "edges": mean_edges
    }).reset_index(drop=True).melt(id_vars="cell_type", var_name="group", value_name="proportion")

    plt.figure(figsize=(10, len(significant_types) * 0.4))
    sns.barplot(data=df, x="proportion", y="cell_type", hue="group")
    plt.title("significantly different cell types (fdr < 0.05)")
    plt.tight_layout()
    plt.savefig("/Users/rohit/Downloads/figures/significant_barplot.png", dpi=300)
    plt.close()

results = test_cell_type_differences(centers, edges)

plot_significant_barplot(centers, edges, results)

# def plot_volcano(results_df, centers, edges):
#     centers = centers.set_index("FINAL_CELL_TYPE")
#     edges = edges.set_index("FINAL_CELL_TYPE")

#     mean_centers = centers.mean(axis=1)
#     mean_edges = edges.mean(axis=1)
#     logfc = np.log2((mean_edges + 1e-6) / (mean_centers + 1e-6))

#     results_df["log2_fc"] = logfc.values
#     results_df["-log10_p"] = -np.log10(results_df["adj_p_value"] + 1e-10)

#     plt.figure(figsize=(10, 7))
#     ax = sns.scatterplot(
#         data=results_df,
#         x="log2_fc",
#         y="-log10_p",
#         hue="significant",
#         s=80
#     )
#     plt.axhline(-np.log10(0.05), color="red", linestyle="--", label="FDR = 0.05")
#     plt.axvline(0, color="gray", linestyle="--")

#     # Add cell type labels to each point
#     for i, row in results_df.iterrows():
#         ax.text(
#             row["log2_fc"], row["-log10_p"], row["cell_type"],
#             fontsize=8, ha="right", va="bottom"
#         )

#     plt.xlabel("Log2 Fold Change (edges / centers)")
#     plt.ylabel("-Log10 Adjusted p-value")
#     plt.title("Differential Abundance by Cell Type")
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig("/Users/rohit/Downloads/figures/volcano_celltype_diff.png", dpi=300)
#     plt.close()

# def plot_pcoa_biplot(centers, edges, num_top_features=10):
#     centers = centers.set_index("FINAL_CELL_TYPE").T.rename_axis("SAMPLE")
#     edges = edges.set_index("FINAL_CELL_TYPE").T.rename_axis("SAMPLE")

#     combined = pd.concat([centers, edges])
#     condition = ["centers"] * len(centers) + ["edges"] * len(edges)

#     # Compute Bray-Curtis distance matrix
#     distances = np.array([braycurtis(u, v) for u in combined.values for v in combined.values])
#     distances = distances.reshape(len(combined), len(combined))
#     distance_matrix = DistanceMatrix(distances, ids=combined.index)

#     # PCoA
#     pcoa_results = pcoa(distance_matrix)
#     coords = pcoa_results.samples.iloc[:, :2].copy()
#     coords["condition"] = condition

#     # Feature (cell type) loadings: correlate with PCoA axes
#     feature_corr = {}
#     for feature in combined.columns:
#         rho1, _ = spearmanr(combined[feature], coords.iloc[:, 0])
#         rho2, _ = spearmanr(combined[feature], coords.iloc[:, 1])
#         feature_corr[feature] = (rho1, rho2)

#     feature_df = pd.DataFrame(feature_corr, index=["PCoA1_corr", "PCoA2_corr"]).T
#     feature_df["magnitude"] = np.sqrt(feature_df["PCoA1_corr"]**2 + feature_df["PCoA2_corr"]**2)
#     feature_df["abs_corr"] = feature_df[["PCoA1_corr", "PCoA2_corr"]].abs().max(axis=1)
#     # top_features = feature_df.nlargest(num_top_features, "magnitude")
#     top_features = feature_df[feature_df["abs_corr"] > 0.3]
#     # top_features = feature_df

#     # Plot
#     fig, ax = plt.subplots(figsize=(7, 6))
#     sns.scatterplot(data=coords, x=coords.columns[0], y=coords.columns[1], hue="condition", s=100, ax=ax, edgecolor="black")

#     for cell_type, row in top_features.iterrows():
#         ax.arrow(0, 0, row["PCoA1_corr"], row["PCoA2_corr"], color="gray", alpha=0.7, head_width=0.02)
#         ax.text(row["PCoA1_corr"] * 1.1, row["PCoA2_corr"] * 1.1, cell_type, fontsize=8, ha="center", va="center")

#     var_explained = pcoa_results.proportion_explained * 100
#     ax.set_xlabel(f"PCoA1 ({var_explained[0]:.2f}%)")
#     ax.set_ylabel(f"PCoA2 ({var_explained[1]:.2f}%)")
#     ax.set_title("PCoA Biplot (Bray-Curtis)")
#     ax.axhline(0, color="black", lw=0.5)
#     ax.axvline(0, color="black", lw=0.5)
#     ax.legend()
#     plt.tight_layout()
#     plt.savefig("/Users/rohit/Downloads/figures/pcoa_biplot.png", dpi=300)
#     plt.close()

# results = test_cell_type_differences(centers, edges)
# plot_volcano(results, centers, edges)
# plot_pcoa_biplot(centers, edges)

# print(results)