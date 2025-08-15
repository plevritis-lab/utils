import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

# ========== 1. LOAD DATA ==========

def load_persistence_data(directory):
    data = []
    for fname in os.listdir(directory):
        if fname.endswith(".json"):
            with open(os.path.join(directory, fname)) as f:
                content = json.load(f)
                for entry in content:
                    cell_type = entry.get("cell_type", "unknown")
                    for pdg in entry.get("persistence", []):
                        birth, death = pdg["birth"], pdg["death"]
                        if death is None:
                            continue
                        data.append({
                            "file": fname,
                            "cell_type": cell_type,
                            "dimension": pdg["dimension"],
                            "birth": birth,
                            "death": death,
                            "persistence": death - birth
                        })
    return pd.DataFrame(data)

center_df = load_persistence_data("/Users/rohit/Downloads/topology/centers")
edge_df = load_persistence_data("/Users/rohit/Downloads/topology/edges")

center_df["condition"] = "center"
edge_df["condition"] = "edge"
full_df = pd.concat([center_df, edge_df], ignore_index=True)

# ========== 2. COMPUTE SUMMARY STATS ==========

summary = (
    full_df.groupby(["condition", "cell_type", "dimension"])
    .agg(
        count=("persistence", "count"),
        mean_persistence=("persistence", "mean"),
        max_persistence=("persistence", "max"),
        total_persistence=("persistence", "sum")
    )
    .reset_index()
)

# ========== 3. STATISTICAL TESTING ==========

def compare_conditions(df, feature):
    results = []
    for cell_type in df["cell_type"].unique():
        for dim in df["dimension"].unique():
            center_vals = df[(df["condition"] == "center") &
                             (df["cell_type"] == cell_type) &
                             (df["dimension"] == dim)][feature]
            edge_vals = df[(df["condition"] == "edge") &
                           (df["cell_type"] == cell_type) &
                           (df["dimension"] == dim)][feature]
            if len(center_vals) >= 5 and len(edge_vals) >= 5:
                stat, pval = mannwhitneyu(center_vals, edge_vals, alternative='two-sided')
                results.append({
                    "cell_type": cell_type,
                    "dimension": dim,
                    "feature": feature,
                    "center_mean": center_vals.mean(),
                    "edge_mean": edge_vals.mean(),
                    "p_value": pval
                })
    return pd.DataFrame(results)

comparison_df = compare_conditions(full_df, "persistence")

# ========== 4. VISUALIZE ALL CELL TYPES & DIMENSIONS ==========

def plot_all_celltypes(df, save_dir="/Users/rohit/Downloads/topology/plots"):
    os.makedirs(save_dir, exist_ok=True)
    for cell_type in df["cell_type"].unique():
        for dim in [0, 1]:
            subset = df[(df["cell_type"] == cell_type) & (df["dimension"] == dim)]
            if subset["condition"].nunique() < 2 or subset.empty:
                continue
            # Clip persistence values at 1st and 99th percentiles to reduce long tails
            lower = subset["persistence"].quantile(0.01)
            upper = subset["persistence"].quantile(0.99)
            clipped = subset.copy()
            clipped["persistence"] = clipped["persistence"].clip(lower, upper)
            plt.figure(figsize=(8, 5))
            sns.violinplot(data=clipped, x="condition", y="persistence", inner="box", cut=0)
            plt.title(f"{cell_type} - H{dim} Persistence Distribution")
            plt.ylabel("Persistence (death - birth)")
            plt.xlabel("Condition")
            plt.tight_layout()
            filename = f"{cell_type.replace(' ', '_').replace('(', '').replace(')', '')}_H{dim}.png"
            plt.savefig(os.path.join(save_dir, filename))
            plt.close()

plot_all_celltypes(full_df)

# ========== 5. FLAG MOST DIVERGENT FEATURES ==========

comparison_df["abs_diff"] = np.abs(comparison_df["center_mean"] - comparison_df["edge_mean"])
top_differences = (
    comparison_df[comparison_df["p_value"] < 0.05]
    .sort_values(by="abs_diff", ascending=False)
    .head(10)
)

print("Top 10 divergent cell type/dimension combinations (significant differences):")
print(top_differences)

# ========== 6. PCA + UMAP ON SUMMARY FEATURES ==========

# Pivot to wide format
wide_summary = summary.pivot_table(
    index=["condition", "cell_type"],
    columns="dimension",
    values=["mean_persistence", "max_persistence", "total_persistence"],
    fill_value=0
).reset_index()

# Flatten multiindex columns
wide_summary.columns = ['_'.join(map(str, col)).strip('_') for col in wide_summary.columns.values]

# Extract metadata and features
metadata_cols = [col for col in wide_summary.columns if col.startswith("condition") or col.startswith("cell_type")]
features = wide_summary.drop(columns=metadata_cols, errors='ignore')
X = StandardScaler().fit_transform(features)

# PCA
pca = PCA(n_components=2)
pca_proj = pca.fit_transform(X)
wide_summary["PCA1"] = pca_proj[:, 0]
wide_summary["PCA2"] = pca_proj[:, 1]

# UMAP
umap_proj = umap.UMAP(n_components=2, random_state=42).fit_transform(X)
wide_summary["UMAP1"] = umap_proj[:, 0]
wide_summary["UMAP2"] = umap_proj[:, 1]
# PCA plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=wide_summary, x="PCA1", y="PCA2", hue="condition", style="cell_type")
plt.title("PCA of Persistence Summary Features")
plt.tight_layout()
plt.savefig("/Users/rohit/Downloads/topology/pca_persistence_summary.png")
plt.close()

# UMAP plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=wide_summary, x="UMAP1", y="UMAP2", hue="condition", style="cell_type")
plt.title("UMAP of Persistence Summary Features")
plt.tight_layout()
plt.savefig("/Users/rohit/Downloads/topology/umap_persistence_summary.png")
plt.close()