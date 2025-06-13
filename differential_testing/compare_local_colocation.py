import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

# Load the data
df = pd.read_csv("/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/local_spatial_summary/microarray_clusters.csv")

# Step 1: Compute proportion of each cluster within each cell type (ignores condition for now)
cluster_counts = (
    df.groupby(['FINAL_CELL_TYPE', 'CLUSTER'])  # count cells per (cell type × cluster)
    .size()
    .unstack(fill_value=0)  # clusters become columns; missing = 0
)

# Step 2: Convert raw counts to proportions
cluster_props = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)

# Step 3: Plot stacked barplot of cluster proportions by cell type
colors = sns.color_palette("pastel", n_colors=cluster_props.shape[1])

cluster_props.plot(
    kind='barh',
    stacked=True,
    color=colors,
    figsize=(12, 6)
)
plt.title("Proportion of Clusters by Cell Type")
plt.xlabel("Proportion of Cells")
plt.ylabel("Cell Type")
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("/Users/rohit/Downloads/figures/cluster_proportion_by_cell_type.png", dpi=300)

results = []
residuals_dict = {}

# Loop through each cell type
for cell_type, group in df.groupby("FINAL_CELL_TYPE"):
    # Contingency table: rows = clusters, columns = conditions (e.g., REGION)
    contingency = pd.crosstab(group['CLUSTER'], group['REGION'])

    # Skip if only one condition present (can't test across conditions)
    if contingency.shape[1] < 2:
        continue

    # 2x2 case: use Fisher’s Exact Test
    if contingency.shape == (2, 2):
        _, p = stats.fisher_exact(contingency)
        test_type = 'fisher'
        residuals = None  # Fisher’s test does not give residuals
    else:
        # Chi-squared test
        chi2, p, dof, expected = stats.chi2_contingency(contingency)
        test_type = 'chi2'

        # Standardized residuals
        residuals = (contingency - expected) / np.sqrt(expected)
        residuals_dict[cell_type] = residuals

    # Save result
    results.append({
        'FINAL_CELL_TYPE': cell_type,
        'p_value': p,
        'n_clusters': contingency.shape[0],
        'n_conditions': contingency.shape[1],
        'test_type': test_type
    })

# Convert to DataFrame and correct for multiple testing
results_df = pd.DataFrame(results)
results_df['p_adj'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
results_df['log10_p'] = -np.log10(results_df['p_adj'] + 1e-10)
results_df['significant'] = results_df['p_adj'] < 0.05

# Plot summary barplot
plt.figure(figsize=(10, 5))
sns.barplot(
    data=results_df.sort_values('p_adj'),
    x='log10_p',
    y='FINAL_CELL_TYPE',
    hue='significant',
    dodge=False,
    palette={True: 'crimson', False: 'gray'}
)
plt.axvline(-np.log10(0.05), color='black', linestyle='--', label='fdr = 0.05')
plt.xlabel('-log10(fdr-adjusted p-value)')
plt.ylabel('cell type ')
plt.title('cluster distribution differences across conditions')
plt.legend(title='significant')
plt.tight_layout()
plt.savefig("/Users/rohit/Downloads/figures/cluster_distribution_differences_across_conditions.png", dpi=300)

# (OPTION 1) - SEQUENTIAL # Plot residuals for significant Chi2-tested cell types
# for cell_type in results_df[(results_df['significant']) & (results_df['test_type'] == 'chi2')]['FINAL_CELL_TYPE']:
#     res = residuals_dict[cell_type]

#     plt.figure(figsize=(8, 4))
#     sns.heatmap(res, annot=True, center=0, cmap='vlag', fmt=".2f")
#     plt.title(f"{cell_type} - Standardized Residuals (Cluster x Condition)")
#     plt.xlabel('Condition')
#     plt.ylabel('Cluster')
#     plt.tight_layout()
#     plt.savefig(f"/Users/rohit/Downloads/figures/residuals_{cell_type.replace(' ', '_')}.png", dpi=300)
#     plt.close()

# (OPTION 2) - SINGLE HEATMAP Combine all significant chi2 residuals into a single heatmap
# sig_chi2_types = results_df[(results_df['significant']) & (results_df['test_type'] == 'chi2')]['FINAL_CELL_TYPE']

# # Collect residuals into a MultiIndex DataFrame: (Cell Type, Cluster) x Condition
# residuals_list = []
# for cell_type in sig_chi2_types:
#     res = residuals_dict[cell_type]
#     # Add cell type as a column for later stacking
#     res = res.copy()
#     res['FINAL_CELL_TYPE'] = cell_type
#     residuals_list.append(res.reset_index())

# if residuals_list:
#     all_residuals = pd.concat(residuals_list, ignore_index=True)
#     # Melt for heatmap: rows = (Cell Type, Cluster), columns = Condition
#     melted = all_residuals.melt(
#         id_vars=['FINAL_CELL_TYPE', 'CLUSTER'],
#         var_name='CONDITION',
#         value_name='RESIDUAL'
#     )
#     # Remove rows where CONDITION is 'FINAL_CELL_TYPE' (from melt)
#     melted = melted[melted['CONDITION'] != 'FINAL_CELL_TYPE']

#     # Pivot for heatmap
#     heatmap_data = melted.pivot_table(
#         index=['FINAL_CELL_TYPE', 'CLUSTER'],
#         columns='CONDITION',
#         values='RESIDUAL'
#     )

#     plt.figure(figsize=(min(20, 2 + 0.5 * heatmap_data.shape[0]), 1 + 0.5 * heatmap_data.shape[0]))
#     sns.heatmap(
#         heatmap_data,
#         annot=False,
#         center=0,
#         cmap='vlag',
#         cbar_kws={'label': 'Standardized Residual'}
#     )
#     plt.title("Standardized Residuals for Significant Cell Types (Cluster x Condition)")
#     plt.xlabel('Condition')
#     plt.ylabel('Cell Type / Cluster')
#     plt.tight_layout()
#     plt.savefig("/Users/rohit/Downloads/figures/all_significant_residuals_heatmap.png", dpi=300)
#     plt.close()

# (OPTION 3) Flatten all residuals into a long DataFrame
residual_rows = []

for cell_type, residual_df in residuals_dict.items():
    for cluster in residual_df.index:
        for region in residual_df.columns:
            residual_rows.append({
                'cell_type': cell_type,
                'cluster': cluster,
                'region': region,
                'residual': residual_df.loc[cluster, region]
            })

residual_long_df = pd.DataFrame(residual_rows)
residual_long_df['cell_cluster'] = residual_long_df['cell_type'] + " | cluster " + residual_long_df['cluster'].astype(str)

# Pivot for heatmap
heatmap_df = residual_long_df.pivot(index='cell_cluster', columns='region', values='residual')

# Sort rows by magnitude (optional)
heatmap_df = heatmap_df.loc[heatmap_df.abs().max(axis=1).sort_values(ascending=False).index]

# Plot all residuals on one heatmap
plt.figure(figsize=(12, max(6, 0.3 * len(heatmap_df))))
sns.heatmap(
    heatmap_df,
    cmap='vlag',
    center=0,
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={'label': 'Standardized Residual'}
)
plt.title("all cell type × cluster × condition residuals")
plt.xlabel("condition")
plt.ylabel("cell type | cluster")
plt.tight_layout()
plt.savefig("/Users/rohit/Downloads/figures/all_celltype_cluster_condition_residuals_heatmap.png", dpi=300)
plt.close()