# metadisco
## Meta-DiSCO: A statistical framework for differential cell-type colocalization across groups in single-cell spatial omics
![alt text](https://github.com/plevritis-lab/metadisco/blob/main/images/metadisco_figure1.png)

## Installation
The metadisco package is currently private within the plevritis-lab GitHub group. For lab members, the easiest way to test and use metadisco is to clone this repository somewhere on your local computer and then call the following at the beginning of your scripts: 
```r
devtools::load_all("<PATH-TO-CLONED-REPOSITORY>")
```

## Usage
metadisco largely relies on two data structures, the spatial omic (spomic) object, and the metadisco object.
### spomic
A spomic object can be created either through an existind data frame or by providing a path to a csv file. The data must contain columns for x, y, cell_type, and sample.
```r
# Initiate spomic via data frame
set.seed(123)
n <- 100
example_df <- data.frame(
  x = rgamma(n, 20, 2),
  y = rgamma(n, 20, 2),
  cell_type = sample(x = c("A", "B", "C"), size = n, replace = TRUE),
  sample = rep("example", n)
)
spomic <- createSpomic(example_df)

# Initiate spomic via path to .csv
# spomic <- createSpomic(<PATH-TO-CSV>)
```

The metadisco package contains useful data visualization functions to explore your data. 
```r
# Scatter plot of your spatial omic data, colored by cell type
plotSpomic(spomic, point_size = 0.5)

# Bar plot of cell type proportions
plotCellProportions(spomic, show_pct = TRUE)
```

Importantly, to run downstream analysis on a spomic object, you must set some crucial hyperparameters. 
```r
# Beware, these hyperparameters are just for the example and are not necessarily defaults for analysis on actual data 
spomic <- setSpomicHypers(
  spomic_obj = spomic,
  k_neigh = 5, # Number of neighbors used in CLQ calculation
  bandwidth = 10, # Max radius for neighbors in CLQ calculation 
  weight_scheme = "linear", # Weight scheme for neighbors in CLQ calculation
  tile_size = 2, # Size of square tiles to use in spatial bootstrap
  window_size = 1, # Number of adjacent tiles to include in spatial bootstrap
  n_bootstrap = 100, # Number of times to bootstrap
  precompute_neighbors = TRUE # Whether to identify nearest neighbors upon instantiation
)
```

After setting hyperparameters, we can call crucial functions such as calculating the CLQ between every unique cell-type pair. 
```r
clq <- getCLQs(spomic)
print(clq)
```

### metadisco
The metadisco class is instantiated by creating lists of spomic objects based on a grouping variable 
```r
set.seed(123)
group1_list <- list()
group2_list <- list()

n <- 100
for(i in 1:5) {
  n <- 100
  group1_df <- data.frame(
    x = rgamma(n, 20, 2),
    y = rgamma(n, 20, 2),
    cell_type = sample(x = c("A", "B", "C"), size = 100, replace = TRUE),
    sample = rep(paste("Group 1, Sample", i), n)
  )
  group1_spomic <- createSpomic(group1_df)
  group1_spomic <- setSpomicHypers(
    spomic_obj = group1_spomic,
    k_neigh = 5, 
    bandwidth = 10, 
    weight_scheme = "linear", 
    tile_size = 2, 
    window_size = 1, 
    n_bootstrap = 10, 
    precompute_neighbors = TRUE 
  )
  group1_list[[i]] <- group1_spomic
  
  group2_df <- data.frame(
    x = rgamma(n, 30, 2),
    y = rgamma(n, 30, 2),
    cell_type = sample(x = c("A", "B", "C"), size = 100, replace = TRUE),
    sample = rep(paste("Group 2, Sample", i), n)
  )
  group2_spomic <- createSpomic(group2_df)
  group2_spomic <- setSpomicHypers(
    spomic_obj = group1_spomic,
    k_neigh = 5, 
    bandwidth = 10, 
    weight_scheme = "linear", 
    tile_size = 2, 
    window_size = 1, 
    n_bootstrap = 10, 
    precompute_neighbors = TRUE 
  )
  group2_list[[i]] <- group2_spomic
}

# Warning: this function does the bulk of the metadisco pipeline so it takes a few minutes
metadisco <- createMetaDisco(
  condition1 = group1_list, 
  condition2 = group2_list, 
  k_neigh = 5, # Number of neighbors used in CLQ calculation
  bandwidth = 10, # Max radius for neighbors in CLQ calculation 
  weight_scheme = "linear", # Weight scheme for neighbors in CLQ calculation
  tile_size = 2, # Size of square tiles to use in spatial bootstrap
  window_size = 1, # Number of adjacent tiles to include in spatial bootstrap
  n_bootstrap = 10, 
  n_cores = parallel::detectCores()
)

```

