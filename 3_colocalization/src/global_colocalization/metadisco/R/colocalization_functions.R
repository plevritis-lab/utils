#' Get nearest neighbors for a Spomic object
#'
#' @description
#' Computes the nearest neighbors for each spatial point in a `Spomic` object.
#'
#' @param spomic_obj A `Spomic` object containing spatial coordinates.
#'
#' @return A list where each element represents the neighbors for a given point,
#' indexed by row number.
#'
#' @details
#' The function retrieves the nearest `k` neighbors for each spatial point in `spomic@coords`.
#'
#' @importFrom spdep knearneigh
#'
#' @export
getNeighbors <- function(spomic_obj){
  # Validate input
  stopifnot(inherits(spomic_obj, "Spomic"))

  coords <- spomic_obj@coords
  k <- spomic_obj@k_neigh

  # Ensure the number of neighbors requested is valid
  if (k >= nrow(coords)) {
    stop("k_neigh must be smaller than the number of data points.")
  }

  knn <- spdep::knearneigh(spomic_obj@coords, k = spomic_obj@k_neigh)
  nb_list <- split(knn$nn, row(knn$nn))

  return(nb_list)
}

#' Search for cell neighbors
#'
#' @description
#' Computes the proportion of neighbors of type `cellB` for each cell of type `cellA `
#' in a `Spomic` object. This function supports different weighting schemes for calculating neighbor proportions.
#'
#' @param spomic_obj A `Spomic` object.
#' @param cellA A character string specifying the type of cell (e.g., "A") to search for as the primary cell.
#' @param cellB A character string specifying the type of cell (e.g., "B") to evaluate as neighbors of `cellA`.
#'
#' @return A matrix where rows correspond to cells of type `cellA` and the single column contains the weighted
#' proportion of neighbors of type `cellB`.
#' If no cells of type `cellA` are found, the function returns `NULL`.
#'
#' @details
#' The function uses the `@nb_list` slot in the `Spomic` object to retrieve the
#' neighbors for each cell of type `cellA`.
#' Proportions are weighted based on the distance to neighbors, using the following schemes:
#' - `"constant"`: All neighbors are weighted equally.
#' - `"linear"`: Neighbors closer to the current cell are weighted higher, and the weights decay linearly.
#'
#' @export
searchCellNeighbors <- function(spomic_obj, cellA, cellB) {
  # Get indices of cellA in cell assignments
  cell_A_indices <- which(spomic_obj@cell_assignments == cellA)

  # Return NULL if no indices are found
  if (length(cell_A_indices) == 0) {
    return(NULL)
  }

  # Initialize the matrix to store proportions
  proportion_nb <- numeric(length(cell_A_indices))
  names(proportion_nb) <- cell_A_indices

  # Iterate over each index to calculate the proportion of neighbors
  for (j in seq_along(cell_A_indices)) {
    current_cell_ID <- cell_A_indices[j]
    neighbors <- spomic_obj@nb_list[[current_cell_ID]]
    neighbor_types <- spomic_obj@cell_assignments[neighbors]

    # Continue only if neighbors exist
    if (length(neighbors) > 0) {
      neighbor_ranks <- seq_along(neighbors)

      # Determine weights based on the scheme
      weights <- switch(
        spomic_obj@weight_scheme,
        "constant" = rep(1, length(neighbor_ranks)),
        "linear" = 1 - (neighbor_ranks - 1) / length(neighbors),
        "squared_exponential" = exp(-0.5*(neighbor_ranks^2)/(spomic_obj@k_neigh^2)),
        stop("Invalid weight_scheme")
      )

      total_weight <- sum(weights)
      cell_type_mask <- (neighbor_types == cellB)

      # Calculate the weighted proportion if any matching neighbors are found
      if (any(cell_type_mask)) {
        proportion_nb[j] <- sum(weights[cell_type_mask]) / total_weight
      } else {
        proportion_nb[j] <- 0
      }
    }
  }

  # Return the result as a matrix
  proportion_nb_matrix <- matrix(proportion_nb, nrow = length(cell_A_indices), ncol = 1, dimnames = list(cell_A_indices, cellB))
  return(proportion_nb_matrix)
}

#' @export
calcCLQ <- function(spomic_obj, cellA, cellB){
  cell_A_indices <- which(spomic_obj@cell_assignments == cellA)
  cell_B_indices <- which(spomic_obj@cell_assignments == cellB)
  if(length(cell_B_indices) <= 0 | length(cell_A_indices) <= 0){
    return(NA)
  } else{
    # total_cell <- length(unique(c(cell_A_indices, cell_B_indices)))
    total_cell <- nrow(spomic_obj@df)
    cell_A_nb_of_B <- searchCellNeighbors(spomic_obj, cellA, cellB)
    Cab <- sum(cell_A_nb_of_B, na.rm=TRUE)
    if(cellA==cellB){
      CLQ_ab <- (Cab/length(cell_A_indices))/((length(cell_B_indices)-1)/(total_cell-1))
    } else{
      CLQ_ab <- (Cab/length(cell_A_indices))/(length(cell_B_indices)/(total_cell-1))
    }
  }
  return(CLQ_ab)
}

#' @export
getCLQs <- function(spomic_obj){
  all_ranges <- c()
  CLQsub_all <- c()
  # all_cell_nb_in_bandwidth <- tryCatch(spdep::dnearneigh(as.matrix(spomic_obj@coords), 0, spomic_obj@bandwidth, longlat = NULL),error=function(e){cat("ERROR :",conditionMessage(e), "\n");skip_to_next <<- TRUE})
  # all_cell_nb_in_bandwidth <- tryCatch({
  #   # Attempt to create neighbor object
  #   nb <- spdep::dnearneigh(
  #     as.matrix(spomic_obj@coords),
  #     0,
  #     spomic_obj@bandwidth,
  #     longlat = NULL
  #   )
  #
  #   # Check if the neighbor object is empty
  #   if (length(nb) == 0 || all(sapply(nb, length) == 0)) {
  #     return(list())  # Return an empty list if no neighbors
  #   }
  #
  #   # Return the neighbor object if valid
  #   return(nb)
  # }, error = function(e) {
  #   cat("ERROR:", conditionMessage(e), "\n")
  #   return(list())  # Return an empty list on error
  # })
  all_cell_nb_in_bandwidth <- tryCatch({
    # Attempt to create neighbor object
    nb <- spdep::dnearneigh(
      as.matrix(spomic_obj@coords),
      0,
      spomic_obj@bandwidth,
      longlat = NULL
    )

  })

  all_cell_nb_in_bandwidth <- lapply(all_cell_nb_in_bandwidth, function(neighbors) {
    if (length(neighbors) == 1) {
      return(NA)  # Replace empty lists with NA
    } else {
      return(neighbors)  # Keep neighbors as is
    }
  })



  int_nb_list <- lapply(seq_along(spomic_obj@nb_list), function(i) {
    intersect(spomic_obj@nb_list[[i]], all_cell_nb_in_bandwidth[[i]])
  })

  cell_types <- unique(spomic_obj@cell_assignments)
  CLQ_matrix <- matrix(0L, nrow = length(cell_types), ncol = length(cell_types))
  rownames(CLQ_matrix) <- cell_types
  colnames(CLQ_matrix) <- cell_types

  for(k in 1:length(cell_types)){
    for(j in 1:length(cell_types)){
      cellA <- cell_types[k]
      cellB <- cell_types[j]
      CLQ_result <- calcCLQ(spomic_obj, cellA, cellB)

      CLQ_matrix[k,j] <- CLQ_result
    }
  }
  CLQsub <- as.data.frame(as.table(CLQ_matrix), stringsAsFactors = FALSE) |>
    dplyr::rename(A = Var1, B = Var2, CLQ = Freq)
  # CLQsub <- dplyr::tibble(Sample = sample, CLQsub)
  return(CLQsub |> tidyr::unite(A_B, A, B, sep="_"))
}

#' @export
#' @param data univariate vector of values to fit Gaussian Mixture Model
#' @param return_n boolean to return the number of components fit, if FALSE the function will return bool on whether the data is univariate or not
determineUnimodal <- function(data, return_n = TRUE){
  mog <- mclust::densityMclust(data, plot=FALSE)
  if(return_n){return(mog$G)}
  return(mog$G == 1)
}

#' @export
fitMoG <- function(colocalization_df){
  # colocalization_df should be a filtered dataframe of colocalization values for a given cell pair
  # it is colocalization of both original and bootstrapped values
  # print(unique(colocalization_df$A_B))
  original_stat <- colocalization_df |> dplyr::filter(sample == "Original") |> tidyr::drop_na() |> dplyr::pull(colocalization_stat)
  if(length(original_stat)==0){
    return(data.frame(mean=NA,
                      se=NA,
                      n_components=NA,
                      proportion=NA,
                      observed_colocalization_stat=NA))
  }

  data <- colocalization_df |> dplyr::filter(sample == "Bootstrap") |> tidyr::drop_na() |> dplyr::pull(colocalization_stat)
  n_boot <- length(data)
  if(length(unique(data)) <= 1){
    # return(data.frame(mean=NA,
    #                   se=NA,
    #                   n_components=NA,
    #                   proportion=NA,
    #                   observed_colocalization_stat=original_stat))
    return(data.frame(mean=mean(data),
                      se=1e-5,
                      n_components=NA,
                      proportion=NA,
                      observed_colocalization_stat=original_stat))
  }


  mog <- mclust::densityMclust(data, plot=FALSE)
  optimal_components <- mog$G

  if(optimal_components == 1){
    return(data.frame(mean=unname(mog$parameters$mean),
             se=sqrt(unname(mog$parameters$variance$sigmasq)/n_boot) + 1e-5,
             n_components=1,
             proportion=unname(mog$parameters$pro),
             observed_colocalization_stat=original_stat)
    )
  } else {
    means <- unname(mog$parameters$mean)
    vars <- unname(mog$parameters$variance$sigmasq)
    n <- n_boot * unname(mog$parameters$pro)
    ses <- vars/n
    return(data.frame(mean=means,
                      se=ses + 1e-5,
                      n_components=mog$G,
                      proportion=unname(mog$parameters$pro),
                      observed_colocalization_stat=original_stat))
  }
}

# postprocessCLQs <- function(colocalization_df, keep_unimodal=TRUE){
#   param_df <- colocalization_df |>
#     dplyr::group_by(A_B) |>
#     dplyr::do(fitMoG(.)) |>
#     dplyr::filter(proportion>0.8)
#     # dplyr::filter(n_components==1)
#   return(param_df)
# }

postprocessCLQs <- function(colocalization_df, n_bootstrap) {

  # Fit Mixture of Gaussians model by A_B group and calculate std
  param_df <- colocalization_df |>
    dplyr::group_by(A_B) |>
    dplyr::do(fitMoG(.)) |>
    dplyr::mutate(std = se * sqrt(proportion * n_bootstrap))

  # Get unique cell pairs
  cell_pairs <- unique(param_df$A_B)

  # Initialize a list to store results
  results_list <- vector("list", length(cell_pairs))

  # Loop over each cell pair to calculate bias and related stats
  for (i in seq_along(cell_pairs)) {
    cell_pair <- cell_pairs[i]
    mog_results <- param_df |> dplyr::filter(A_B == cell_pair)

    # Extract parameters
    means <- mog_results$mean
    stds <- mog_results$std
    ses <- mog_results$se
    proportions <- mog_results$proportion
    observed_colocalization_stat <- unique(mog_results$observed_colocalization_stat)

    # Check for NA values in means
    if (any(is.na(means)) | any(is.na(mog_results$n_components))) {
      results_list[[i]] <- data.frame(
        A_B = cell_pair,
        bias = NA,
        mean = NA,
        se = NA,
        n_components = NA,
        proportion = NA
        # ,
        # comment = "Spatial bootstrap unsuccessful"
      )
    } else {
      if (length(means) == 1) {
        ################ testing
        probs <- sapply(seq_along(means), function(j) {
          prob <- pnorm(observed_colocalization_stat, mean = means[j], sd = stds[j])
          min(prob, 1 - prob)  # Choose the less extreme probability
        })
        if (probs > 0.05) {
          results_list[[i]] <- data.frame(
            A_B = cell_pair,
            bias = means - observed_colocalization_stat,
            mean = observed_colocalization_stat,
            se = ses,
            n_components = length(probs),
            proportion = proportions
          )
        } else {
          results_list[[i]] <- data.frame(
            A_B = cell_pair,
            bias = means - observed_colocalization_stat,
            mean = observed_colocalization_stat,
            se = ses + sqrt(abs(means - observed_colocalization_stat)),
            n_components = length(probs),
            proportion = proportions
          )
          ############
        }} else {
        # Calculate probabilities for multi-component case
        probs <- sapply(seq_along(means), function(j) {
          prob <- pnorm(observed_colocalization_stat, mean = means[j], sd = stds[j])
          min(prob, 1 - prob)  # Choose the less extreme probability
        })


        # Find the best distribution (least extreme probability)
        best_distribution <- which.max(probs)

        if (probs[best_distribution] > 0.05) {
          results_list[[i]] <- data.frame(
            A_B = cell_pair,
            bias = means[best_distribution] - observed_colocalization_stat,
            mean = observed_colocalization_stat,
            se = ses[best_distribution],
            n_components = length(probs),
            proportion = proportions[best_distribution]
          )
        } else {
          results_list[[i]] <- data.frame(
            A_B = cell_pair,
            bias = means[best_distribution] - observed_colocalization_stat,
            mean = observed_colocalization_stat,
            se = ses[best_distribution] + sqrt(abs(means[best_distribution] - observed_colocalization_stat)),
            n_components = length(probs),
            proportion = proportions[best_distribution]
          )
          # # If no valid distribution is found (all probs < 0.05), return NA
          # results_list[[i]] <- data.frame(
          #   A_B = cell_pair,
          #   bias = NA,
          #   mean = NA,
          #   se = NA,
          #   n_components = length(probs),
          #   proportion = NA
          # )
        }
        #   ################
        #   # Case with only one component
        #   results_list[[i]] <- data.frame(
        #     A_B = cell_pair,
        #     bias = means - observed_colocalization_stat,
        #     mean = observed_colocalization_stat,
        #     se = mog_results$se,
        #     n_components = 1,
        #     proportion = proportions
        #     # ,
        #     # comment = ""
        #   )
        # }
      }
    }
  }

  # Combine all results into a single data frame
  results_df <- dplyr::bind_rows(results_list)

  return(results_df)
}


getCLQs <- function(spomic_obj){
  all_ranges <- c()
  CLQsub_all <- c()
  all_cell_nb_in_bandwidth <- tryCatch(spdep::dnearneigh(as.matrix(spomic_obj@coords), 0, spomic_obj@bandwidth, longlat = NULL),error=function(e){cat("ERROR :",conditionMessage(e), "\n");skip_to_next <<- TRUE})
  int_nb_list <- lapply(seq_along(spomic_obj@nb_list), function(i) {
    intersect(spomic_obj@nb_list[[i]], all_cell_nb_in_bandwidth[[i]])
  })

  cell_types <- unique(spomic_obj@cell_assignments)
  CLQ_matrix <- matrix(0L, nrow = length(cell_types), ncol = length(cell_types))
  rownames(CLQ_matrix) <- cell_types
  colnames(CLQ_matrix) <- cell_types
  for(k in 1:length(cell_types)){
    for(j in 1:length(cell_types)){
      cellA <- cell_types[k]
      cellB <- cell_types[j]
      CLQ_result <- calcCLQ(spomic_obj, cellA, cellB)
      CLQ_matrix[k,j] <- CLQ_result
    }
  }
  CLQsub <- as.data.frame(as.table(CLQ_matrix), stringsAsFactors = FALSE) |>
    dplyr::rename(A = Var1, B = Var2, colocalization_stat = Freq)
  return(CLQsub |> tidyr::unite(A_B, A, B, sep="_"))
}


export_functions <- function(cl){
  require(parallel)
  require(doParallel)
  clusterExport(cl, varlist = c(
    # "blockBootstrap",
                                "windowBootstrap",
                                "tileImage",
                                "filterTiles",
                                "createSpomic",
                                "setSpomicHypers",
                                "getNeighbors",
                                "getCLQs",
                                "calcCLQ",
                                "searchCellNeighbors"
  ))
  # Ensure the package is loaded in each worker
  clusterEvalQ(cl, devtools::load_all(path = "/Users/jacobchang/Lab/metadisco")) # this is going to need to change eventually?
}

runCLQ <- function(spomic_obj, bootstrap_method = "window", n_cores = parallel::detectCores()){
  require(parallel)
  require(doParallel)
  original_CLQs <- getCLQs(spomic_obj) |> dplyr::mutate(sample = "Original")
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  export_functions(cl)
  clusterExport(cl, varlist = c(
    # "blockBootstrap",
    "windowBootstrap", "setSpomicHypers", "getCLQs"), envir = environment())

  set.seed(123)
  CLQs_list <- foreach(b = 1:spomic_obj@n_bootstrap, .packages = c("dplyr")) %dopar% {
    # if(bootstrap_method == "block") bootstrap_obj <- blockBootstrap(spomic_obj)
    # else if(bootstrap_method == "window")
      bootstrap_obj <- windowBootstrap(spomic_obj)

    bootstrap_obj <- setSpomicHypers(bootstrap_obj,
                                     tile_size = spomic_obj@tile_size,
                                     window_size = spomic_obj@window_size,
                                     k_neigh = spomic_obj@k_neigh, bandwidth = spomic_obj@bandwidth,
                                     n_bootstrap = spomic_obj@n_bootstrap,
                                     weight_scheme = spomic_obj@weight_scheme)
    bootstrap_CLQs <- getCLQs(bootstrap_obj) |> dplyr::mutate(sample = "Bootstrap")
    list(bootstrap_CLQs)
  }
  all_CLQs <- rbind(original_CLQs, dplyr::bind_rows(CLQs_list))
  stopCluster(cl)

  spomic_obj@colocalization <- all_CLQs
  spomic_obj@processed_colocalization <- postprocessCLQs(spomic_obj@colocalization,spomic_obj@n_bootstrap)

  return(spomic_obj)
}
