#' MetaDisco Class
#'
#' @description
#' Defines the structure for the `MetaDisco` class, which is used to represent metadata
#' and analysis parameters for comparative spatial omics studies. This class
#' supports the storage of confitions, spatial tiling parameters, windowing,
#' bootstrap configurations, and analysis results
#'
#' @slot condition1 A list of `Spomic` classes for samples within the first group/condition.
#' @slot condition2 A list of `Spomic` classes for samples within the second group/condition.
#' @slot tile_size A numeric value indicating the size of tiles to be used in the analysis and bootstrapping.
#' @slot window_size An integer specifying the size of the spatial window for analysis.
#' @slot k_neigh An integer representing the number of nearest neighbors to consider
#' in the colocalization calculations.
#' @slot bandwidth An integer specifying the maximum distance from a point to consider in colocalization calculations.
#' @slot n_bootstrap An integer denoting the number of bootstrapping iterations to perform.
#' @slot weight_scheme A character string describing the weighting scheme to be applied for colocalization calculations.
#' @slot results A data frame containing the results of differential colocalization testing.
#'
#'
#' @export
setClass("MetaDisco", slots = list(
  condition1 = "list",
  condition2 = "list",
  tile_size = "numeric",
  window_size = "integer",
  k_neigh = "integer",
  bandwidth = "integer",
  n_bootstrap = "integer",
  weight_scheme = "character",
  results = "data.frame"
  )
)


#' @export
setMetaDiscoHypers <- function(metadisco_obj,
                               tile_size,
                               k_neigh,
                               window_size,
                               bandwidth,
                               n_bootstrap,
                               weight_scheme){
  metadisco_obj@tile_size <- as.numeric(tile_size)
  metadisco_obj@window_size <- as.integer(window_size)
  metadisco_obj@k_neigh <- as.integer(k_neigh)
  metadisco_obj@bandwidth <- as.integer(bandwidth)
  metadisco_obj@n_bootstrap <- as.integer(n_bootstrap)
  metadisco_obj@weight_scheme <- as.character(weight_scheme)

  return(metadisco_obj)
}

#' #' @export
#' createMetaDisco <- function(condition1,
#'                             condition2,
#'                             n_tile,
#'                             window_size,
#'                             k_neigh,
#'                             bandwidth,
#'                             n_bootstrap,
#'                             bootstrap_method,
#'                             weight_scheme,
#'                             n_cores = parallel::detectCores()){
#'   list1 <- list()
#'   for(i in 1:length(condition1)){
#'     print(condition1[i])
#'     if(is.list(condition1)){
#'       list1[[i]] <- createSpomic(condition1[[i]])
#'     } else{
#'       list1[[i]] <- createSpomic(condition1[i])
#'     }
#'     # list1[[i]] <- createSpomic(condition1[i])
#'     list1[[i]] <- setSpomicHypers(spomic_obj=list1[[i]],
#'                                   n_tile=n_tile,
#'                                   window_size=window_size,
#'                                   k_neigh=k_neigh,
#'                                   bandwidth=bandwidth,
#'                                   n_bootstrap=n_bootstrap,
#'                                   weight_scheme=weight_scheme)
#'     list1[[i]] <- runCLQ(spomic_obj=list1[[i]],
#'                          bootstrap_method=bootstrap_method,
#'                          n_cores=n_cores)
#'   }
#'
#'
#'   list2 <- list()
#'   for(i in 1:length(condition2)){
#'     print(condition2[i])
#'     if(is.list(condition2)){
#'       list2[[i]] <- createSpomic(condition2[[i]])
#'     } else{
#'       list2[[i]] <- createSpomic(condition2[i])
#'     }
#'     # list2[[i]] <- createSpomic(condition2[i])
#'     list2[[i]] <- setSpomicHypers(spomic_obj=list2[[i]],
#'                                   n_tile=n_tile,
#'                                   window_size=window_size,
#'                                   k_neigh=k_neigh,
#'                                   bandwidth=bandwidth,
#'                                   n_bootstrap=n_bootstrap,
#'                                   weight_scheme=weight_scheme)
#'     list2[[i]] <- runCLQ(spomic_obj=list2[[i]],
#'                          bootstrap_method=bootstrap_method,
#'                          n_cores=n_cores)
#'   }
#'
#'
#'
#'   object <- new("MetaDisco", condition1 = list1, condition2 = list2)
#'   names(object@condition1) <- sapply(object@condition1, function(x) x@sample_id)
#'   names(object@condition2) <- sapply(object@condition2, function(x) x@sample_id)
#'
#'   object <- setMetaDiscoHypers(object,
#'                                n_tile=n_tile,
#'                                k_neigh=k_neigh,
#'                                window_size=window_size,
#'                                bandwidth=bandwidth,
#'                                n_bootstrap=n_bootstrap,
#'                                weight_scheme=weight_scheme)
#'   return(object)
#' }


#' Create a MetaDisco Object
#'
#' @description
#' Constructs a `MetaDisco` object from two conditions, each containing one or
#' more `Spomic` objects or data respresentations that can be converted into `Spomic`
#' objects. This function processes spatial omics data for both conditions, sets
#' analysis parameters (hypers), and performs colocalization analysis with
#' optional bootstrapping.
#'
#' @param condition1 A list of `Spomic` objects, data frames or paths to CSV files for the first condition.
#' @param condition1 A list of `Spomic` objects, data frames of paths to CSV files for the second condition.
#' @param tile_size A numeric value specifying the size of tiles to be used for spatial bootstrapping.
#' @param window_size An integer specifying the size of the spatial window for analysis.
#' @param k_neigh An integer representing the number of nearest neighbors to consider in colocalization analysis.
#' @param bandwidth An integer defining the bandwidth used to restrict colocalization analysis.
#' @param n_bootstrap An integer specifying the number of bootstrapping iterations to perform.
#' @param weight_scheme A character string describing the weighting scheme for colocalization.
#' @param n_cores An integer specifying the number of CPU cores to use for parallel processing.
#' Defaults to the number of available cores.
#'
#' @return A `MetaDisco` object containing processed spatial omics data for both conditions, along with the analysis parameters.
#'
#' @details
#' The function processes the input data for each condition using the following steps:
#' 1. Converts data representations (e.g., data frames or CSV paths) into `Spomic` objects if necessary.
#' 2. Sets hyperparameters (`tile_size`, `window_size`, `k_neigh`, etc.) for each `Spomic` object.
#' 3. Runs colocalization analysis (`runCLQ`) for each `Spomic` object.
#' 4. Combines processed `Spomic` objects into `MetaDisco` object with all analysis parameters set.
#'
#' @export
createMetaDisco <- function(condition1,
                            condition2,
                            tile_size,
                            window_size,
                            k_neigh,
                            bandwidth,
                            n_bootstrap,
                            weight_scheme,
                            n_cores = parallel::detectCores()) {
  # Helper function to process a condition
  process_condition <- function(condition) {
    condition_list <- list()
    for (i in seq_along(condition)) {
      print(condition[[i]])
      if (inherits(condition[[i]], "Spomic")) {
        spomic_obj <- condition[[i]]
      } else if (is.list(condition[[i]])) {
        spomic_obj <- createSpomic(condition[[i]])
      } else {
        spomic_obj <- createSpomic(condition[i])
      }

      spomic_obj <- setSpomicHypers(spomic_obj = spomic_obj,
                                    tile_size = tile_size,
                                    window_size = window_size,
                                    k_neigh = k_neigh,
                                    bandwidth = bandwidth,
                                    n_bootstrap = n_bootstrap,
                                    weight_scheme = weight_scheme)
      spomic_obj <- runCLQ(spomic_obj = spomic_obj,
                           bootstrap_method = bootstrap_method,
                           n_cores = n_cores)
      condition_list[[i]] <- spomic_obj
    }
    return(condition_list)
  }

  # Process both conditions
  list1 <- process_condition(condition1)
  list2 <- process_condition(condition2)

  # Create MetaDisco object
  object <- new("MetaDisco", condition1 = list1, condition2 = list2)
  names(object@condition1) <- sapply(object@condition1, function(x) x@sample_id)
  names(object@condition2) <- sapply(object@condition2, function(x) x@sample_id)

  # Set MetaDisco hypers
  object <- setMetaDiscoHypers(object,
                               tile_size = tile_size,
                               k_neigh = k_neigh,
                               window_size = window_size,
                               bandwidth = bandwidth,
                               n_bootstrap = n_bootstrap,
                               weight_scheme = weight_scheme)
  return(object)
}

#' @export
saveMetaDisco <- function(metadisco_obj, condition1_name, condition2_name, path, saveRDS = TRUE){
  path_w_hypers <- file.path(path, paste0("tilesize", metadisco_obj@tile_size,
                                          "_nbors", metadisco_obj@k_neigh,
                                          "_bandwidth", metadisco_obj@bandwidth,
                                          "_bootstrap", metadisco_obj@n_bootstrap,
                                          "_window", metadisco_obj@window_size),
                             paste0(metadisco_obj@weight_scheme, "_weights")
                             )
  # make_dir(path_w_hypers)
  dir.create(path_w_hypers, recursive = TRUE)
  print(path_w_hypers)
  if(saveRDS) saveRDS(metadisco_obj, file.path(path_w_hypers, paste0(condition1_name, "_vs_", condition2_name, ".Rds")))
}

getMetaOutputs <- function(metadisco_obj, condition = c("condition1", "condition2"), plot_forest = TRUE, output_path = "output"){
  # Temporary (need to fix this problem upstream of the function)
  names(metadisco_obj@condition1) <- paste("sample_", as.character(1:length(metadisco_obj@condition1)), sep = "")
  names(metadisco_obj@condition2) <- paste("sample_", as.character(1:length(metadisco_obj@condition2)), sep = "")
  ###

  condition1_df <- data.frame()
  for(i in 1:length(metadisco_obj@condition1)){
    condition1_df <- rbind(condition1_df, metadisco_obj@condition1[[i]]@processed_colocalization |> dplyr::mutate(sample = names(metadisco_obj@condition1)[i]))
  }

  cell_pairs <- unique(condition1_df$A_B)
  coef <- list()
  problematic_pairs <- list()
  for(cell_pair in cell_pairs){
    df <- condition1_df |> dplyr::filter(A_B == cell_pair) |> tidyr::drop_na() |> dplyr::filter(proportion > 0.5)
    if(nrow(df) > 1) { # Check if df has any rows
      # meta_model <- metafor::rma.uni(yi=df$mean, sei=df$se, method="SJ")
      warning_message <- NULL

      # Use tryCatch to capture the model and handle warnings
      meta_model <- tryCatch(
        {
          withCallingHandlers(
            {
              model <- metafor::rma.uni(yi=df$mean, sei=df$se, method="SJ")
              model
            },
            warning = function(w) {
              # Store the warning message
              warning_message <<- conditionMessage(w)
              # Continue running the model even if there's a warning
              invokeRestart("muffleWarning")
            }
          )
        },
        error = function(e) {
          # Store error messages if necessary
          warning_message <<- conditionMessage(e)
          NA  # Return NA for the model if there's an error
        }
      )

      # Store the problematic pair if there's a warning message
      if (!is.null(warning_message)) {
        problematic_pairs[[cell_pair]] <- warning_message
      }
      if(plot_forest){
        fp <- file.path(output_path,
                        paste0("tilesize", metadisco_obj@tile_size,
                               "_nbors", metadisco_obj@k_neigh,
                               "_bandwidth", metadisco_obj@bandwidth,
                               "_bootstrap", metadisco_obj@n_bootstrap,
                               "_window", metadisco_obj@window_size),
                        paste0(metadisco_obj@weight_scheme, "_weights"),
                        "forest_plots")
        dir.create(fp, recursive = TRUE)
        png(file=file.path(fp,
                           paste0(condition[1], "_", cell_pair, '_forestplot.png')))
            metafor::forest.rma(meta_model)
            dev.off()
      }
      coef[[cell_pair]] <- coef(summary(meta_model))
    } else {
      coef[[cell_pair]] <- NA # Assign NA or some other indicator for no data
    }
  }
  print("Cell pairs with large variance ratio issues:")
  print(problematic_pairs)
  condition1_coef <- purrr::map_df(coef, ~dplyr::as_tibble(.x), .id = "cell_pair") |>
    dplyr::mutate(condition = "condition1")|>
    dplyr::select(cell_pair, estimate, se, zval, pval, ci.lb, ci.ub, condition)
  # |>
  #   dplyr::filter(!cell_pair %in% names(problematic_pairs))


  condition2_df <- data.frame()
  for(i in 1:length(metadisco_obj@condition1)){
    condition2_df <- rbind(condition2_df, metadisco_obj@condition2[[i]]@processed_colocalization |> dplyr::mutate(sample = names(metadisco_obj@condition2)[i]))
  }

  cell_pairs <- unique(condition2_df$A_B)
  coef <- list()
  problematic_pairs <- list()
  for(cell_pair in cell_pairs){
    df <- condition2_df |> dplyr::filter(A_B == cell_pair) |> tidyr::drop_na() |> dplyr::filter(proportion > 0.5)
    if(nrow(df) > 1) { # Check if df has any rows
      # meta_model <- metafor::rma.uni(yi=df$mean, sei=df$se, method="SJ")
      warning_message <- NULL

      # Use tryCatch to capture the model and handle warnings
      meta_model <- tryCatch(
        {
          withCallingHandlers(
            {
              model <- metafor::rma.uni(yi=df$mean, sei=df$se, method="SJ")
              model
            },
            warning = function(w) {
              # Store the warning message
              warning_message <<- conditionMessage(w)
              # Continue running the model even if there's a warning
              invokeRestart("muffleWarning")
            }
          )
        },
        error = function(e) {
          # Store error messages if necessary
          warning_message <<- conditionMessage(e)
          NA  # Return NA for the model if there's an error
        }
      )

      # Store the problematic pair if there's a warning message
      if (!is.null(warning_message)) {
        problematic_pairs[[cell_pair]] <- warning_message
      }
      if(plot_forest){
        png(file=file.path(output_path,
                           paste0("tilesize", metadisco_obj@tile_size,
                                  "_nbors", metadisco_obj@k_neigh,
                                  "_bandwidth", metadisco_obj@bandwidth,
                                  "_bootstrap", metadisco_obj@n_bootstrap,
                                  "_window", metadisco_obj@window_size),
                           paste0(metadisco_obj@weight_scheme, "_weights"),
                           "forest_plots",
                           paste0(condition[2], "_", cell_pair, '_forestplot.png')))
        metafor::forest.rma(meta_model)
        dev.off()
      }
      coef[[cell_pair]] <- coef(summary(meta_model))
    } else {
      coef[[cell_pair]] <- NA # Assign NA or some other indicator for no data
    }
  }
  print("Cell pairs with large variance ratio issues:")
  print(problematic_pairs)
  condition2_coef <- purrr::map_df(coef, ~dplyr::as_tibble(.x), .id = "cell_pair") |>
    dplyr::mutate(condition = "condition2")|>
    dplyr::select(cell_pair, estimate, se, zval, pval, ci.lb, ci.ub, condition)
  # |>
  #   dplyr::filter(!cell_pair %in% names(problematic_pairs))

  return(rbind(condition1_coef, condition2_coef))
}

#' @export
calcZ <- function(coef1, coef2){
  if(is.null(coef1$estimate) | is.null(coef2$estimate) | is.null(coef1$se) | is.null(coef2$se)){
    break
  } else{
    return((coef2$estimate - coef1$estimate)/sqrt((coef2$se)^2 + (coef1$se)^2))
  }
}


#' @export
calcLog2 <- function(coef1, coef2){
  if(is.null(coef1$estimate) | is.null(coef2$estimate)){
    return(NA)
    break
  } else{
    # return(log2((coef1$estimate+1)/ (coef2$estimate+1)))
    return(log2(coef2$estimate+1) - log2(coef1$estimate+1))

  }
}

#' @export
runMeta <- function(metadisco_obj, plot_forest=TRUE, output_path = "output"){
  randeff_coef <- getMetaOutputs(metadisco_obj, condition = c("condition1", "condition2"), plot_forest=plot_forest, output_path = output_path)
  cell_pairs <- unique(randeff_coef$cell_pair)
  cond1 <- randeff_coef |> dplyr::filter(condition == "condition1")
  cond2 <- randeff_coef |> dplyr::filter(condition == "condition2")

  zscores <- list()
  for(pair in cell_pairs){
    zscores[[pair]] <- dplyr::tibble(A_B = pair,
                                     zscore = calcZ(cond1 |> dplyr::filter(cell_pair == pair),
                                                    cond2 |> dplyr::filter(cell_pair == pair)),
                                     log2fc = calcLog2(cond1 |> dplyr::filter(cell_pair == pair),
                                                       cond2 |> dplyr::filter(cell_pair == pair)))
  }
  results <- dplyr::bind_rows(zscores)
  results <- results |> dplyr::mutate(pval = (1-pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = TRUE)) * 2,
                                      FDR = p.adjust(pval, method="BH"),
                                      Bonferroni = p.adjust(pval, method="bonferroni")
  )

  metadisco_obj@results <- results

  return(metadisco_obj)
}
