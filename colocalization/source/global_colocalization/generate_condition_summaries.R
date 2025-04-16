library(argparse)
library(devtools)
library(dplyr)
library(tidyr)

load_all("metadisco")

#' computes cell counts for all samples in a directory
#'
#' @param files list of cell assignment .csv files from celesta
#' @param save_path path to save cell counts of form ("{base_directory}/counts/{condition}")
#' @param cell_types list of cell types to include in the count matrix excluding unknown
#'
#' @examples
#' compute_cell_counts(
#'   c("./condition_one/assignments/521_S1_reg024_assignments.csv",
#'     "./condition_one/assignments/521_S1_reg028_assignments.csv"),
#'   "./counts/condition_one",
#'   c("epithelial cell (e-cadherin+), endothelial cell (CD31+))
#' )
compute_cell_counts <- function(files, save_path, cell_types) {
    cell_types <- c("unknown", cell_types)
    
    count_matrix <- lapply(files, function(file) {
        sample_name <- sub("_assignments\\.csv$", "", basename(file))

        sample_assignments <- read.csv(file) %>%
            count(`FINAL_CELL_TYPE`) %>%
            complete(`FINAL_CELL_TYPE` = factor(`FINAL_CELL_TYPE`, levels = cell_types), fill = list(n = 0)) %>%
            rename(!!sample_name := n)
        
        return(sample_assignments)
    })
    
    count_matrix <- Reduce(function(x, y) {
        full_join(x, y, by = "FINAL_CELL_TYPE")
    }, count_matrix)
    
    write.csv(count_matrix, sprintf("%s_count_matrix.csv", save_path), row.names = FALSE)
}

#' computes global colocalization statistics for all samples in a directory
#'
#' @param files list of cell assignment .csv files from celesta
#' @param save_path path to save colocalization matrix of form ("{base_directory}/global_colocalizations/{condition}")
#' @param n_neighbors number of nearest neighbors to use for colocalization analysis; \
#'                    defaults to 15
#' @param bandwidth radius of circular neighborhood to consider; \
#'                  defaults to 100
#' @param weight_scheme weight scheme to use for colocalization analysis; \
#'                      options include (constant, linear, squared_exponential);
#'                      defaults to linear
#' @param log_transform indicates whether to log-plus-one transform the colocalization statistics; \
#'                      defaults to FALSE
#'
#' @examples
#' compute_global_colocalization(
#'   c("./condition_one/assignments/521_S1_reg024_assignments.csv",
#'     "./condition_one/assignments/521_S1_reg028_assignments.csv"),
#'   "./global_colocalizations/condition_one"
#' )
compute_global_colocalization <- function(files, save_path, n_neighbors = 15, bandwidth = 100, weight_scheme = "linear", log_transform = FALSE) {
    colocalization_matrix <- lapply(files, function(file) {
        sample_name <- sub("_assignments\\.csv$", "", basename(file))
        
        sample_assignments <- read.csv(file) %>% 
            rename(x = `X`, y = `Y`, cell_type = `FINAL_CELL_TYPE`) %>% 
            mutate(sample = sample_name)
        
        spomic <- createSpomic(sample_assignments) %>%
            setSpomicHypers(tile_size = 250, window_size = 10, k_neigh = n_neighbors, 
                                bandwidth = bandwidth, n_bootstrap = 100, weight_scheme = weight_scheme)
        
        colocalization <- getCLQs(spomic) %>%
            mutate(!!sample_name := if (log_transform) log(`colocalization_stat` + 1) else `colocalization_stat`) %>%
            select(A_B, !!sample_name)
        
        return(colocalization)
    })
    
    colocalization_matrix <- Reduce(function(x, y) {
        return(full_join(x, y, by = "A_B"))
    }, colocalization_matrix)
    
    transformation_label <- if (log_transform) "transformed" else "untransformed"
    write.csv(colocalization_matrix, sprintf("%s_global_colocalization_matrix_%s_%s_weighting_%d_neighbors_%d_pixel_bandwidth.csv", 
                                                 save_path, transformation_label, weight_scheme, n_neighbors, bandwidth), row.names = FALSE)
}

#' parses command line arguments
#'
#' @returns a list of parsed arguments
#'
#' @examples
#' parse_arguments()
parse_arguments <- function() {
    parser <- ArgumentParser(description = "batch global CLQ processing of cell assignment files; \
                                            recommended to use (weight_scheme = constant, n_neighbors = 10) or (weight_scheme = linear / squared_exponential, n_neighbors = 15)")
    
    parser$add_argument("--bandwidth", type = "integer", default = 100, help = "bandwidth radius of circular neighborhood to consider; \
                                                                                defaults to 100")
    parser$add_argument("--data_directory", help = "path to a data directory of cell assignment .csv files", required = TRUE)
    parser$add_argument("--filter", default = "all", help = "comma-separated list of sample names to process, or 'all' to process everything; \
                                                             defaults to 'all'")
    parser$add_argument("--log_transform", action = "store_true", help = "whether to log-transform the colocalization statistics")
    parser$add_argument("--n_neighbors", type = "integer", default = 15, help = "number of nearest neighbors to use for global colocalization analysis; \
                                                                                 defaults to 15")
    parser$add_argument("--save_path", help = "path to save conglomerate pairwise global CLQs and cell counts", required = TRUE)
    parser$add_argument("--signature_matrix", help = "path to a signature matrix file", required = TRUE)
    parser$add_argument("--weight_scheme", default = "linear", help = "weight scheme to use for colocalization analysis; \
                                                                       options include (constant, linear, squared_exponential); \
                                                                       defaults to linear")
    
    args <- parser$parse_args()
    
    return(args)
}

#' reads in arguments and generates pairwise global CLQ values across cell types
#'
#' @examples
#' main()
main <- function() {
    args <- parse_arguments()
    
    bandwidth <- args$bandwidth
    data_directory <- args$data_directory
    filter <- args$filter
    log_transform <- args$log_transform
    n_neighbors <- args$n_neighbors
    save_path <- args$save_path
    signature_matrix <- args$signature_matrix
    weight_scheme <- args$weight_scheme

    sample_assignments <- list.files(data_directory, pattern = "\\.csv$", full.names = TRUE)
    signature_matrix <- read.csv(signature_matrix)
    
    if (filter != "all") {
        filter <- strsplit(filter, ",")[[1]]
        sample_assignments <- sample_assignments[basename(sample_assignments) %in%
                                                     paste0(filter, "_assignments.csv")]
    }
    
    lapply(c("counts", "global_colocalizations"), function(dir) {
        dir_path <- file.path(save_path, dir)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
    })
    
    compute_cell_counts(sample_assignments, file.path(save_path, "counts", basename(data_directory)), signature_matrix$CELL_TYPE)
    compute_global_colocalization(sample_assignments, file.path(save_path, "global_colocalizations", basename(data_directory)), n_neighbors, bandwidth, weight_scheme, log_transform)
}

main()