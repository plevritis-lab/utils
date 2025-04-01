suppressMessages(library(Rmixmod))
suppressMessages(library(spData))
suppressMessages(library(sf))

library(argparse)
library(Rmixmod)
library(spdep)
library(zeallot)

library(CELESTA)

#' identifies cell types in spatial proteomics data using celesta
#'
#' @param quantified_imaging_data cell measurements in .csv format
#' @param signature_matrix celesta signature matrix in .csv format
#' @param thresholds sample-specific cell type thresholds in .csv format
#' @param save_path path to save celesta assignments of form ("{base_directory}/{sample_name}")
#' @param cofactor asinh cofactor for marker expression; \
#'                 defaults to 10
#'
#' @examples
#' identify_cell_types(
#'   "./521_S1_reg024_cell_measurements.csv",
#'   "./center_negative.csv",
#'   "./521_S1_reg024_thresholds.csv",
#'   "./521_S1_reg024"
#' )
identify_cell_types <- function(quantified_imaging_data, signature_matrix, thresholds, save_path, cofactor = 10) {
    celesta <- CreateCelestaObject(project_title = save_path,
                                   prior_marker_info = signature_matrix,
                                   imaging_data_file = quantified_imaging_data)

    celesta <- FilterCells(celesta)

    celesta <- AssignCells(celesta,
                           high_expression_threshold_anchor = thresholds$ANCHOR,
                           high_expression_threshold_index = thresholds$INDEX,
                           save_result = F)
    
    markers <- colnames(quantified_imaging_data)[!colnames(quantified_imaging_data) %in% c("CELL_IDENTIFIER", "MAJOR_AXIS_LENGTH", "MINOR_AXIS_LENGTH",
                                                                                               "X", "Y", "SIZE", "ECCENTRICITY", "ORIENTATION")]
    transformed_marker_expressions <- asinh(data.matrix(quantified_imaging_data[, markers]) / cofactor)
    marker_probabilities <- CalcMarkerActivationProbability(transformed_marker_expressions)
    colnames(marker_probabilities) <- paste0(colnames(marker_probabilities), "_PROBABILITY")
    
    assignments <- data.frame(celesta@final_cell_type_assignment, check.names = FALSE)
    colnames(assignments) <- gsub(" ", "_", toupper(colnames(assignments)))
    assignments <- assignments[, c("CELL_TYPE_NUMBER", "FINAL_CELL_TYPE")]
    
    assignments$FINAL_CELL_TYPE <- gsub("Unknown", "unknown", assignments$FINAL_CELL_TYPE)
    
    write.csv(cbind(quantified_imaging_data, marker_probabilities, assignments), sprintf("%s_assignments.csv", save_path), row.names = FALSE)
}

#' parses command line arguments
#'
#' @returns a list of parsed arguments
#'
#' @examples
#' parse_arguments()
parse_arguments <- function() {
    parser <- ArgumentParser(description = "batch celesta processing of spatial proteomics files")
    
    parser$add_argument("--data_directory", help = "path to a data directory of .csv files", required = TRUE)
    parser$add_argument("--filter", help = "comma-separated list of sample names to process, or 'all' to process everything; \
                                            defaults to all", default = "all")
    parser$add_argument("--save_path", help = "path to save celesta's output files", required = TRUE)
    parser$add_argument("--signature_matrix", help = "path to the signature matrix", required = TRUE)
    parser$add_argument("--thresholds_directory", help = "path to a thresholds directory of .csv files", required = TRUE)
    
    args <- parser$parse_args()

    return(args)
}

#' reads in arguments and applies celesta to each data file
#'
#' @examples
#' main()
main <- function() {
    args <- parse_arguments()

    data_directory <- args$data_directory
    filter <- args$filter
    save_path <- args$save_path
    signature_matrix <- args$signature_matrix
    thresholds_directory <- args$thresholds_directory

    signature_matrix <- read.csv(signature_matrix)
    sample_quantifications <- list.files(data_directory, pattern = "\\.csv$", full.names = TRUE)
    
    if (filter != "all") {
        filter <- strsplit(filter, ",")[[1]]
        sample_quantifications <- sample_quantifications[basename(sample_quantifications) %in% 
                                                             paste0(filter, "_cell_measurements.csv")]
    }
    
    if (!dir.exists(save_path)) {
        dir.create(save_path, recursive = T)
    }
    
    for (sample in sample_quantifications) {
        sample_data <- read.csv(sample)
        sample_name <- sub("\\_cell_measurements.csv$", "", basename(sample))
        sample_thresholds <- read.csv(file.path(thresholds_directory, paste0(sample_name, "_thresholds.csv")))
        
        identify_cell_types(sample_data, signature_matrix, sample_thresholds, file.path(save_path, sample_name))
    }
}

main()