#' Spomic (Spatial omic) Class
#'
#' @description
#' Defines the structure for Spomic objects, used to store spatial and colocalization data.
#'
#' @slot sample_id A string with the sample name / ID
#' @slot df A data frame containing columns for coordinates, sample, and cell type.
#' @slot coords A data frame specifying spatial coordinates (TODO: in the future, this can probably be removed...)
#' @slot cell_assignments A vector assigning cells to groups or types.
#' @slot tile_size Numeric, tile size to be used in spatial bootstrapping.
#' @slot window_size Integer, window size for spatial analysis.
#' @slot k_neigh Integer, number of nearest neighbors to consider for colocalization calculation
#' @slot bandwidth Integer, bandwidth for smoothing operations.
#' @slot n_bootstrap Integer, number of spatial bootstrapping iterations to perform.
#' @slot nb_list A list of indices for the k_neigh nearest neighbors of each cell.
#' @slot colocalization A data frame with colocalization data processed for Gaussianity and bias.
#' @slot weight_scheme A character string describing the weighting scheme used in colocalization calculation.
#'
#' @export
setClass(
  "Spomic",
  slots = list(
    sample_id = "character",
    df = "data.frame",
    coords = "data.frame",
    cell_assignments = "vector",
    tile_size = "numeric",
    window_size = "integer",
    k_neigh = "integer",
    bandwidth = "integer",
    n_bootstrap = "integer",
    nb_list = "list",
    colocalization = "data.frame",
    processed_colocalization = "data.frame",
    weight_scheme = "character"
  )
)

# Helper function to validate columns of the input data frame needed to initiate the Spomic class
validate_columns <- function(df, required_columns) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_columns, collapse = ", ")))
  }
}

#' Create a Spomic object
#'
#' @description
#' Constructs a 'Spomic' object from either a data frame or a path to a CSV file.
#'
#'
#' @param p A data frame or a path to a CSV file.
#' @param drop_na A boolean for whether to remove rows containing NAs. Defaults to TRUE.
#' @return An object of class `Spomic`.
#' @details
#' The input data must include the following columns:
#' - `x`: Numeric values representing sptial x-coordinates
#' - `y`: Numeric values representing spatial y-coordinates.
#' - `cell_type`: A character vector denoting cell type assignments.
#' - `sample`: A character vector specifying the sample name / ID.
#'
#' @importFrom tidyr drop_na
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   x = runif(10),
#'   y = runif(10),
#'   cell_type = sample(c("A", "B", "C"), 10, replace = TRUE),
#'   sample = "Sample_1"
#' )
#' spomic_obj <- createSpomic(df)
#'
#' # Example with a CSV file
#' # write.csv(df, "data.csv", row.names = FALSE)
#' # spomic_obj <- createSpomic("data.csv")
#'
#' @export
createSpomic <- function(p, drop_na = TRUE) {
  # Check input type and read data
  if (is.character(p)) {
    if (!file.exists(p)) stop("File does not exist: ", p)
    df <- read.csv(p)
  } else if (is.data.frame(p)) {
    df <- p
  } else {
    stop("Input must be a data.frame or path to a CSV file.")
  }

  if (drop_na) {
    df <- tidyr::drop_na(df)
  }

  if (nrow(df) == 0) stop("Data frame is empty.")

  # Validate required columns
  required_columns <- c("x", "y", "cell_type", "sample")
  validate_columns(df, required_columns)

  # Ensure correct column types
  if (!is.numeric(df$x) || !is.numeric(df$y)) {
    stop("Columns 'x' and 'y' must contain numeric values.")
  }

  # Ensure `sample` column has a single unique value
  unique_samples <- unique(df$sample)
  if (length(unique_samples) > 1) {
    stop("Data frame contains multiple unique sample IDs. Ensure only one sample per Spomic object.")
  }

  # Create Spomic object
  object <- new("Spomic",
                sample_id = as.character(unique_samples),
                df = df,
                coords = df[, c("x", "y"), drop = FALSE],
                cell_assignments = as.character(df$cell_type),
                tile_size = 0,  # Default or configurable value
                window_size = 0L,  # Default or configurable value
                k_neigh = 0L,  # Default or configurable value
                bandwidth = 0L,  # Default or configurable value
                n_bootstrap = 0L,  # Default or configurable value
                nb_list = list(),  # Default empty list
                colocalization = data.frame(),  # Default empty data frame
                processed_colocalization = data.frame(),  # Default empty data frame
                weight_scheme = ""
  )

  return(object)
}

#' Set hyperparameters for a Spomic object
#'
#' @param spomic_obj An object of class Spomic.
#' @param tile_size Tile size.
#' @param window_size Window size.
#' @param k_neigh Number of neighbors.
#' @param bandwidth Bandwidth parameter.
#' @param n_bootstrap Number of bootstrap samples.
#' @param weight_scheme Weight scheme for colocation quotient calculation ("constant", "linear", "squared_exponential")
#' @param precompute_neighbors Boolean indicating whether to precompute neighbors.
#' @return Updated Spomic object.
#' @export
setSpomicHypers <- function(spomic_obj,
                            tile_size,
                            window_size,
                            k_neigh,
                            bandwidth,
                            n_bootstrap,
                            weight_scheme,
                            precompute_neighbors = TRUE) {
  # Validate object type
  if (!inherits(spomic_obj, "Spomic")) stop("Input must be a Spomic object.")

  # Validate input parameters
  if (any(c(tile_size, window_size, k_neigh, bandwidth, n_bootstrap) <= 0)) {
    stop("Numerical parameters must be positive.")
  }
  if (!is.character(weight_scheme) || weight_scheme == "") {
    stop("Weight scheme must be a non-empty character string.")
  }

  spomic_obj@tile_size <- as.numeric(tile_size)
  spomic_obj@window_size <- as.integer(window_size)
  spomic_obj@k_neigh <- as.integer(k_neigh)
  spomic_obj@bandwidth <- as.integer(bandwidth)
  spomic_obj@n_bootstrap <- as.integer(n_bootstrap)
  spomic_obj@weight_scheme <- as.character(weight_scheme)

  # Optionally precompute neighbors
  if (precompute_neighbors) {
    spomic_obj@nb_list <- getNeighbors(spomic_obj)
  } else {
    spomic_obj@nb_list <- list() # Initialize to an empty list if not precomputing
  }

  return(spomic_obj)
}


#' Save a Spomic object
#'
#' @description
#' Saves a `Spomic` object as an RDS file in a specified directory. The RDS will save in a subdirectory containing hyperparameter information.
#'
#' @param spomic_obj An object of class `Spomic`.
#' @param path A character string specifying the directory where the object should be saved.
#' @details
#' The function constructs a directory path using metadata from the `Spomic` object,
#' including parameters such as `tilesize`, `nbors` (neighbors), `bandwidth`,
#' `bootstrap`, and `window`. The directory and file are created dynamically if
#' they do not exist.
#'
#' The final file is saved with the naming convention `<sample_id>_spomic.Rds`.
#'
#' @importFrom tools file_path_sans_ext
#' @export
saveSpomic <- function(spomic_obj, path){
  # Validate inputs
  stopifnot(inherits(spomic_obj, "Spomic"))
  if (!dir.exists(path)) {
    stop("The specified directory does not exist: ", path)
  }

  # Construct dynamic subdirectory based on Spomic metadata
  subdir <- file.path(
    path,
    paste0(
      "tilesize", spomic_obj@tile_size,
      "_nbors", spomic_obj@k_neigh,
      "_bandwidth", spomic_obj@bandwidth,
      "_bootstrap", spomic_obj@n_bootstrap,
      "_window", spomic_obj@window_size
    ),
    paste0(spomic_obj@weight_scheme, "_weights")
  )

  # Create directory if it does not exist
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
  }

  # Save the Spomic object as an RDS file
  rds_path <- file.path(subdir, paste0(spomic_obj@sample_id, "_spomic.Rds"))
  saveRDS(spomic_obj, rds_path)
  message("Spomic object saved to: ", rds_path)
}
