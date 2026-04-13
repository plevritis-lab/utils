#' Tile a 2D image
#'
#' @description
#' Splits a 2D dataset with `x` and `y` coordinates into tiles of specified size
#' and rescales the coordinates within each tile to the range [0, 1].
#'
#' @param img A data frame with numeric `x` and `y` columns representing coordinates.
#' @param tile_size A positive numeric value specifying the size of each tile.
#' @return A list of tibbles, where each tibble represents a tile with rescaled coordinates.
#' @details
#' The function calculates tile row and column indices based on the `tile_size`,
#' assigns each data point to a corresponding tile, and rescales te coordinates
#' within each tile  to a normalized scale between 0 and 1.
#'
#' @importFrom stringr str_pad
#' @importFrom dplyr mutate
#' @export
tileImage <- function(img, tile_size) {
  # Input validation
  if (!is.data.frame(img)) stop("`img` must be a data frame.")
  if (!all(c("x", "y") %in% names(img))) stop("`img` must contain `x` and `y` columns.")
  if (!is.numeric(tile_size) || tile_size <= 0) stop("`tile_size` must be a positive numeric value.")

  # Calculate number of tiles in each dimension
  n_tile_x <- ceiling(diff(range(img$x)) / tile_size)
  n_tile_y <- ceiling(diff(range(img$y)) / tile_size)

  # Define steps (tile width and height) can probably remove
  x_step <- tile_size
  y_step <- tile_size

  # Split the image into tiles and rescale the coordinates within each tile
  tiled_image <- img |>
    dplyr::mutate(
      # Calculate row and column positions based on tile size
      grid_row = as.numeric(cut(x, breaks = seq(min(x), max(x), length.out = n_tile_x + 1), include.lowest = TRUE, labels = 1:n_tile_x)),
      grid_col = as.numeric(cut(y, breaks = seq(min(y), max(y), length.out = n_tile_y + 1), include.lowest = TRUE, labels = 1:n_tile_y)),
      grid = paste0(
        stringr::str_pad(grid_row, width = nchar(as.character(n_tile_x)), pad = "0"),
        "_",
        stringr::str_pad(grid_col, width = nchar(as.character(n_tile_y)), pad = "0")
      )
    ) |>
  dplyr::mutate(
    x_rescaled = (x - (min(x) + (grid_row - 1) * x_step)) / x_step,
    y_rescaled = (y - (min(y) + (grid_col - 1) * y_step)) / y_step
  )

  # Split into list of tiles
  tiled_list <- split(tiled_image, tiled_image$grid)

  return(tiled_list)
}

#' Filter tiled data
#'
#' @description
#' Filters a list of tibbles based on the presence or absence of `NA` values.
#'
#' @param tibbles Filters a list of tibbles based on the presence or absence of `NA` values.
#' @param find_empty Logical. If `TRUE`, returns only tiles that contain `NA` values.
#' If `FALSE`, returns only tiles without any `NA` values. Default is `FALSE`.
#' @return A filtered list of tibbles based on the specified condition.
#'
#' @export
filterTiles <- function(tibbles, find_empty=FALSE) {
  # Ensure input is a list
  if (!is.list(tibbles)) stop("`tibbles` must be a list of tibbles.")

  # Define a function to check for NA values in a tibble
  has_na <- function(tbl) any(is.na(tbl))

  # Filter based on `find_empty`
  filtered_tibbles <- if (find_empty) {
    tibbles[sapply(tibbles, has_na)]
  } else {
    tibbles[!sapply(tibbles, has_na)]
  }

  return(filtered_tibbles)
}

#' Perform Window Bootstrap
#'
#' @@description
#' This function performs a spatial bootstrapping operation on a `Spomic` object by dividing the data into tiles,
#' applying a moving wondow of specified tile radius, and sampling within the non-empty tiles. The function resamples
#' tiles within the specified window size and constructs a bootstrapped image.
#'
#' @param spomic_obj A `Spomic` object containing the spatial data and tiling parameters set using `metadisco::setSpomicHypers()`
#' @return A new bootstrapped `Spomic` object where the data has been resampled.
#' @details
#' The function performs the following steps:
#' 1. Tiles the data into smaller chunks based on the tilesize parameter.
#' 2. Identifies empty and non-empty tiles.
#' 3. For each non-empty tile, selects neighboring tiles within the specified `window_size`.
#' 4. Samples tile from the neighboring non-empty tiles.
#' 5. Generate new image and returns a new bootstrapped `Spomic` object.
#'
#' @export
windowBootstrap <- function(spomic_obj) {
  # Input validation
  stopifnot(inherits(spomic_obj, "Spomic"))

  # Tiling and parameters
  x_step <- spomic_obj@tile_size
  y_step <- spomic_obj@tile_size
  tiles <- tileImage(spomic_obj@df, tile_size = spomic_obj@tile_size)

  # Separate empty and non-empty tiles
  empty_tiles <- filterTiles(tiles, find_empty=TRUE)
  non_empty_tiles <- filterTiles(tiles, find_empty=FALSE)

  # Prepare for bootstrapping
  bootstrapped_tiles <- list()
  tile_names_tibble <- dplyr::tibble(names(tiles)) |>
    tidyr::separate('names(tiles)', into = c("row", "col"), sep = "_") |>
    dplyr::mutate(row = as.integer(row), col = as.integer(col))

  # Iterate over non-empty tiles
  for(tile in seq_along(non_empty_tiles)) {
    current_tile <- non_empty_tiles[tile]
    indx_a <- as.numeric(strsplit(names(current_tile), "_")[[1]][1])
    indx_b <- as.numeric(strsplit(names(current_tile), "_")[[1]][2])

    # Define window range
    combos <- tidyr::crossing(
      row = (indx_a - spomic_obj@window_size):(indx_a + spomic_obj@window_size),
      col = (indx_b - spomic_obj@window_size):(indx_b + spomic_obj@window_size)
    )

    # Identify targets within the window
    targets <- which(tile_names_tibble$row %in% combos$row & tile_names_tibble$col %in% combos$col)

    # Filter tiles within the window
    window_tiles <- tiles[targets]
    # window_empty_tiles <- filterTiles(window_tiles, find_empty = TRUE)
    window_non_empty_tiles <- filterTiles(window_tiles, find_empty = FALSE)

    # Randomly sample a tile
    shuffled_tiles <- sample(x = window_non_empty_tiles, size = 1, replace = TRUE)
    names(shuffled_tiles) <- names(non_empty_tiles[tile])

    # Append to bootstrapped tiles
    bootstrapped_tiles <- append(bootstrapped_tiles, shuffled_tiles)
  }

  # Combine bootstrapped tiles with empty tiles
  bootstrapped_tiles <- bootstrapped_tiles[unique(names(bootstrapped_tiles))]
  bootstrapped_tileImage <- append(bootstrapped_tiles, empty_tiles)

  # Rescale the tile coordinates to reconstruct the full image
  for (k in seq_along(bootstrapped_tileImage)) {
    idx <- strsplit(names(bootstrapped_tileImage)[k], "_")[[1]]
    i <- as.numeric(idx[1])
    j <- as.numeric(idx[2])
    bootstrapped_tileImage[[k]] <- bootstrapped_tileImage[[k]] |>
      # Rescale the tile positions to create a cohesive image
      dplyr::mutate(x = (x_rescaled*x_step + ((i-1) * x_step)) + min(spomic_obj@coords$x),
                    y = (y_rescaled*y_step + ((j-1) * y_step)) + min(spomic_obj@coords$y))
  }

  # Combine all tiles into a single data frame
  bootstrapped_entireImage <- dplyr::bind_rows(bootstrapped_tileImage, .id = "tile") |>
    tidyr::drop_na()

  # TODO: Can fix the return to automatically inherit the hyperparameters from the original spomic,
  # since in other parts of the code you have to manually respecify them which is annoying.

  # Return a new Spomic object
  return(createSpomic(bootstrapped_entireImage))
}


# ==============================================================================
# Deprecated
# ==============================================================================
# No longer used. Retaining for documentation.
# blockBootstrap <- function(spomic_obj){
# x_step <- diff(range(spomic_obj@df$x)) / spomic_obj@n_tile  # y_step <- diff(range(spomic_obj@df$y)) / spomic_obj@n_tile
#
#   n_tile_x <- ceiling(diff(range(img$x)) / tile_size)
#   n_tile_y <- ceiling(diff(range(img$y)) / tile_size)
#
#   tiles <- tileImage(spomic_obj@df, n_tile = spomic_obj@n_tile)
#   empty_tiles <- filterTiles(tiles, find_empty = TRUE)
#   non_empty_tiles <- filterTiles(tiles, find_empty = FALSE)
#
#   # Sample the non-empty tiles (with replacement)
#   shuffled_tiles <- sample(x = non_empty_tiles,
#                            size = length(non_empty_tiles),
#                            replace = TRUE)
#   # Overwrite their position names
#   # This basically assigns the resampled tiles to their new locations
#   names(shuffled_tiles) <- names(non_empty_tiles)
#
#   for(k in 1:length(shuffled_tiles)){
#     idx <- strsplit(names(shuffled_tiles)[k], "_")[[1]]
#     i <- as.numeric(idx[1])
#     j <- as.numeric(idx[2])
#     shuffled_tiles[k][[1]] <- shuffled_tiles[k][[1]] |>
#       # Rescale the tile positions to create a cohesive image
#       dplyr::mutate(x = (x_rescaled*x_step + ((i-1) * x_step)) + min(spomic_obj@coords$x),
#                     y = (y_rescaled*y_step + ((j-1) * y_step)) + min(spomic_obj@coords$y))
#   }
#   # Concatenate with the nonempty tiles to get all tiles in a full image
#   # bootstrapped_image <- c(shuffled_tiles, empty_tiles)
#   bootstrapped_image <- shuffled_tiles[order(names(shuffled_tiles))]
#
#   bootstrapped_image <- dplyr::bind_rows(bootstrapped_image, .id="tile")
#   return(createSpomic(bootstrapped_image))
# }
