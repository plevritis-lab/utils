make_dir <- function(dir){
  if(!dir.exists(dir)){
    dir.create(dir, showWarnings = TRUE)
    print(paste("Created directory: ", dir))
  }
}

#' @export
plotCellProportions <- function(spomic_obj, stack = TRUE) {
  proportions <- spomic_obj@df |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(count = dplyr::n()) |>
    dplyr::mutate(proportion = count / sum(count),
                  label = paste0(round(proportion * 100, 1), "%"),
                  cell_type = forcats::fct_reorder(cell_type, proportion, .desc = TRUE))

  if(stack) {
    p <- ggplot(proportions, aes(x = "", y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", width = 0.25) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Proportions of Cell Types", x = NULL, y = "Proportion") +
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank())
  }

  if(!stack) {
    p <- ggplot(proportions, aes(x = reorder(cell_type, proportion), y = proportion)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(x = "Cell Type", y = "Proportion", title = "Cell Type Proportion") +
        theme_minimal()
  }

  return(p)
}


#' @export
plotSpomic <- function(spomic_obj, point_size = 0.75){
  g <- ggplot(spomic_obj@df, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(color = "white"),
          legend.spacing.y = unit(0.1, 'cm'),    # Reduce vertical spacing
          legend.key.size = unit(0.5, 'cm'),       # Increase color legend size
          legend.background = element_rect(fill = "black", color = NA),  # Black legend background
          legend.box.background = element_rect(fill = "black", color = NA)  # Black legend box background
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))  # Control color legend point size

  return(g)
}

#' @export
plotCellPair <- function(spomic_obj, cellA, cellB, point_size = 0.75){
  df <- spomic_obj@df |>
    dplyr::mutate(cells = dplyr::case_when(cell_type == !!cellA | cell_type == !!cellB ~ as.character(cell_type),
                                           TRUE ~ "Other"))

  unique_cells <- unique(df$cells)
  default_colors <- c("red", "blue")  # Default colors for cellA and cellB
  other_color <- "gray50"

  # Constructing color mapping dynamically
  color_values <- c(setNames(default_colors, c(cellA, cellB)), Other = other_color)

  g <- ggplot(df, aes(x = x, y = y, color = cells)) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    scale_color_manual(values = color_values) +
    theme(plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(color = "white"),
          legend.spacing.y = unit(0.1, 'cm'),    # Reduce vertical spacing
          legend.key.size = unit(0.5, 'cm'),       # Increase color legend size
          legend.background = element_rect(fill = "black", color = NA),  # Black legend background
          legend.box.background = element_rect(fill = "black", color = NA)  # Black legend box background
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))
  return(g)
}

#' @export
plotKDE <- function(spomic, cell) {
  require(viridis)
  full_data <- spomic@df
  # Compute convex hull indices
  hull_indices <- chull(full_data$x, full_data$y)
  # Extract the convex hull points in order
  hull_data <- full_data[hull_indices, ]

  # Filtered dataset for density estimation
  filtered_data <- spomic@df |>
    dplyr::filter(cell_type == cell)

  # Create the plot with convex hull outline
  p <- ggplot() +
    coord_fixed(ratio = 1) +
    geom_density_2d_filled(data = filtered_data, aes(x, y, fill = as.numeric(after_stat(level)))) +
    scale_fill_viridis_c(option = "magma") +  # Smooth gradient color
    geom_polygon(data = hull_data, aes(x, y), fill = NA, color = "white", linewidth = 1) +

    # Title & Theme
    ggtitle(paste(cell, "Kernel Density with Convex Hull Outline")) +
    ggpubr::theme_pubclean() +
    labs(fill = "Kernel Density")

  return(p)
}


 #' @export
require(ggrepel)
# require(ggtext)
plotVolcano <- function(metadisco_obj, condition1_name, condition2_name){
  collapsed_tibble <- metadisco_obj@results |>
    dplyr::arrange(desc(abs(zscore))) |>
    tidyr::separate(A_B, into = c("A", "B"), sep = "_", remove = FALSE) |>
    dplyr::mutate(differential = dplyr::case_when(
      zscore > 0 & FDR < 0.05 ~ "Up",
      zscore < 0 & FDR < 0.05 ~ "Down",
      TRUE ~ "Insignificant"),
      bar_color = dplyr::case_when(log2fc <= 0 ~ 0,
                                   TRUE ~ 1))

  significant_indices <- which(collapsed_tibble$FDR < 0.05)
  volcano_plot <- ggplot(collapsed_tibble, aes(x = log2fc, y = -log10(pval))) +

    geom_point(size = 2) +
    geom_point(data = collapsed_tibble[significant_indices, ], aes(x = log2fc, y = -log10(pval), color = differential), size = 3) +
    scale_color_manual(values = c("Down" = "darkblue", "Up" = "darkred", "Insignificant" = "gray"), guide = FALSE) +
    geom_label_repel(data = collapsed_tibble[significant_indices, ],
                     aes(label = A_B, color = differential),
                     min.segment.length = unit(1, "lines"),
                     xlim = c(-3, 3), ylim = c(0, 6),
                     box.padding = 0.5, point.padding = 0.5,
                     segment.size = 0.1, segment.color = "grey50",
                     force = 3,
                     show.legend = FALSE,
                     size = 3.0,
                     max.overlaps = Inf) +
    labs(title = paste0("<span style='color:darkblue;'>",
                        condition1_name,
                        "</span> vs. <span style='color:darkred;'>",
                        condition2_name,
                        "</span>"),
         y = "-log10(adj-pval)") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    # theme(plot.title = ggtext::element_markdown()) +
    # theme(text=element_text(size=12), #change font size of all text
    #       axis.text=element_text(size=8), #change font size of axis text
    #       axis.title=element_text(size=10), #change font size of axis titles
    #       plot.title=element_text(size=10), #change font size of plot title
    #       legend.text=element_text(size=8), #change font size of legend text
    #       legend.title=element_text(size=10))
  theme(plot.margin = margin(2, 2, 2, 2, "cm"))

  return(volcano_plot)
}

#' @export
plotMetaBarplot <- function(metadisco_obj, condition1_name, condition2_name, barwidth = 0.1, drop_insignificant = TRUE){
  require(ggtext)
  collapsed_tibble <- metadisco_obj@results |>
    dplyr::arrange(desc(abs(zscore))) |>
    tidyr::separate(A_B, into = c("A", "B"), sep = "_", remove = FALSE) |>
    dplyr::mutate(differential = dplyr::case_when(
      zscore > 0 & FDR < 0.05 ~ "Up",
      zscore < 0 & FDR < 0.05 ~ "Down",
      TRUE ~ "Insignificant"),
      bar_color = dplyr::case_when(log2fc <= 0 ~ 0,
                                   TRUE ~ 1))

  if(drop_insignificant){
    collapsed_tibble <- collapsed_tibble |>
      dplyr::filter(differential != "Insignificant") |>
      dplyr::filter(FDR < 0.05)
  }

  # barplot <- collapsed_tibble |>
  #   dplyr::arrange(log2fc) |>
  #   ggplot(aes(x = log2fc, y = reorder(A_B, -log2fc), fill=differential)) +
  #   geom_bar(stat="identity", width = barwidth) +
  #   scale_fill_manual(values = c("Down" = "darkblue", "Up" = "darkred", "Insignificant" = "gray"), guide = FALSE) +
  #   theme_minimal() +
  #   theme(plot.title = element_markdown()) + # Use element_markdown to enable HTML styling
  #   # theme(text=element_text(size=10), #change font size of all text
  #   #       axis.text=element_text(size=8), #change font size of axis text
  #   #       axis.title=element_text(size=10), #change font size of axis titles
  #   #       legend.text=element_text(size=8), #change font size of legend text
  #   #       legend.title=element_text(size=10)) +
  #   labs(title = paste0("<span style='color:darkblue;'>",
  #                       condition1_name,
  #                       "</span> vs. <span style='color:darkred;'>",
  #                       condition2_name,
  #                       "</span>"),
  #        y = "Cell pair")

  barplot <- collapsed_tibble |>
    dplyr::arrange(log2fc) |>
    ggplot(aes(x = log2fc, y = reorder(A_B, -log2fc), fill=differential)) +
    geom_bar(stat="identity", width = barwidth) +
    scale_fill_manual(name = "Significant colocalizations",
                      values = c("Down" = "darkblue", "Up" = "darkred", "Insignificant" = "gray"),
                      labels = c("Down" = condition1_name, "Up" = condition2_name, "Insignificant" = "Not significant")) +
    labs(x = expression(log[2] ~ "FC"),
         y = "Cell type pairs") +
    ggpubr::theme_pubr()
  return(barplot)
}

plotNetwork <- function(metadisco_obj, keep_significant = TRUE, display_weights = FALSE){
  require(qgraph)
  adjacency_mat <- metadisco_obj@results |>
    tidyr::separate(A_B, into = c("A", "B"), sep = "_") |>
    dplyr::select(A, B, zscore) |>
    tidyr::pivot_wider(names_from = B, values_from = zscore, names_repair = "minimal")

  adjacency_mat <- as.matrix(adjacency_mat[-1]) # Remove the first column if it's the row names
  rownames(adjacency_mat) <- colnames(adjacency_mat)
  if(keep_significant) adjacency_mat[abs(adjacency_mat) < 1.96] <- 0
  if(display_weights) return(qgraph(abs(adjacency_mat), labels = colnames(adjacency_mat), label.norm = "OOOO", edge.labels = TRUE, posCol = "darkred", negCol = "darkblue"))
  else return(qgraph(adjacency_mat, labels = colnames(adjacency_mat), label.norm = "OOOO", edge.labels = FALSE, posCol = "darkred", negCol = "darkblue"))
}
