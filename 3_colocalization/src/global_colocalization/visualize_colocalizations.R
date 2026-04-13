library(rlang)

generate_colocalization_heatmap <- function(data, cell_types, counts_one, counts_two) {
    # counts <- rbind(data.frame(CELL_TYPE = cell_types, GLOBAL_PROPORTION = rowSums(counts_two[, -1]) / sum(as.matrix(counts_two[, -1])), ROW = 2),
    #                 data.frame(CELL_TYPE = cell_types, GLOBAL_PROPORTION = rowSums(counts_one[, -1]) / sum(as.matrix(counts_one[, -1])), ROW = 1))

    # proportions <- ggplot2::ggplot(counts, ggplot2::aes(
    #     x = factor(counts$CELL_TYPE, levels = cell_types),
    #     y = counts$ROW,
    #     fill = counts$GLOBAL_PROPORTION
    # )) +
    #     ggplot2::geom_tile(color = "black", width = 1, height = 1) +
    #     ggplot2::scale_fill_gradient(low = "white", high = "darkgreen", name = "proportion") +
    #     ggplot2::theme_minimal() +
    #     ggplot2::coord_fixed(ratio = 1, ylim = c(0.5, 2.5), xlim = c(0.5, length(cell_types) + 0.5)) +
    #     ggplot2::theme(
    #         axis.title.x = ggplot2::element_blank(),
    #         axis.title.y = ggplot2::element_blank(),
    #         axis.text.x = ggplot2::element_blank(),
    #         axis.text.y = ggplot2::element_blank(),
    #         plot.margin = ggplot2::margin(0, 0, 0, 0)
    #     ) +
    #     ggplot2::scale_y_continuous(
    #         breaks = c(2, 1),
    #         labels = c("condition two", "condition one"),
    #         expand = c(0, 0)
    #     ) +
    #     ggplot2::theme(
    #         panel.grid.major = ggplot2::element_blank(),
    #         panel.grid.minor = ggplot2::element_blank()
    #     )

    all_pairs <- expand.grid(CELL_TYPE_ONE = cell_types, CELL_TYPE_TWO = cell_types) |>
        dplyr::mutate(CELL_TYPE_PAIR = paste(.data$CELL_TYPE_ONE, .data$CELL_TYPE_TWO, sep = "_"))

    heatmap <- all_pairs |> dplyr::left_join(data, by = "CELL_TYPE_PAIR")
    heatmap$FILL_COLOR <- ifelse(is.na(heatmap$FOLD_CHANGE) | is.na(heatmap$P_VALUE) | heatmap$P_VALUE >= 0.05, NA, heatmap$FOLD_CHANGE)

    heatmap <- ggplot2::ggplot(heatmap, ggplot2::aes(x = CELL_TYPE_ONE, y = CELL_TYPE_TWO, fill = FILL_COLOR)) +
        ggplot2::geom_tile(color = "black", width = 1, height = 1) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                      midpoint = 0, na.value = "gray", name = "fold change") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            axis.title = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0)
        ) +

        ggplot2::coord_fixed()

    # combined <- proportions / heatmap + patchwork::plot_layout(heights = c(1, 5))
    # combined

    ggplot2::ggsave("/Users/rohit/Downloads/figures/colocalization_heatmap_negative_positive_wilcox_test_with_correction.pdf", heatmap, width = 10, height = 10, dpi = 300)
    heatmap
}

centers_negative <- read.csv("/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative/colocalizations/global_colocalizations/assignments_global_colocalization_matrix_untransformed_linear_weighting_15_neighbors_100_pixel_bandwidth.csv")
centers_positive <- read.csv("/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/centers/node_positive/colocalizations/global_colocalizations/assignments_global_colocalization_matrix_untransformed_linear_weighting_15_neighbors_100_pixel_bandwidth.csv")
edges_negative <- read.csv("/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/edges/node_negative/colocalizations/global_colocalizations/assignments_global_colocalization_matrix_untransformed_linear_weighting_15_neighbors_100_pixel_bandwidth.csv")
edges_positive <- read.csv("/Users/rohit/Desktop/plevritis_analysis/data/primary_tumors/edges/node_positive/colocalizations/global_colocalizations/assignments_global_colocalization_matrix_untransformed_linear_weighting_15_neighbors_100_pixel_bandwidth.csv")

colocalization_one <- merge(centers_negative, centers_positive, by = "A_B", suffixes = c("_negative", "_positive"))
colocalization_one <- colocalization_one[, !duplicated(names(colocalization_one))]
colocalization_two <- merge(edges_negative, edges_positive, by = "A_B", suffixes = c("_negative", "_positive"))
colocalization_two <- colocalization_two[, !duplicated(names(colocalization_two))]

# colocalization_one <- merge(centers_negative, edges_negative, by = "A_B", suffixes = c("_centers", "_edges"))
# colocalization_one <- colocalization_one[, !duplicated(names(colocalization_one))]
# colocalization_two <- merge(centers_positive, edges_positive, by = "A_B", suffixes = c("_centers", "_edges"))
# colocalization_two <- colocalization_two[, !duplicated(names(colocalization_two))]

# counts_one <- read.csv("/Users/rohit/Downloads/counts/innpd_count_matrix.csv")
# counts_two <- read.csv("/Users/rohit/Downloads/counts/nnpd_count_matrix.csv")

signature_matrix <- read.csv("/Users/rohit/Desktop/plevritis_analysis/celesta/signature_matrices/primary_tumors/centers/node_negative/center_negative_signature_matrix.csv")

colocalization_one <- colocalization_one[!apply(is.na(colocalization_one[, !(names(colocalization_one) == "A_B")]), 1, all), ]
colocalization_two <- colocalization_two[!apply(is.na(colocalization_two[, !(names(colocalization_two) == "A_B")]), 1, all), ]

shared_rows <- intersect(colocalization_one$A_B, colocalization_two$A_B)

colocalization_one <- colocalization_one[colocalization_one$A_B %in% shared_rows, ]
colocalization_two <- colocalization_two[colocalization_two$A_B %in% shared_rows, ]

colocalization_one <- colocalization_one[match(shared_rows, colocalization_one$A_B), ]
colocalization_two <- colocalization_two[match(shared_rows, colocalization_two$A_B), ]

one_data <- colocalization_one[, -which(names(colocalization_one) == "A_B")]
two_data <- colocalization_two[, -which(names(colocalization_two) == "A_B")]

t_test_results <- lapply(seq_along(shared_rows), function(i) {
    x <- as.numeric(one_data[i, ])
    y <- as.numeric(two_data[i, ])

    x <- x[!is.na(x)]
    y <- y[!is.na(y)]

    if (length(x) >= 2 && length(y) >= 2) {
        wilcox.test(x, y)
    } else {
        list(p.value = NA, error = "insufficient data")
    }
})

p_values <- sapply(t_test_results, function(colocation) colocation$p.value)
p_values <- p.adjust(p_values, method = "fdr")

condition_one_mean <- rowMeans(one_data, na.rm = TRUE)
condition_two_mean <- rowMeans(two_data, na.rm = TRUE)

fold_change <- log2((condition_one_mean) / (condition_two_mean))

data <- data.frame(shared_rows, fold_change, p_values)
names(data)[names(data) == "shared_rows"] <- "CELL_TYPE_PAIR"
names(data)[names(data) == "p_values"] <- "P_VALUE"
names(data)[names(data) == "fold_change"] <- "FOLD_CHANGE"

generate_colocalization_heatmap(data, c("unknown", signature_matrix$CELL_TYPE))
