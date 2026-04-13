library(devtools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(stringr)

load_all("meta-disco")

gather_CLQs <- function(condition) {
    merged_CLQs <- data.frame(A_B = numeric(0))

    for (core in list.files(condition, recursive = FALSE, full.names = TRUE)) {
        core_data <- read.csv(core)
        core_name <- paste(strsplit(basename(core), "_")[[1]][1:3], collapse = "_")

        core_data$sample <- core_name

        names(core_data)[names(core_data) == "X"] <- "x"
        names(core_data)[names(core_data) == "Y"] <- "y"
        names(core_data)[names(core_data) == "FINAL_CELL_TYPE"] <- "cell_type"

        # TODO - comment out! (START)

        # core_data$cell_type <- gsub("Fibroblast \\(FAP\\+CD90-\\)", "Fibroblast (FAP+)", core_data$cell_type)
        # core_data$cell_type <-  gsub("Fibroblast \\(FAP\\+CD90\\+\\)", "Fibroblast (FAP+)", core_data$cell_type)

        # core_data$cell_type <- gsub("CD90\\+ stroma cell", "Stromal cell (CD90+)", core_data$cell_type)
        # core_data$cell_type <- gsub("Regulatory T cell \\(CD4\\+ FOXP3\\+\\)", "Regulatory T cell (CD4+FOXP3+)", core_data$cell_type)
        # core_data$cell_type <- gsub("PDPN\\+ cell", "Mesenchymal cell (podoplanin+)", core_data$cell_type)

        # TODO - comment out! (END)

        spomic_object <- createSpomic(core_data)

        spomic_object <- setSpomicHypers(spomic_object,
                                         tile_size = 250,
                                         k_neigh = 15,
                                         bandwidth = 100,
                                         n_bootstrap = 100,
                                         window_size = 10,
                                         weight_scheme = "linear")

        CLQ = getCLQs(spomic_object)

        # browser()

        # TODO - comment out! (START)

        # CLQ$CLQ <- log(CLQ$CLQ + 1)

        # TODO - comment out! (END)

        names(CLQ)[names(CLQ) == "colocalization_stat"] <- core_name

        merged_CLQs <- merge(merged_CLQs, CLQ, by = "A_B", all = TRUE)
    }

    return(merged_CLQs)
}

gather_cell_counts <- function(condition, classes) {
    columns = c("sample", classes)

    merged_counts <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(merged_counts) <- columns

    for (core in list.files(condition, recursive = FALSE, full.names = TRUE)) {

        # browser()

        core_name <- paste(strsplit(basename(core), "_")[[1]][1:3], collapse = "_")
        core_data <- read.csv(core)

        local_counts <- as.data.frame.matrix(t(table(core_data$FINAL_CELL_TYPE)))
        local_counts$sample <- core_name

        missing_cells_one <- setdiff(names(local_counts), columns)
        missing_cells_two <- setdiff(columns, names(local_counts))

        for (cell in c(missing_cells_one, missing_cells_two)) {
            local_counts[[cell]] <- 0
        }

        stopifnot(length(columns) == ncol(local_counts))

        merged_counts <- rbind(merged_counts, local_counts)
    }

    return(merged_counts)
}

filter_by_counts <- function(CLQs, counts, threshold = 5) {
    for (row_index in 1 : nrow(counts)) {
        row <- counts[row_index, ]

        sample_name <- as.character(row["sample"])
        cell_type_counts <- row[colnames(row) != "sample"]

        drop_cell_types <- names(cell_type_counts)[cell_type_counts <= threshold]

        for (cell_type in drop_cell_types) {
            escaped_cell_type <- gsub("\\(", "\\\\(", cell_type, fixed = FALSE)
            escaped_cell_type <- gsub("\\)", "\\\\)", escaped_cell_type, fixed = FALSE)
            escaped_cell_type <- gsub("\\+", "\\\\+", escaped_cell_type, fixed = FALSE)

            subset_colocalizations <- which(grepl(escaped_cell_type, CLQs$A_B))

            CLQs[subset_colocalizations, sample_name] <- NA
        }
    }

    return(CLQs)
}

filter_by_epithelial <- function(CLQs, counts, threshold = 0.1, direction = "forward") {
    samples <- counts$sample

    counts <- counts[, !names(counts) == "sample"]
    counts$sum <- rowSums(counts)
    counts$epithelial_proportion <- counts$`Epithelial cell (Cytokeratin+)` / counts$sum

    if (direction == "forward") {
        keep_indices <- which(counts$epithelial_proportion <= threshold)
    } else {
        keep_indices <- which(counts$epithelial_proportion >= threshold)
    }

    keep_columns <- c("A_B", samples[keep_indices])

    CLQs <- CLQs[, names(CLQs) %in% keep_columns]
    return(CLQs)
}

generate_volcano_plot <- function(node_negative, node_positive, t_test_results, mode, alpha = 0.05) {
    p_values <- sapply(t_test_results, function(colocation) colocation$p.value)
    # p_values <- p.adjust(p_values, method = "fdr")

    significant_indices <- which(p_values < alpha)

    node_negative_mean <- rowMeans(node_negative[, !(names(node_negative) %in% "A_B")], na.rm = TRUE)
    node_positive_mean <- rowMeans(node_positive[, !(names(node_positive) %in% "A_B")], na.rm = TRUE)

    fold_change <- log2(node_positive_mean / node_negative_mean)

    populations <- str_wrap(node_negative$A_B, width = 20)

    data <- data.frame(fold_change, p_values, populations)

    # TODO - comment out (START)

    # significant_data <- data[significant_indices, ]
    #
    # negative_fold_change <- significant_data[significant_data$fold_change < 0, ]
    # positive_fold_change <- significant_data[significant_data$fold_change > 0, ]
    #
    # negative_fold_change <- negative_fold_change[!grepl("Unknown", negative_fold_change$populations), ]
    # positive_fold_change <- positive_fold_change[!grepl("Unknown", positive_fold_change$populations), ]
    #
    # negative_fold_change <- negative_fold_change[order(negative_fold_change$p_values), ]
    # positive_fold_change <- positive_fold_change[order(positive_fold_change$p_values), ]
    #
    # sampled_negative_fold_change = data.frame()
    # sampled_positive_fold_change = data.frame()
    #
    # CD90_negative_count = 0
    # CD90_positive_count = 0
    #
    # for (i in 1:nrow(negative_fold_change)) {
    #     negative_row <- negative_fold_change[i, ]
    #     positive_row <- positive_fold_change[i, ]
    #
    #     if (CD90_negative_count != 2 && grepl("CD90", negative_row["populations"])) {
    #         CD90_negative_count <- CD90_negative_count + 1
    #         sampled_negative_fold_change <- rbind(sampled_negative_fold_change, negative_row)
    #     }
    #
    #     if (CD90_positive_count != 2 && grepl("CD90", positive_row["populations"])) {
    #         CD90_positive_count <- CD90_positive_count + 1
    #         sampled_positive_fold_change <- rbind(sampled_positive_fold_change, positive_row)
    #     }
    #
    #     if (nrow(sampled_negative_fold_change) != 5 && !grepl("CD90", negative_row["populations"])) {
    #         sampled_negative_fold_change <- rbind(sampled_negative_fold_change, negative_row)
    #     }
    #
    #     if (nrow(sampled_positive_fold_change) != 5 && !grepl("CD90", positive_row["populations"])) {
    #         sampled_positive_fold_change <- rbind(sampled_positive_fold_change, positive_row)
    #     }
    # }
    #
    # # sampled_negative_fold_change <- negative_fold_change[1:5, ]
    # # sampled_positive_fold_change <- positive_fold_change[1:5, ]
    #
    # filtered_data <- rbind(sampled_negative_fold_change, sampled_positive_fold_change)
    # filtered_populations <- str_wrap(filtered_data$populations, width = 20)

    # TODO - comment out (END)

    volcano_plot <- ggplot(data, aes(x = fold_change, y = -log10(p_values))) +

        geom_point(size = 2) +
        geom_point(data = data[significant_indices, ], aes(x = fold_change, y = -log10(p_values)), color = "red", size = 2) +

        geom_label_repel(data = data[significant_indices, ], aes(label = populations),
        # geom_label_repel(data = filtered_data, aes(label = filtered_populations),
                         xlim = c(min(data$fold_change) - 1, max(data$fold_change) + 1),
                         ylim = c(min(-log10(p_values)) - 1, max(-log10(p_values)) + 1),
                         box.padding = 0.8, point.padding = 0.5,
                         segment.size = 0.5, segment.color = "grey50",
                         force = 3,
                         color = "red", show.legend = FALSE,
                         size = 2.5,
                         max.overlaps = Inf) +

        theme_minimal() +
        theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +

        coord_cartesian(clip = "off") +

        labs(x = expression(log[2] ~ "(FoldChange)"), y = expression(-log[10] ~ "(" * italic(P) ~ "value" * ")"), title = "Cell type A _ Cell type B") +

        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.border = element_rect(color = "black", fill = NA),
              axis.ticks = element_line(color = "black"),
              axis.line = element_line(color = "black")) +

        expand_limits(x = c(-2.5, 2.5), y = c(0, 4))

    ggsave(sprintf("../output/%s/volcano_plot.pdf", mode), plot = volcano_plot, width = 8, height = 6)
}

mode <- "lymph_nodes"

if (!dir.exists(sprintf("../output/%s", mode))) {
    dir.create(sprintf("../output/%s", mode), recursive = TRUE)
}

signature_matrix <- read.csv("../../celesta/signature_matrices/uninvolved_lymph_nodes.csv")

node_negative_counts <- gather_cell_counts(sprintf("../input/%s/negative_benign", mode), c(signature_matrix$CELL_TYPE, "unknown"))
node_positive_counts <- gather_cell_counts(sprintf("../input/%s/positive_benign", mode), c(signature_matrix$CELL_TYPE, "unknown"))

node_negative_CLQs <- gather_CLQs(sprintf("../input/%s/negative_benign", mode))
node_positive_CLQs <- gather_CLQs(sprintf("../input/%s/positive_benign", mode))

# TODO - comment out (START)

# filtering strategy 1

# forward pass

# data <- data.frame(CLQ = numeric(), threshold = character(), condition = character(), stringsAsFactors = FALSE)
#
# for (threshold in c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0)) {
#     filtered_negative_CLQs <- filter_by_epithelial(node_negative_CLQs, node_negative_counts, threshold)
#     filtered_positive_CLQs <- filter_by_epithelial(node_positive_CLQs, node_positive_counts, threshold)
#
#     interaction <- "Dendritic cell (CD11c+)_Regulatory T cell (CD4+ FOXP3+)"
#
#     if (!is.vector(filtered_negative_CLQs)) {
#         node_negative_interaction_CLQs <- filtered_negative_CLQs[filtered_negative_CLQs$A_B == interaction, ]
#         node_negative_interaction_CLQs <- node_negative_interaction_CLQs[, !names(node_negative_interaction_CLQs) == "A_B"]
#
#         data <- bind_rows(data, data.frame(CLQ = as.numeric(unlist(node_negative_interaction_CLQs)),
#                                            threshold = as.character(threshold), condition = "N0"))
#     }
#
#     if (!is.vector(filtered_positive_CLQs)) {
#         node_positive_interaction_CLQs <- filtered_positive_CLQs[filtered_positive_CLQs$A_B == interaction, ]
#         node_positive_interaction_CLQs <- node_positive_interaction_CLQs[, !names(node_positive_interaction_CLQs) == "A_B"]
#
#         data <- bind_rows(data, data.frame(CLQ = as.numeric(unlist(node_positive_interaction_CLQs)),
#                                            threshold = as.character(threshold), condition = "N+"))
#     }
# }
#
# boxplot <- ggplot(data, aes(x = threshold, y = CLQ, fill = condition)) +
#     geom_boxplot(width = 0.75) +
#
#     theme(
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1),
#     ) +
#
#     theme(plot.title = element_text(hjust = 0.5)) +
#
#     stat_summary(fun.y = mean, geom = "point", shape = 0, size = 2, color = "black",
#                  position = position_dodge2(width = 0.75,
#                                             preserve = "single")) +
#
#     labs(title = "forward threshold iteration") +
#
#     scale_fill_brewer(palette="RdBu") +
#
#     stat_compare_means(method = "t.test", aes(label = ..p.format..), label.y = 0.15) +
#     stat_compare_means(method = "wilcox.test", aes(label = ..p.format..), label.y = 0.16)
#
# ggsave(sprintf("../output/%s/forward_boxplot.pdf", mode), boxplot, width = 8, height = 5)

# backward pass

# data <- data.frame(CLQ = numeric(), threshold = character(), condition = character(), stringsAsFactors = FALSE)
#
# for (threshold in c(0.5, 0.4, 0.3, 0.2, 0.1, 0.0)) {
#     filtered_negative_CLQs <- filter_by_epithelial(node_negative_CLQs, node_negative_counts, threshold, "backward")
#     filtered_positive_CLQs <- filter_by_epithelial(node_positive_CLQs, node_positive_counts, threshold, "backward")
#
#     interaction <- "Dendritic cell (CD11c+)_Regulatory T cell (CD4+ FOXP3+)"
#
#     if (!is.vector(filtered_negative_CLQs)) {
#         node_negative_interaction_CLQs <- filtered_negative_CLQs[filtered_negative_CLQs$A_B == interaction, ]
#         node_negative_interaction_CLQs <- node_negative_interaction_CLQs[, !names(node_negative_interaction_CLQs) == "A_B"]
#
#         if (length(node_negative_interaction_CLQs) >= 2 &
#             length(node_negative_interaction_CLQs) >= 2) {
#
#             data <- bind_rows(data, data.frame(CLQ = as.numeric(unlist(node_negative_interaction_CLQs)),
#                                                threshold = as.character(threshold), condition = "N0"))
#         }
#     }
#
#     if (!is.vector(filtered_positive_CLQs)) {
#         node_positive_interaction_CLQs <- filtered_positive_CLQs[filtered_positive_CLQs$A_B == interaction, ]
#         node_positive_interaction_CLQs <- node_positive_interaction_CLQs[, !names(node_positive_interaction_CLQs) == "A_B"]
#
#         if (length(node_positive_interaction_CLQs) >= 2 &
#             length(node_positive_interaction_CLQs) >= 2) {
#
#             data <- bind_rows(data, data.frame(CLQ = as.numeric(unlist(node_positive_interaction_CLQs)),
#                                                threshold = as.character(threshold), condition = "N+"))
#         }
#     }
# }
#
# boxplot <- ggplot(data, aes(x = threshold, y = CLQ, fill = condition)) +
#     geom_boxplot(width=0.75) +
#
#     theme(
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1),
#     ) +
#
#     theme(plot.title = element_text(hjust = 0.5)) +
#
#     stat_summary(fun.y = mean, geom = "point", shape = 0, size = 2, color = "black",
#                  position = position_dodge2(width = 0.75,
#                                             preserve = "single")) +
#
#     labs(title = "backward threshold iteration") +
#
#     scale_fill_brewer(palette="RdBu") +
#
#     stat_compare_means(method = "t.test", aes(label = ..p.format..), label.y = 0.15) +
#     stat_compare_means(method = "wilcox.test", aes(label = ..p.format..), label.y = 0.16)
#
# ggsave(sprintf("../output/%s/backwards_boxplot.pdf", mode), boxplot, width = 8, height = 5)

# filtering strategy 2

# node_negative_CLQs <- filter_by_counts(node_negative_CLQs, node_negative_counts)
# node_positive_CLQs <- filter_by_counts(node_positive_CLQs, node_positive_counts)

# TODO - comment out (END)

stopifnot(node_negative_CLQs$A_B == node_positive_CLQs$A_B)

write.csv(node_negative_CLQs, sprintf("../output/%s/node_negative_CLQ.csv", mode))
write.csv(node_positive_CLQs, sprintf("../output/%s/node_positive_CLQ.csv", mode))

write.csv(node_negative_counts, sprintf("../output/%s/node_negative_counts.csv", mode))
write.csv(node_positive_counts, sprintf("../output/%s/node_positive_counts.csv", mode))

t_test_results <- lapply(1:nrow(node_negative_CLQs), function(row) {
    t_test_result <- t.test(node_negative_CLQs[, -which(names(node_negative_CLQs) == "A_B")][row, ],
                            node_positive_CLQs[, -which(names(node_positive_CLQs) == "A_B")][row, ])

    return(t_test_result)
})

p_values <- sapply(t_test_results, function(colocation) colocation$p.value)
# p_values <- p.adjust(p_values, method = "fdr")

p_values <- data.frame(node_negative_CLQs$A_B, p_values)

write.csv(p_values, sprintf("../output/%s/p_values.csv", mode))

generate_volcano_plot(node_negative_CLQs, node_positive_CLQs, t_test_results, mode)
