library(devtools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(stringr)

load_all("/scratch/wsspaces/shyunlee-utils/scripts/colocalization/colocalization/global_colocalization/metadisco-main")

compute_global_colocalization <- function(assignments_directory, log_transform = FALSE) {
  composite_colocalizations <- data.frame(A_B = numeric(0))
  
  for (sample in list.files(assignments_directory, full.names = TRUE)) {
    sample_name <- gsub("\\.csv$", "", basename(sample))
    
    sample_assignments <- read.csv(sample)
    sample_assignments$sample <- sample_name
    
    names(sample_assignments)[names(sample_assignments) == "X"] <- "x"
    names(sample_assignments)[names(sample_assignments) == "Y"] <- "y"
    names(sample_assignments)[names(sample_assignments) == "FINAL_CELL_TYPE"] <- "cell_type"
    
    spomic <- createSpomic(sample_assignments)
    spomic <- setSpomicHypers(spomic, tile_size = 250, window_size = 10, k_neigh = 15, 
                              bandwidth = 100, n_bootstrap = 100, weight_scheme = "linear")
    
    global_colocalization <- getCLQs(spomic)
    if (log_transform) {
      global_colocalization$colocalization_stat <- log(global_colocalization$colocalization_stat + 1)
    }
    
    names(global_colocalization)[names(global_colocalization) == "colocalization_stat"] <- sample_name
    composite_colocalizations <- merge(composite_colocalizations, global_colocalization, by = "A_B", all = TRUE)
  }
  
  return(composite_colocalizations)
}

generate_volcano_plot <- function(condition_one, condition_two, t_test_results, base_save_directory, mode, correction = FALSE, clip = FALSE, alpha = 0.05) {
  p_values <- sapply(t_test_results, function(colocation) colocation$p.value)
  if (correction) {
    p_values <- p.adjust(p_values, method = "fdr")
  }
  
  significant_indices <- which(p_values < alpha) 
  
  condition_one_mean <- rowMeans(condition_one[, !(names(condition_one) %in% "A_B")], na.rm = TRUE)
  condition_two_mean <- rowMeans(condition_two[, !(names(condition_two) %in% "A_B")], na.rm = TRUE)
  
  fold_change <- log2((condition_one_mean + 1) / (condition_two_mean + 1))
  
  populations <- str_wrap(condition_one$A_B, width = 20)
  
  data <- data.frame(fold_change, p_values, populations)
  
  if (clip) {
    significant_data <- data[significant_indices, ]
    
    condition_one_significant_data <- significant_data[significant_data$fold_change > 0, ]
    condition_two_significant_data <- significant_data[significant_data$fold_change < 0, ]
    
    condition_one_significant_data <- condition_one_significant_data[!grepl("unknown", condition_one_significant_data$populations), ]
    condition_two_significant_data <- condition_two_significant_data[!grepl("unknown", condition_two_significant_data$populations), ]
    
    condition_one_significant_data <- condition_one_significant_data[order(condition_one_significant_data$p_values), ]
    condition_two_significant_data <- condition_two_significant_data[order(condition_two_significant_data$p_values), ]
    
    sampled_condition_one <- condition_one_significant_data[1:10, ]
    sampled_condition_two <- condition_two_significant_data[1:10, ]
    
    filtered_data <- rbind(sampled_condition_one, sampled_condition_two)
    filtered_populations <- str_wrap(filtered_data$populations, width = 20)
  }
  
  volcano_plot <- ggplot(data, aes(x = fold_change, y = -log10(p_values))) + 
    
    geom_point(size = 2) +
    geom_point(data = data[significant_indices, ], aes(x = fold_change, y = -log10(p_values)), color = "red", size = 2) +
    
    # geom_label_repel(data = data[significant_indices, ], aes(label = populations), # TODO: CHANGE (- CLIP)
    geom_label_repel(data = filtered_data, aes(label = filtered_populations), # TODO: CHANGE (+ CLIP)
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
    
    expand_limits(x = c(-8, 5), y = c(0, 5)) # TODO: CHANGE
  
  ggsave(sprintf("%s/%s/%s_volcano_plot_correction_%d.pdf", base_save_directory, mode, mode, correction), plot = volcano_plot, width = 8, height = 6)
  write.csv(data.frame(condition_one$A_B, p_values, fold_change), sprintf("%s/%s/%s_p_values_correction_%d.csv", base_save_directory, mode, mode, correction))
}

#' reads in arguments and generates pairwise global CLQ values across cell types
#'
#' @examples
#' main()
main <- function() {
  condition_one <- "panin" # TODO: CHANGE
  condition_two <- "pdac" # TODO: CHANGE
  
  base_save_directory <- "../output/TMA1/CLQs_retry" # TODO: CHANGE
  mode <- sprintf("%s_%s", condition_one, condition_two)
  
  if (!dir.exists(sprintf("%s/%s", base_save_directory, mode))) {
    dir.create(sprintf("%s/%s", base_save_directory, mode), recursive = TRUE)
  }
  
  condition_one_assignments_directory <- sprintf("TMA1/%s", condition_one) # TODO: CHANGE
  condition_two_assignments_directory <- sprintf("TMA1/%s", condition_two) # TODO: CHANGE
  
  log_transformation <- FALSE # TODO: CHANGE
  
  condition_one_colocalizations <- compute_global_colocalization(condition_one_assignments_directory, log_transformation) # TODO: CHANGE (+/- LOG TRANSFORMATION)
  condition_two_colocalizations <- compute_global_colocalization(condition_two_assignments_directory, log_transformation) # TODO: CHANGE (+/- LOG TRANSFORMATION)
  
  write.csv(condition_one_colocalizations, sprintf("%s/%s/%s_colocalizations_transformed_%d.csv", base_save_directory, mode, condition_one, log_transformation))
  write.csv(condition_two_colocalizations, sprintf("%s/%s/%s_colocalizations_transformed_%d.csv", base_save_directory, mode, condition_two, log_transformation))
  
  condition_one_colocalizations <- condition_one_colocalizations[rowSums(!is.na(condition_one_colocalizations[, !names(condition_one_colocalizations) %in% "A_B"])) >= 2, ]
  condition_two_colocalizations <- condition_two_colocalizations[rowSums(!is.na(condition_two_colocalizations[, !names(condition_two_colocalizations) %in% "A_B"])) >= 2, ]
  
  shared_rows <- intersect(condition_one_colocalizations$A_B, condition_two_colocalizations$A_B)
  condition_one_colocalizations <- condition_one_colocalizations[condition_one_colocalizations$A_B %in% shared_rows, ]
  condition_two_colocalizations <- condition_two_colocalizations[condition_two_colocalizations$A_B %in% shared_rows, ]
  
  write.csv(condition_one_colocalizations, sprintf("%s/%s/%s_colocalizations_filtered.csv", base_save_directory, mode, condition_one))
  write.csv(condition_two_colocalizations, sprintf("%s/%s/%s_colocalizations_filtered.csv", base_save_directory, mode, condition_two))
  
  t_test_results <- lapply(1:length(shared_rows), function(row) {
    t_test_result <- t.test(condition_one_colocalizations[, -which(names(condition_one_colocalizations) == "A_B")][row, ],
                            condition_two_colocalizations[, -which(names(condition_two_colocalizations) == "A_B")][row, ])
    
    return(t_test_result)
  })

  generate_volcano_plot(condition_one_colocalizations, condition_two_colocalizations, t_test_results, base_save_directory, mode, clip = TRUE) # TODO: CHANGE (+/- CORRECTION, +/- CLIP)
}

main()