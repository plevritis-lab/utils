library(ggplot2)

mode <- "edge"

node_negative_counts <- read.csv(sprintf("../output/%s/node_negative_counts.csv", mode), check.names = FALSE, row.names = 1)
node_positive_counts <- read.csv(sprintf("../output/%s/node_positive_counts.csv", mode), check.names = FALSE, row.names = 1)

node_negative_counts$nodal_status <- "N0"
node_positive_counts$nodal_status <- "N+"

cumulative_counts <- rbind(node_negative_counts, node_positive_counts)
nodal_status <- cumulative_counts$nodal_status

cumulative_counts <- cumulative_counts[, !names(cumulative_counts) == "sample"]
cumulative_counts <- cumulative_counts[, !names(cumulative_counts) == "nodal_status"]

cumulative_counts$sum <- rowSums(cumulative_counts)
cumulative_counts$epithelial_proportion <- cumulative_counts$`Epithelial cell (Cytokeratin+)` / cumulative_counts$sum

data <- cumulative_counts$epithelial_proportion

breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
midpoints <- (breaks[-length(breaks)] + breaks[-1]) / 2

data_cut <- cut(data, breaks = breaks, right = FALSE)
nodal_cut <- split(nodal_status, data_cut)

data_labels <- data.frame(matrix(0, nrow = length(nodal_cut), ncol = 4))
colnames(data_labels) <- c("centers", "N0", "N+", "y")
rownames(data_labels) <- names(nodal_cut)

data_labels$centers <- midpoints
data_labels$y <- 0

for (interval in names(nodal_cut)) {
    for (label in c("N0", "N+")) {
        data_labels[interval, label] <- sum(nodal_cut[[interval]] == label)
    }
}

data_labels$counts = data_labels$`N+` + data_labels$N0

data <- data.frame(x = breaks, y = rep(0, length(breaks)))

number_line <- ggplot(data, aes(x = x, y = y)) +
    geom_line(size = 1) +
    geom_point(shape = 108, size = 5) +

    geom_text(data = data, aes(label = x), hjust = .5, vjust = 2.5) +
    geom_text(data = data_labels, aes(x = centers, label = `N+`, y = y), hjust = 0.5, vjust = -2) +
    geom_text(data = data_labels, aes(x = centers, label = N0, y = y), hjust = 0.5, vjust = -4) +

    theme_void() +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank())

ggsave(sprintf("../output/%s/epithelial_counts.pdf", mode), plot = number_line, width = 5, height = 3)
