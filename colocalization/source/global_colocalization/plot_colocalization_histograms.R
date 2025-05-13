library(ggplot2)
library(mclust)

mode <- "lymph_nodes"
condition <- "node_positive"

CLQs <- read.csv(sprintf("../output/%s/%s_CLQ.csv", mode, condition), check.names = FALSE, row.names = 1)

dir.create(sprintf("../output/%s/%s_histograms", mode, condition), showWarnings = FALSE)

set.seed(123)

for (i in 1:nrow(CLQs)) {
    row <- data.frame(CLQs[i, ])
    label <- row$A_B
    
    metrics <- row[, -which(names(row) == "A_B")]
    
    fit <- Mclust(na.omit(unlist(metrics)))
    
    metrics_seq <- data.frame(
        x = seq(min(metrics, na.rm = TRUE), max(metrics, na.rm = TRUE), length.out = 1000)
    )
    
    density <- matrix(0, nrow = nrow(metrics_seq), ncol = fit$G)
    
    for (i in 1:fit$G) {
        density[, i] <- dnorm(unlist(metrics_seq), mean = fit$parameters$mean[i], sd = sqrt(fit$parameters$variance$sigmasq[i]))
    }
    
    mixture_density <- rowSums(density * fit$parameters$pro)
    plot_data <- data.frame(x = metrics_seq, density = mixture_density)
    
    sturges_bins <- ceiling(log2(length(metrics)) + 1)
    histogram <- ggplot() +

        geom_histogram(aes(x = na.omit(unlist(metrics)), 
                           y = after_stat(density)),
                       color = "black",
                       fill = "skyblue",
                       alpha = 0.7,
                       bins = sturges_bins) +

        geom_density(data = na.omit(t(metrics)), aes(x = na.omit(t(metrics))), alpha = 0.2, fill = "blue") +
        
        # geom_line(data = plot_data, aes(x = x, y = density), color = "red", alpha = 0.2) + # toggle plot controlling number of distributions

        labs(title = label, x = "colocalization", y = "density") +

        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),

            plot.title = element_text(size = 9, hjust = 0.5),
            axis.title.x = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            axis.text.x = element_text(size = 5),
            axis.text.y = element_text(size = 5)
        ) +
        
        ylim(0, max(density(na.omit(unlist(metrics)))$y) + 0.25) +
        
        annotate("text", x = Inf, y = Inf, label = sprintf("G = %d", fit$G), 
                    hjust = 1.2, vjust = 1.5, size = 3,
                    color = "red", fontface = "italic")

    ggsave(sprintf("../output/%s/%s_histograms/%s.png", mode, condition, label), plot = histogram, width = 6, height = 4)
}