library(ggplot2)
library(mclust)

if (!dir.exists("../output/TMA1/innpd_colocalization_histograms")) {
  dir.create("../output/TMA1/innpd_colocalization_histograms", recursive = TRUE)
}

# below is reading colocation that isn't filtered
global_colocalizations <- read.csv("../output/TMA1/CLQs_retry/innpd_nnpd/innpd_colocalizations_transformed_0.csv", check.names = FALSE, row.names = 1)

set.seed(123)

# iterate through rows of global colocalization (cell pairs i.e. 900 in this case). 
for (i in 1:nrow(global_colocalizations)) {
  row <- data.frame(global_colocalizations[i, ]) # getting cell pairs
  
  cell_pair_label <- row$A_B #getting labels of cell pairs
  metrics <- row[, -which(names(row) == "A_B")] # getting CLQs across sample
  
  # rowSums(is.na(metrics))[[1]] <= length(metrics) - 5 && 
  if (length(unique(na.omit(unlist(metrics)))) > 1) { #looking at CLQs and checking at least 2 unique values
    fit <- Mclust(na.omit(unlist(metrics))) # fitting how many Gaussians are in my data and fitting it
    
    metrics_sequence <- data.frame( #order sequence as per x values
      x = seq(min(metrics, na.rm = TRUE), 
              max(metrics, na.rm = TRUE), 
              length.out = 1000)
    )
    
    density <- matrix(0, nrow = nrow(metrics_sequence), ncol = fit$G) # x by y plotting
    
    for (i in 1:fit$G) { #calculate Gaussian density acc to formula - Gaussian distribution formula
      density[, i] <- dnorm(unlist(metrics_sequence), mean = fit$parameters$mean[i], sd = sqrt(fit$parameters$variance$sigmasq[i]))
    }
    
    mixture_density <- rowSums(density * fit$parameters$pro)
    plot_data <- data.frame(x = metrics_sequence$x, density = mixture_density) # x values from above y values from all Gaussians - making df
    
    #below does formatting
    ymax <- tryCatch({ 
      max(density(na.omit(unlist(metrics)))$y)
    }, error = function(e) {
      1
    })
    
    #sturges formula to calculate optimal bin for histogram
    sturges_bins <- ceiling(log2(length(metrics)) + 1)
    # manual binning below
    # sturges_bins <- 20
    histogram <- ggplot() +
      
      geom_histogram(aes(x = na.omit(unlist(metrics)), 
                         y = after_stat(density)),
                     color = "black",
                     fill = "skyblue",
                     alpha = 0.7,
                     bins = sturges_bins) +
      
      geom_density(data = data.frame(x = na.omit(unlist(metrics))), aes(x = x), alpha = 0.2, fill = "blue") +
      
      # geom_line(data = plot_data, aes(x = x, y = density), color = "red", alpha = 0.2) + # toggle plot controlling number of distributions
      
      labs(title = cell_pair_label, x = "colocalization", y = "density") +
      
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
      
      ylim(0, ymax + 0.25) +
      
      annotate("text", x = Inf, y = Inf, label = sprintf("G = %d", fit$G), 
               hjust = 1.2, vjust = 1.5, size = 3,
               color = "red", fontface = "italic")
    
    ggsave(sprintf("../output/TMA1/innpd_colocalization_histograms/innpd_%s.png", cell_pair_label), plot = histogram, width = 6, height = 4)
  }
}