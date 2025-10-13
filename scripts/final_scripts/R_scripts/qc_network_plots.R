#!/usr/bin/env Rscript

# Usage: Rscript qc_network_plots.R input_igraph.rds output_pdf

suppressMessages({
  library(igraph)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript qc_network_plots.R input_igraph.rds output_pdf")
}

input_rds <- args[1]
output_pdf <- args[2]

# Load igraph object
g <- readRDS(input_rds)
if (!inherits(g, "igraph")) stop("Input is not an igraph object!")

sample_name <- gsub("\\.rds$", "", basename(input_rds))

# DEGREE & STRENGTH
deg <- degree(g)

# Be robust if edge weights are missing
w_all <- E(g)$weight
if (is.null(w_all)) w_all <- rep(1, ecount(g))
str <- strength(g, weights = w_all)

# COMPONENTS
comp <- components(g)
num_comp <- comp$no
comp_sizes <- comp$csize
top10_idx <- order(comp_sizes, decreasing = TRUE)[1:min(10, length(comp_sizes))]
top10_sizes <- comp_sizes[top10_idx]
top10_table <- data.frame(Rank = 1:length(top10_sizes), Size = top10_sizes)

# Prepare summary statistics for moments
deg_stats <- sprintf(
  "Mean=%.2f | Median=%d | SD=%.2f | Max=%d | Min=%d",
  mean(deg), median(deg), sd(deg), max(deg), min(deg)
)
str_stats <- sprintf(
  "Mean=%.2f | Median=%.2f | SD=%.2f | Max=%.2f | Min=%.2f",
  mean(str), median(str), sd(str), max(str), min(str)
)

# --- Histogram binning ---

# Degree: fixed max number of bins for clarity
max_deg_bins <- 200  # adjust if you want fewer/more
if (length(deg) > 0) {
  dmin <- min(deg, na.rm = TRUE)
  dmax <- max(deg, na.rm = TRUE)
  if (is.finite(dmin) && is.finite(dmax) && dmax > dmin) {
    deg_breaks <- seq(dmin, dmax, length.out = max_deg_bins + 1)
  } else {
    # fallback when all degrees identical or not finite
    deg_breaks <- max_deg_bins
  }
} else {
  deg_breaks <- max_deg_bins
}

# Strength: many bins using a capped Freedman–Diaconis rule
fd_bins <- function(x, min_bins = 80, max_bins = 300) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(min_bins)
  iqr <- IQR(x)
  bw <- 2 * iqr / (n^(1/3))
  rng <- range(x)
  if (!is.finite(bw) || bw <= 0 || !is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] == rng[1]) {
    return(max(min_bins, min(max_bins, length(unique(round(x, 6))))))
  }
  k <- ceiling((rng[2] - rng[1]) / bw)
  k <- max(min_bins, min(max_bins, k))
  return(k)
}
str_breaks <- fd_bins(str)

# Open PDF device for multi-page report
pdf(output_pdf, width = 9, height = 6)

# DEGREE DISTRIBUTION
hist(deg, breaks = deg_breaks, col = "#396AB1", border = NA,
     main = sprintf("Degree Distribution\n%s", sample_name),
     xlab = "Degree", ylab = "Frequency", sub = deg_stats)
mtext(sprintf("Num. components: %d | Largest: %d", num_comp, top10_sizes[1]),
      side = 3, line = 0.3, cex = 0.8)

# STRENGTH DISTRIBUTION
hist(str, breaks = str_breaks, col = "#DA7C30", border = NA,
     main = sprintf("Strength Distribution\n%s", sample_name),
     xlab = "Strength", ylab = "Frequency", sub = str_stats)
mtext(sprintf("Num. components: %d | Largest: %d", num_comp, top10_sizes[1]),
      side = 3, line = 0.3, cex = 0.8)

# TABLE OF TOP 10 COMPONENT SIZES
plot.new()
title(main = sprintf("Top 10 Largest Components\n%s", sample_name))
table_str <- capture.output(print(top10_table, row.names = FALSE))
text(0, 1, paste(table_str, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 1.2)

# PLOT FIRST 4 LARGEST COMPONENTS (subgraphs)
largest <- order(comp$csize, decreasing = TRUE)[1:min(4, comp$no)]
par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
for (i in largest) {
  subg <- induced_subgraph(g, which(comp$membership == i))
  # Fruchterman-Reingold layout with positive weights
  w <- abs(E(subg)$weight)
  plot(
    subg,
    main = sprintf("Component %d (n=%d)", i, vcount(subg)),
    vertex.label = NA, vertex.size = 5, edge.arrow.size = 0.3,
    layout = layout_with_fr(subg, weights = w)
  )
}
par(mfrow = c(1, 1)) # Reset

dev.off()

cat(sprintf("Plots for %s saved to %s\n", sample_name, output_pdf))
cat(sprintf("Number of components: %d\n", num_comp))
cat("Top 10 largest components:\n")
print(top10_table)
