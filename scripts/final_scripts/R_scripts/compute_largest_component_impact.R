#!/usr/bin/env Rscript

# Usage: Rscript compute_largest_component_impact.R input_graph.rds output_impact.tsv

suppressMessages({
  library(igraph)
  library(tibble)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_tsv <- args[2]

cat("Input graph (RDS):", input_rds, "\n")
cat("Output impact TSV:", output_tsv, "\n")

# Load igraph object
g <- readRDS(input_rds)
stopifnot(inherits(g, "igraph"))

# Custom largest component impact function
largest_component_impact <- function(graph) {
  n <- vcount(graph)
  impact_vals <- numeric(n)
  names(impact_vals) <- V(graph)$name

  # Use absolute weights for percolation/impact
  weights <- abs(E(graph)$weight)
  E(graph)$weight <- weights

  # Baseline largest component
  comp <- components(graph)
  largest_cc_idx <- which(comp$membership == which.max(comp$csize))
  baseline_cc_size <- length(largest_cc_idx)
  baseline_cc <- induced_subgraph(graph, largest_cc_idx)
  baseline_weight <- sum(E(baseline_cc)$weight)

  for (v in seq_len(n)) {
    # Remove node v
    g_temp <- delete_vertices(graph, v)
    if (vcount(g_temp) > 0 && ecount(g_temp) > 0) {
      comp_temp <- components(g_temp)
      new_cc_size <- max(comp_temp$csize)
      new_cc_idx <- which(comp_temp$membership == which.max(comp_temp$csize))
      new_cc <- induced_subgraph(g_temp, new_cc_idx)
      new_weight <- sum(E(new_cc)$weight)
      size_impact <- (baseline_cc_size - new_cc_size) / baseline_cc_size
      weight_impact <- if (baseline_weight > 0) (baseline_weight - new_weight) / baseline_weight else 0
      impact_vals[v] <- 0.5 * size_impact + 0.5 * weight_impact
    } else {
      # If removal causes no edges, assign max impact
      impact_vals[v] <- 1.0
    }
    if (v %% 1000 == 0) cat("Processed", v, "nodes\n")
  }
  return(impact_vals)
}

cat("Computing largest component impact centrality...\n")
impact_vec <- largest_component_impact(g)

out_tbl <- tibble(
  gene = names(impact_vec),
  largest_component_impact = impact_vec
)
write_tsv(out_tbl, output_tsv)
cat("✓ Done. Output written to", output_tsv, "\n")
