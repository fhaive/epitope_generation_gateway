#!/usr/bin/env Rscript

suppressMessages({
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop("Usage: Enrich_Membrane_Edges.R <raw_network_rds> <ppi_filtered_rds> <output_rds> <output_mm> <output_genes> <qc_tsv> <config_file> <humannet_file> <mem_genes_file>")
}

raw_network_rds  <- args[1]
ppi_filtered_rds <- args[2]
output_rds       <- args[3]
output_mm        <- args[4]
output_genes     <- args[5]
qc_tsv           <- args[6]
config_file      <- args[7]
humannet_file    <- args[8]
mem_genes_file   <- args[9]

cat("Raw network RDS:", raw_network_rds, "\n")
cat("PPI-filtered RDS:", ppi_filtered_rds, "\n")
cat("Output final filtered RDS:", output_rds, "\n")
cat("Matrix Market:", output_mm, "\n")
cat("Gene list:", output_genes, "\n")
cat("QC TSV:", qc_tsv, "\n")
cat("Config File:", config_file, "\n")
cat("HumanNet File:", humannet_file, "\n")
cat("Membrane Genes File:", mem_genes_file, "\n")

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_mm), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_genes), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(qc_tsv), recursive = TRUE, showWarnings = FALSE)

test_network <- readRDS(raw_network_rds)
g_ppi_top <- readRDS(ppi_filtered_rds)
HumanNet <- readRDS(humannet_file)
mem_genes <- readRDS(mem_genes_file)
config <- yaml::read_yaml(config_file)

enabled <- TRUE
add_prop <- 0.05
denominator <- "unfiltered_ppi_edges"
require_membrane_endpoint <- TRUE
edge_selection_metric <- "absolute_weight"

if (!is.null(config$membrane_enrichment)) {
  if (!is.null(config$membrane_enrichment$enabled)) {
    enabled <- isTRUE(config$membrane_enrichment$enabled)
  }
  if (!is.null(config$membrane_enrichment$add_prop)) {
    add_prop <- as.numeric(config$membrane_enrichment$add_prop)
  }
  if (!is.null(config$membrane_enrichment$denominator)) {
    denominator <- as.character(config$membrane_enrichment$denominator)
  }
  if (!is.null(config$membrane_enrichment$require_membrane_endpoint)) {
    require_membrane_endpoint <- isTRUE(config$membrane_enrichment$require_membrane_endpoint)
  }
  if (!is.null(config$membrane_enrichment$edge_selection_metric)) {
    edge_selection_metric <- as.character(config$membrane_enrichment$edge_selection_metric)
  }
} else if (!is.null(config$mem_prop)) {
  # Backward-compatible interpretation of old reverse-style mem_prop.
  # Example: old mem_prop = 0.95 means add_prop = 0.05.
  add_prop <- 1 - as.numeric(config$mem_prop)
}

if (is.na(add_prop) || add_prop < 0 || add_prop > 1) {
  stop("membrane_enrichment.add_prop must be >= 0 and <= 1")
}

cat(sprintf("Membrane enrichment enabled: %s\n", enabled))
cat(sprintf("Membrane add_prop: %.4f\n", add_prop))
cat(sprintf("Denominator: %s\n", denominator))
cat(sprintf("Require membrane endpoint: %s\n", require_membrane_endpoint))
cat(sprintf("Edge selection metric: %s\n", edge_selection_metric))

if (edge_selection_metric != "absolute_weight") {
  stop("Only membrane_enrichment.edge_selection_metric = 'absolute_weight' is currently supported")
}

# -------- 1. Rebuild weighted unfiltered PPI graph --------
# This preserves the old denominator logic:
# n_to_add = ceiling(0.05 * number_of_weighted_unfiltered_PPI_edges)

ppi_edges <- HumanNet %>%
  select(from = ensembl1, to = ensembl2) %>%
  filter(from != to)

g_ppi <- igraph::graph_from_data_frame(ppi_edges, directed = FALSE)
edge_df <- igraph::as_data_frame(g_ppi, what = "edges")

idx_from <- match(edge_df$from, rownames(test_network))
idx_to   <- match(edge_df$to, colnames(test_network))

edge_df$weight <- NA_real_

valid1 <- !is.na(idx_from) & !is.na(idx_to)
edge_df$weight[valid1] <- test_network[cbind(idx_from[valid1], idx_to[valid1])]

missing <- is.na(edge_df$weight)

if (any(missing)) {
  idx_from2 <- match(edge_df$to[missing], rownames(test_network))
  idx_to2   <- match(edge_df$from[missing], colnames(test_network))
  valid2 <- !is.na(idx_from2) & !is.na(idx_to2)
  edge_df$weight[which(missing)[valid2]] <- test_network[cbind(idx_from2[valid2], idx_to2[valid2])]
}

E(g_ppi)$weight <- edge_df$weight
g_ppi_weighted <- delete_edges(g_ppi, which(is.na(E(g_ppi)$weight)))

cat(sprintf(
  "Weighted unfiltered PPI graph: %d nodes, %d edges\n",
  vcount(g_ppi_weighted),
  ecount(g_ppi_weighted)
))

# -------- 2. Optionally add membrane-touched edges --------

ppi_top_df <- igraph::as_data_frame(g_ppi_top, what = "edges")

n_candidates <- 0
n_to_add <- 0
n_added <- 0

if (!enabled || add_prop == 0) {

  cat("Membrane enrichment disabled. Returning PPI-only graph.\n")
  final_graph <- g_ppi_top

} else {

  all_genes <- rownames(test_network)

  if (!identical(rownames(test_network), colnames(test_network))) {
    stop("Raw network matrix rownames and colnames are not identical")
  }

  mem_genes_present <- intersect(mem_genes, all_genes)

  if (length(mem_genes_present) == 0) {

    warning("No membrane genes found in this sample network. Returning PPI-only graph.")
    final_graph <- g_ppi_top

  } else {

    if (!require_membrane_endpoint) {
      stop("Only require_membrane_endpoint = true is currently supported")
    }

    # Old logic: candidate edges are all raw-network edges touching at least one membrane gene.
    mem_idx <- which(all_genes %in% mem_genes_present)

    mat_idx_mem_row <- expand.grid(i = mem_idx, j = seq_len(nrow(test_network)))
    mat_idx_mem_col <- expand.grid(i = seq_len(nrow(test_network)), j = mem_idx)

    mem_edges_idx <- rbind(mat_idx_mem_row, mat_idx_mem_col)
    mem_edges_idx <- mem_edges_idx[mem_edges_idx$i != mem_edges_idx$j, ]

    # Deduplicate undirected edges.
    a <- pmin(mem_edges_idx$i, mem_edges_idx$j)
    b <- pmax(mem_edges_idx$i, mem_edges_idx$j)
    edges_unique <- unique(cbind(i = a, j = b))

    mem_edge_df <- data.frame(
      from = all_genes[edges_unique[, "i"]],
      to = all_genes[edges_unique[, "j"]],
      weight = test_network[edges_unique],
      stringsAsFactors = FALSE
    )

    mem_edge_df <- mem_edge_df[mem_edge_df$from != mem_edge_df$to, ]

    # Remove candidate membrane edges already present in the PPI-filtered graph.
    ppi_top_edge_keys <- paste(
      pmin(ppi_top_df$from, ppi_top_df$to),
      pmax(ppi_top_df$from, ppi_top_df$to),
      sep = "_"
    )

    mem_edge_keys <- paste(
      pmin(mem_edge_df$from, mem_edge_df$to),
      pmax(mem_edge_df$from, mem_edge_df$to),
      sep = "_"
    )

    new_mem_edges <- mem_edge_df[!mem_edge_keys %in% ppi_top_edge_keys, ]
    n_candidates <- nrow(new_mem_edges)

    if (denominator == "unfiltered_ppi_edges") {
      denom_n <- ecount(g_ppi_weighted)
    } else if (denominator == "ppi_filtered_edges") {
      denom_n <- ecount(g_ppi_top)
    } else {
      stop("Unsupported membrane_enrichment.denominator. Use 'unfiltered_ppi_edges' or 'ppi_filtered_edges'.")
    }

    n_to_add <- ceiling(add_prop * denom_n)

    # Safer than the old script: avoids selecting beyond available candidate edges.
    n_to_add <- min(n_to_add, n_candidates)

    cat(sprintf(
      "Will add %d membrane-touched edges using add_prop %.4f and denominator %s = %d\n",
      n_to_add,
      add_prop,
      denominator,
      denom_n
    ))

    if (n_to_add > 0) {

      top_mem_edges <- new_mem_edges[
        order(abs(new_mem_edges$weight), decreasing = TRUE),
      ][seq_len(n_to_add), , drop = FALSE]

      all_edges <- rbind(
        ppi_top_df[, c("from", "to", "weight")],
        top_mem_edges[, c("from", "to", "weight")]
      )

      all_edge_keys <- paste(
        pmin(all_edges$from, all_edges$to),
        pmax(all_edges$from, all_edges$to),
        sep = "_"
      )

      all_edges <- all_edges[!duplicated(all_edge_keys), ]

      # Preserve old behavior:
      # include all membrane genes as vertices, even if some are isolated.
      all_vertex_names <- unique(c(
        V(g_ppi_top)$name,
        mem_genes,
        top_mem_edges$from,
        top_mem_edges$to
      ))

      final_graph <- igraph::graph_from_data_frame(
        all_edges,
        directed = FALSE,
        vertices = all_vertex_names
      )

      n_added <- nrow(top_mem_edges)

    } else {

      cat("No membrane edges added. Returning PPI-only graph.\n")
      final_graph <- g_ppi_top

    }
  }
}

cat(sprintf("Final graph: %d nodes, %d edges\n", vcount(final_graph), ecount(final_graph)))

# -------- 3. Save outputs --------

saveRDS(final_graph, file = output_rds)

adjmat <- as_adjacency_matrix(final_graph, attr = "weight", sparse = TRUE, type = "both")
writeMM(adjmat, output_mm)
writeLines(rownames(adjmat), output_genes)

qc <- data.frame(
  sample = sub("_filtered\\.rds$|\\.rds$", "", basename(output_rds)),
  membrane_enrichment_enabled = enabled,
  add_prop = add_prop,
  denominator = denominator,
  require_membrane_endpoint = require_membrane_endpoint,
  edge_selection_metric = edge_selection_metric,
  n_edges_unfiltered_ppi_weighted = ecount(g_ppi_weighted),
  n_edges_ppi_only = ecount(g_ppi_top),
  n_membrane_candidates = n_candidates,
  n_membrane_edges_requested = n_to_add,
  n_membrane_edges_added = n_added,
  n_edges_final = ecount(final_graph),
  n_vertices_ppi_only = vcount(g_ppi_top),
  n_vertices_final = vcount(final_graph)
)

write.table(qc, file = qc_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✓ Done membrane enrichment\n")
