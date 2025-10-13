#!/usr/bin/env Rscript

suppressMessages({
  library(igraph)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
rds_file    <- args[1]
out_betw    <- args[2]
out_degree  <- args[3]
out_strength<- args[4]

cat("Reading igraph object from RDS...\n")
g <- readRDS(rds_file)
genes <- V(g)$name

# Pull weights (default to 1 if missing)
w <- E(g)$weight
if (is.null(w)) w <- rep(1, ecount(g))

cat("Calculating degree (unweighted)...\n")
deg <- degree(g, mode = "all")
write_tsv(data.frame(gene = genes, degree = as.numeric(deg)), out_degree)

cat("Calculating strength with ABSOLUTE weights...\n")
# Strength = sum of incident edge weights; here we use |weight|
stren_abs <- strength(g, mode = "all", weights = abs(w))
write_tsv(data.frame(gene = genes, strength = as.numeric(stren_abs)), out_strength)

cat("Calculating betweenness centrality (can be slow)...\n")
# Betweenness in igraph interprets 'weights' as edge *lengths* (larger = longer path).
# We want strong-magnitude edges to be "short", irrespective of sign:
#   length = 1 / |weight|, guarding against zeros to avoid Inf.
w_len <- 1 / pmax(abs(w), .Machine$double.eps)
betw <- betweenness(g, directed = FALSE, weights = w_len, normalized = TRUE)
write_tsv(data.frame(gene = genes, betweenness = as.numeric(betw)), out_betw)

cat("All metrics computed and saved!\n")
