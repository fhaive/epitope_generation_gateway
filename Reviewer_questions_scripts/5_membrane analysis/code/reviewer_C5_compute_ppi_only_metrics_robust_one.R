#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

sample <- Sys.getenv("SAMPLE", unset = "")
if (!nzchar(sample)) {
  stop("Set SAMPLE environment variable")
}

IN_DIR <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks/filtered_networks_rds"
)

OUT_BASE <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks"
)

DEG_DIR <- file.path(OUT_BASE, "Network_Metrics_Degree")
STR_DIR <- file.path(OUT_BASE, "Network_Metrics_Strength")
BET_DIR <- file.path(OUT_BASE, "Network_Metrics_Betweenness")
LCI_DIR <- file.path(OUT_BASE, "Network_Metrics_LargestComponentImpact")
QC_DIR  <- file.path(OUT_BASE, "metric_qc")
LOG_DIR <- file.path(OUT_BASE, "logs_metrics_robust_parallel")

dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STR_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LCI_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

sanitize_graph_weights <- function(g) {
  w <- E(g)$weight

  if (is.null(w)) {
    w <- rep(1, ecount(g))
  }

  n_bad <- sum(!is.finite(w))
  n_zero <- sum(is.finite(w) & abs(w) == 0)

  w[!is.finite(w)] <- 0
  E(g)$weight <- as.numeric(w)

  list(
    graph = g,
    n_bad_weights = n_bad,
    n_zero_weights = n_zero
  )
}

largest_component_impact <- function(graph) {
  n <- vcount(graph)
  impact_vals <- numeric(n)
  names(impact_vals) <- V(graph)$name

  if (n == 0) {
    return(impact_vals)
  }

  if (ecount(graph) == 0) {
    impact_vals[] <- 1.0
    return(impact_vals)
  }

  w <- E(graph)$weight
  w[!is.finite(w)] <- 0
  E(graph)$weight <- abs(w)

  comp <- components(graph)
  largest_cc_idx <- which(comp$membership == which.max(comp$csize))
  baseline_cc_size <- length(largest_cc_idx)

  baseline_cc <- induced_subgraph(graph, largest_cc_idx)
  baseline_weight <- sum(E(baseline_cc)$weight, na.rm = TRUE)

  for (v in seq_len(n)) {
    g_temp <- delete_vertices(graph, v)

    if (vcount(g_temp) > 0 && ecount(g_temp) > 0) {
      comp_temp <- components(g_temp)
      new_cc_size <- max(comp_temp$csize)
      new_cc_idx <- which(comp_temp$membership == which.max(comp_temp$csize))
      new_cc <- induced_subgraph(g_temp, new_cc_idx)
      new_weight <- sum(E(new_cc)$weight, na.rm = TRUE)

      size_impact <- (baseline_cc_size - new_cc_size) / baseline_cc_size

      weight_impact <- if (baseline_weight > 0) {
        (baseline_weight - new_weight) / baseline_weight
      } else {
        0
      }

      impact_vals[v] <- 0.5 * size_impact + 0.5 * weight_impact
    } else {
      impact_vals[v] <- 1.0
    }
  }

  impact_vals
}

in_rds <- file.path(IN_DIR, paste0(sample, "_ppi_only_filtered.rds"))

degree_file <- file.path(DEG_DIR, paste0(sample, "_degree.tsv"))
strength_file <- file.path(STR_DIR, paste0(sample, "_strength.tsv"))
betweenness_file <- file.path(BET_DIR, paste0(sample, "_betweenness.tsv"))
impact_file <- file.path(LCI_DIR, paste0(sample, "_impact.tsv"))
status_file <- file.path(QC_DIR, paste0(sample, "_ppi_only_metric_status.tsv"))

if (!file.exists(in_rds)) {
  stop("Missing input graph: ", in_rds)
}

if (
  file.exists(degree_file) &&
  file.exists(strength_file) &&
  file.exists(betweenness_file) &&
  file.exists(impact_file)
) {
  message2("Skipping existing metrics for ", sample)
  quit(save = "no")
}

task_start <- Sys.time()

message2("Reading graph: ", in_rds)
g <- readRDS(in_rds)

clean <- sanitize_graph_weights(g)
g <- clean$graph

genes <- V(g)$name

w <- E(g)$weight
w[!is.finite(w)] <- 0

message2("Computing degree")
degree_dt <- data.table(
  gene = genes,
  degree = as.numeric(degree(g, mode = "all"))
)
fwrite(degree_dt, degree_file, sep = "\t")

message2("Computing strength")
strength_dt <- data.table(
  gene = genes,
  strength = as.numeric(strength(g, mode = "all", weights = abs(w)))
)
fwrite(strength_dt, strength_file, sep = "\t")

message2("Computing betweenness")
eps <- .Machine$double.eps
w_len <- 1 / pmax(abs(w), eps)
w_len[!is.finite(w_len)] <- 1 / eps

betweenness_dt <- data.table(
  gene = genes,
  betweenness = as.numeric(betweenness(
    g,
    directed = FALSE,
    weights = w_len,
    normalized = TRUE
  ))
)
fwrite(betweenness_dt, betweenness_file, sep = "\t")

message2("Computing LCI")
impact <- largest_component_impact(g)

impact_dt <- data.table(
  gene = names(impact),
  largest_component_impact = as.numeric(impact)
)
fwrite(impact_dt, impact_file, sep = "\t")

status <- data.table(
  sample = sample,
  status = "success",
  vertices = vcount(g),
  edges = ecount(g),
  nonfinite_weights_replaced_with_zero = clean$n_bad_weights,
  finite_zero_weights = clean$n_zero_weights,
  elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs"))
)

fwrite(status, status_file, sep = "\t")

message2("Done ", sample)
print(status)
