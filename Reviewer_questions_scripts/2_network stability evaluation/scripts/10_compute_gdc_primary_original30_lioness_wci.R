#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"
WORK <- file.path(ROOT, "reviewer_C_external_benchmark")

VST_FILE <- file.path(WORK, "gdc_primary_original30_expression_matrix", "gdc_primary_original30_vst_mat.tsv")

OUT <- file.path(WORK, "gdc_primary_original30_lioness_wci")
NET_DIR <- file.path(OUT, "Sample_Specific_Networks")
WCI_DIR <- file.path(OUT, "Network_Metrics_Full", "WCI")
STR_DIR <- file.path(OUT, "Network_Metrics_Full", "Strength")

dir.create(NET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WCI_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STR_DIR, recursive = TRUE, showWarnings = FALSE)

TARGETS <- c(
  "Sample_TCGA-BF-A1PV",
  "Sample_TCGA-BF-A3DJ",
  "Sample_TCGA-BF-A3DL",
  "Sample_TCGA-BF-A3DM",
  "Sample_TCGA-BF-A3DN",
  "Sample_TCGA-BF-A5EO",
  "Sample_TCGA-BF-A5EP",
  "Sample_TCGA-BF-A5EQ",
  "Sample_TCGA-BF-A5ER",
  "Sample_TCGA-BF-A5ES",
  "Sample_TCGA-BF-AAOU",
  "Sample_TCGA-BF-AAOX",
  "Sample_TCGA-BF-AAP1",
  "Sample_TCGA-BF-AAP2",
  "Sample_TCGA-BF-AAP4",
  "Sample_TCGA-BF-AAP6",
  "Sample_TCGA-BF-AAP7",
  "Sample_TCGA-BF-AAP8",
  "Sample_TCGA-D3-A5GT",
  "Sample_TCGA-D9-A3Z4",
  "Sample_TCGA-D9-A4Z2",
  "Sample_TCGA-D9-A4Z3",
  "Sample_TCGA-EB-A1NK",
  "Sample_TCGA-EB-A3HV",
  "Sample_TCGA-EB-A3XB",
  "Sample_TCGA-EB-A3XC",
  "Sample_TCGA-EB-A3XD",
  "Sample_TCGA-EB-A3XE",
  "Sample_TCGA-EB-A3XF",
  "Sample_TCGA-EB-A41A"
)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

read_vst <- function(f) {
  dt <- fread(f)
  gene_col <- names(dt)[1]
  genes <- dt[[gene_col]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

cor_from_stats <- function(cross, sums, sumsqs, n) {
  centered_cross <- cross - tcrossprod(sums) / n

  ss <- sumsqs - (sums^2) / n
  ss[ss <= 0] <- NA_real_

  denom <- sqrt(ss)

  centered_cross <- sweep(centered_cross, 1, denom, "/")
  centered_cross <- sweep(centered_cross, 2, denom, "/")

  centered_cross[!is.finite(centered_cross)] <- 0
  diag(centered_cross) <- 1

  centered_cross
}

compute_wci_and_strength <- function(S, sample_id) {
  genes <- rownames(S)

  degree <- rowSums(S != 0)
  strength <- rowSums(abs(S))
  total_weight <- sum(abs(S))

  wci <- if (total_weight > 0) {
    strength / total_weight
  } else {
    rep(NA_real_, length(strength))
  }

  fwrite(
    data.table(
      gene = genes,
      degree = as.numeric(degree),
      strength = as.numeric(strength)
    ),
    file.path(STR_DIR, paste0(sample_id, "_strength.tsv")),
    sep = "\t"
  )

  fwrite(
    data.table(
      gene = genes,
      wci = as.numeric(wci)
    ),
    file.path(WCI_DIR, paste0(sample_id, "_wci.tsv")),
    sep = "\t"
  )
}

message2("Reading VST: ", VST_FILE)
mat <- read_vst(VST_FILE)

missing_targets <- setdiff(TARGETS, colnames(mat))
if (length(missing_targets) > 0) {
  stop("Missing target columns: ", paste(missing_targets, collapse = ", "))
}

mat <- mat[, TARGETS, drop = FALSE]
mat <- mat[order(rownames(mat)), , drop = FALSE]

N <- ncol(mat)
G <- nrow(mat)

message2("Using GDC-primary original-30 matrix: ", G, " genes x ", N, " samples")

message2("Computing full n=30 crossproduct")
cross_all <- tcrossprod(mat)

sums_all <- rowSums(mat)
sumsqs_all <- rowSums(mat^2)

message2("Computing full n=30 Pearson network")
G_alpha <- cor_from_stats(cross_all, sums_all, sumsqs_all, N)
rownames(G_alpha) <- rownames(mat)
colnames(G_alpha) <- rownames(mat)

manifest <- list()

for (i in seq_along(TARGETS)) {
  sample_id <- TARGETS[i]

  out_net <- file.path(NET_DIR, paste0(sample_id, ".rds"))
  out_wci <- file.path(WCI_DIR, paste0(sample_id, "_wci.tsv"))

  if (file.exists(out_net) && file.exists(out_wci)) {
    message2("Skipping existing ", i, "/", length(TARGETS), ": ", sample_id)
    manifest[[length(manifest) + 1]] <- data.table(
      target_index = i,
      sample = sample_id,
      status = "skipped_existing",
      network_file = out_net,
      wci_file = out_wci
    )
    next
  }

  message2("Computing target ", i, "/", length(TARGETS), ": ", sample_id)

  task_start <- Sys.time()

  xq <- mat[, sample_id]

  sums_minus <- sums_all - xq
  sumsqs_minus <- sumsqs_all - xq^2
  cross_minus <- cross_all - tcrossprod(xq)

  G_minus <- cor_from_stats(cross_minus, sums_minus, sumsqs_minus, N - 1)
  rownames(G_minus) <- rownames(mat)
  colnames(G_minus) <- rownames(mat)

  E_q <- N * (G_alpha - G_minus) + G_minus
  rownames(E_q) <- rownames(mat)
  colnames(E_q) <- rownames(mat)

  saveRDS(E_q, out_net, compress = FALSE)
  compute_wci_and_strength(E_q, sample_id)

  elapsed <- as.numeric(difftime(Sys.time(), task_start, units = "secs"))

  manifest[[length(manifest) + 1]] <- data.table(
    target_index = i,
    sample = sample_id,
    status = "success",
    elapsed_seconds = elapsed,
    network_file = out_net,
    wci_file = out_wci
  )

  rm(xq, sums_minus, sumsqs_minus, cross_minus, G_minus, E_q)
  gc(verbose = FALSE)
}

manifest_dt <- rbindlist(manifest, use.names = TRUE, fill = TRUE)

fwrite(manifest_dt, file.path(OUT, "gdc_primary_original30_lioness_wci_manifest.tsv"), sep = "\t")

qc <- data.table(
  benchmark_type = "GDC-primary original 30 only",
  samples = N,
  genes = G,
  target_samples = length(TARGETS),
  output_folder = OUT
)

fwrite(qc, file.path(OUT, "gdc_primary_original30_lioness_wci_qc.tsv"), sep = "\t")

message2("Done. Output: ", OUT)
print(qc)
print(manifest_dt[, .N, by = status])
