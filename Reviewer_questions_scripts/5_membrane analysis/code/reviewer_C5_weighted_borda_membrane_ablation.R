#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

FINAL_DIR <- file.path(ROOT, "TCGA_melanoma/epitopes_prioritisation_complete_rescue_with_vaf/final_epitopes")

PPI_ONLY_DIR <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks"
)

MEMBRANE_RDS <- file.path(ROOT, "TCGA_melanoma/sample_specific_networks/gene_lists/membrane_ensembl.rds")

OUT <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/weighted_borda_ppi_only_ablation"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "ablation_final_epitopes"), recursive = TRUE, showWarnings = FALSE)

TOP_N <- c(10, 20, 50, 100)

FEATURES <- data.table(
  feature = c(
    "Median.MT.IC50.Score",
    "Depmap_survivability_score",
    "Net_Betweenness",
    "Net_Degree",
    "Net_Impact",
    "Net_Strength",
    "Net_WCI"
  ),
  weight = c(0.25, 0.25, 0.10, 0.10, 0.10, 0.10, 0.10),
  direction = c("lower", "lower", "higher", "higher", "higher", "higher", "higher")
)

PPI_REPLACEMENTS <- c(
  Net_Betweenness = "ppi_only_betweenness",
  Net_Degree = "ppi_only_degree",
  Net_Impact = "ppi_only_impact",
  Net_Strength = "ppi_only_strength"
)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

std_ens <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

std_sym <- function(x) {
  trimws(toupper(as.character(x)))
}

sample_from_file <- function(f) {
  x <- basename(f)
  x <- sub("\\.csv$", "", x)
  x <- sub("\\.tsv$", "", x)
  x <- sub("_epitopes_final$", "", x)
  x
}

make_epitope_id <- function(dt) {
  paste(
    dt[["Gene.Name"]],
    dt[["HLA.Allele"]],
    dt[["Peptide.Length"]],
    dt[["MT.Epitope.Seq"]],
    dt[["WT.Epitope.Seq"]],
    sep = "|"
  )
}

read_membrane_symbols <- function() {
  x <- unique(as.character(readRDS(MEMBRANE_RDS)))

  if (mean(grepl("^ENSG", x)) > 0.5) {
    ens <- unique(std_ens(x))

    map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ens,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "SYMBOL")
    )

    map <- as.data.table(map)
    map <- map[!is.na(SYMBOL) & nzchar(SYMBOL)]
    unique(std_sym(map$SYMBOL))
  } else {
    unique(std_sym(x))
  }
}

read_metric_file <- function(sample_id, metric_name) {
  if (metric_name == "degree") {
    f <- file.path(PPI_ONLY_DIR, "Network_Metrics_Degree", paste0(sample_id, "_degree.tsv"))
    value_col <- "degree"
  } else if (metric_name == "strength") {
    f <- file.path(PPI_ONLY_DIR, "Network_Metrics_Strength", paste0(sample_id, "_strength.tsv"))
    value_col <- "strength"
  } else if (metric_name == "betweenness") {
    f <- file.path(PPI_ONLY_DIR, "Network_Metrics_Betweenness", paste0(sample_id, "_betweenness.tsv"))
    value_col <- "betweenness"
  } else if (metric_name == "impact") {
    f <- file.path(PPI_ONLY_DIR, "Network_Metrics_LargestComponentImpact", paste0(sample_id, "_impact.tsv"))
    value_col <- "largest_component_impact"
  } else {
    stop("Unknown metric: ", metric_name)
  }

  if (!file.exists(f)) {
    stop("Missing metric file: ", f)
  }

  dt <- fread(f)

  if (!all(c("gene", value_col) %in% names(dt))) {
    stop("Missing columns in metric file: ", f)
  }

  dt[, ENSEMBL := std_ens(gene)]

  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(dt$ENSEMBL),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )

  map <- as.data.table(map)
  map <- map[!is.na(SYMBOL) & nzchar(SYMBOL)]
  map[, gene_symbol := std_sym(SYMBOL)]

  out <- merge(
    dt[, .(ENSEMBL, value = as.numeric(get(value_col)))],
    map[, .(ENSEMBL, gene_symbol)],
    by = "ENSEMBL",
    allow.cartesian = TRUE
  )

  out[, .(value = max(value, na.rm = TRUE)), by = gene_symbol]
}

read_ppi_metrics <- function(sample_id) {
  x1 <- read_metric_file(sample_id, "degree")
  setnames(x1, "value", "ppi_only_degree")

  x2 <- read_metric_file(sample_id, "strength")
  setnames(x2, "value", "ppi_only_strength")

  x3 <- read_metric_file(sample_id, "betweenness")
  setnames(x3, "value", "ppi_only_betweenness")

  x4 <- read_metric_file(sample_id, "impact")
  setnames(x4, "value", "ppi_only_impact")

  Reduce(function(a, b) merge(a, b, by = "gene_symbol", all = TRUE), list(x1, x2, x3, x4))
}

borda_points <- function(x, direction) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  if (direction == "higher") {
    rank(x, ties.method = "average", na.last = "keep")
  } else if (direction == "lower") {
    rank(-x, ties.method = "average", na.last = "keep")
  } else {
    stop("Bad direction: ", direction)
  }
}

compute_borda <- function(dt, feature_map = NULL) {
  numerator <- rep(0, nrow(dt))
  denominator <- rep(0, nrow(dt))

  for (i in seq_len(nrow(FEATURES))) {
    original_feature <- FEATURES$feature[i]
    use_feature <- original_feature

    if (!is.null(feature_map) && original_feature %in% names(feature_map)) {
      use_feature <- unname(feature_map[[original_feature]])
    }

    if (!(use_feature %in% names(dt))) {
      stop("Missing feature column: ", use_feature)
    }

    pts <- borda_points(dt[[use_feature]], FEATURES$direction[i])
    ok <- !is.na(pts)

    numerator[ok] <- numerator[ok] + FEATURES$weight[i] * pts[ok]
    denominator[ok] <- denominator[ok] + FEATURES$weight[i]
  }

  ifelse(denominator > 0, numerator / denominator, NA_real_)
}

top_overlap <- function(dt, k, rank_a, rank_b, sample_id, label) {
  a_top <- dt[order(get(rank_a))][seq_len(min(k, .N)), epitope_id]
  b_top <- dt[order(get(rank_b))][seq_len(min(k, .N)), epitope_id]

  inter <- length(intersect(a_top, b_top))
  union <- length(union(a_top, b_top))
  min_top <- min(length(a_top), length(b_top))

  a_dt <- dt[epitope_id %in% a_top]
  b_dt <- dt[epitope_id %in% b_top]

  data.table(
    sample = sample_id,
    comparison = label,
    top_n = k,
    n_epitopes = nrow(dt),
    top_overlap_n = inter,
    top_overlap_fraction = ifelse(min_top > 0, inter / min_top, NA_real_),
    top_jaccard = ifelse(union > 0, inter / union, NA_real_),
    a_top_membrane_n = sum(a_dt$is_membrane_gene, na.rm = TRUE),
    b_top_membrane_n = sum(b_dt$is_membrane_gene, na.rm = TRUE)
  )
}

files <- sort(list.files(FINAL_DIR, pattern = "_epitopes_final\\.csv$", full.names = TRUE))

if (length(files) == 0) {
  stop("No final epitope CSV files found")
}

membrane_symbols <- read_membrane_symbols()

msg("Final files: ", length(files))
msg("Membrane symbols: ", length(membrane_symbols))

status_list <- list()
validation_list <- list()
ablation_list <- list()
top_list <- list()
shift_list <- list()

for (f in files) {
  sample_id <- sample_from_file(f)
  msg("Processing ", sample_id)

  start_time <- Sys.time()

  res <- tryCatch({
    dt <- fread(f)

    required <- unique(c(
      "Gene.Name",
      "HLA.Allele",
      "Peptide.Length",
      "MT.Epitope.Seq",
      "WT.Epitope.Seq",
      "Borda_Score",
      "Borda_Rank",
      FEATURES$feature
    ))

    missing <- setdiff(required, names(dt))
    if (length(missing) > 0) {
      stop("Missing columns: ", paste(missing, collapse = ", "))
    }

    dt[, gene_symbol := std_sym(get("Gene.Name"))]
    dt[, epitope_id := make_epitope_id(.SD)]
    dt[, final_borda_score := as.numeric(get("Borda_Score"))]
    dt[, final_borda_rank := as.numeric(get("Borda_Rank"))]
    dt[, is_membrane_gene := gene_symbol %in% membrane_symbols]

    ppi <- read_ppi_metrics(sample_id)

    dt <- merge(dt, ppi, by = "gene_symbol", all.x = TRUE, sort = FALSE)

    dt[is.na(ppi_only_degree), ppi_only_degree := as.numeric(Net_Degree)]
    dt[is.na(ppi_only_strength), ppi_only_strength := as.numeric(Net_Strength)]
    dt[is.na(ppi_only_betweenness), ppi_only_betweenness := as.numeric(Net_Betweenness)]
    dt[is.na(ppi_only_impact), ppi_only_impact := as.numeric(Net_Impact)]

    dt[, reconstructed_original_borda_score := compute_borda(.SD, feature_map = NULL)]
    dt[, reconstructed_original_borda_rank := rank(-reconstructed_original_borda_score, ties.method = "first", na.last = "keep")]

    fmap <- as.list(PPI_REPLACEMENTS)

    dt[, ppi_only_weighted_borda_score := compute_borda(.SD, feature_map = fmap)]
    dt[, ppi_only_weighted_borda_rank := rank(-ppi_only_weighted_borda_score, ties.method = "first", na.last = "keep")]

    dt[, rank_gain_with_membrane := ppi_only_weighted_borda_rank - final_borda_rank]

    validation <- data.table(
      sample = sample_id,
      n_epitopes = nrow(dt),
      comparison = "final_table_vs_recomputed_original_borda",
      spearman = suppressWarnings(cor(dt$final_borda_rank, dt$reconstructed_original_borda_rank, method = "spearman", use = "complete.obs")),
      kendall = suppressWarnings(cor(dt$final_borda_rank, dt$reconstructed_original_borda_rank, method = "kendall", use = "complete.obs"))
    )

    ablation <- data.table(
      sample = sample_id,
      n_epitopes = nrow(dt),
      comparison = "final_membrane_vs_ppi_only_weighted_borda",
      spearman = suppressWarnings(cor(dt$final_borda_rank, dt$ppi_only_weighted_borda_rank, method = "spearman", use = "complete.obs")),
      kendall = suppressWarnings(cor(dt$final_borda_rank, dt$ppi_only_weighted_borda_rank, method = "kendall", use = "complete.obs"))
    )

    top_dt <- rbindlist(lapply(TOP_N, function(k) {
      rbind(
        top_overlap(dt, k, "final_borda_rank", "reconstructed_original_borda_rank", sample_id, "validation_final_vs_recomputed_original"),
        top_overlap(dt, k, "final_borda_rank", "ppi_only_weighted_borda_rank", sample_id, "ablation_final_membrane_vs_ppi_only")
      )
    }), use.names = TRUE, fill = TRUE)

    shift <- data.table(
      sample = sample_id,
      n_epitopes = nrow(dt),
      median_rank_gain_all = median(dt$rank_gain_with_membrane, na.rm = TRUE),
      median_rank_gain_membrane = median(dt[is_membrane_gene == TRUE, rank_gain_with_membrane], na.rm = TRUE),
      median_rank_gain_nonmembrane = median(dt[is_membrane_gene == FALSE, rank_gain_with_membrane], na.rm = TRUE),
      mean_rank_gain_membrane = mean(dt[is_membrane_gene == TRUE, rank_gain_with_membrane], na.rm = TRUE),
      mean_rank_gain_nonmembrane = mean(dt[is_membrane_gene == FALSE, rank_gain_with_membrane], na.rm = TRUE)
    )

    out_cols <- c(
      "Gene.Name",
      "HLA.Allele",
      "Peptide.Length",
      "MT.Epitope.Seq",
      "WT.Epitope.Seq",
      "Median.MT.IC50.Score",
      "Tumor DNA VAF",
      "Depmap_survivability_score",
      "Net_Betweenness",
      "Net_Degree",
      "Net_Impact",
      "Net_Strength",
      "Net_WCI",
      "ppi_only_betweenness",
      "ppi_only_degree",
      "ppi_only_impact",
      "ppi_only_strength",
      "Borda_Score",
      "Borda_Rank",
      "reconstructed_original_borda_score",
      "reconstructed_original_borda_rank",
      "ppi_only_weighted_borda_score",
      "ppi_only_weighted_borda_rank",
      "rank_gain_with_membrane",
      "is_membrane_gene"
    )

    out_cols <- out_cols[out_cols %in% names(dt)]

    fwrite(
      dt[order(ppi_only_weighted_borda_rank), ..out_cols],
      file.path(OUT, "ablation_final_epitopes", paste0(sample_id, "_weighted_ppi_only_ablation_epitopes.csv"))
    )

    list(
      status = data.table(
        sample = sample_id,
        status = "success",
        elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        message = NA_character_
      ),
      validation = validation,
      ablation = ablation,
      top = top_dt,
      shift = shift
    )
  }, error = function(e) {
    list(
      status = data.table(
        sample = sample_id,
        status = "failed",
        elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        message = conditionMessage(e)
      ),
      validation = data.table(),
      ablation = data.table(),
      top = data.table(),
      shift = data.table()
    )
  })

  status_list[[length(status_list) + 1]] <- res$status
  if (nrow(res$validation) > 0) validation_list[[length(validation_list) + 1]] <- res$validation
  if (nrow(res$ablation) > 0) ablation_list[[length(ablation_list) + 1]] <- res$ablation
  if (nrow(res$top) > 0) top_list[[length(top_list) + 1]] <- res$top
  if (nrow(res$shift) > 0) shift_list[[length(shift_list) + 1]] <- res$shift
}

status_dt <- rbindlist(status_list, use.names = TRUE, fill = TRUE)
validation_dt <- rbindlist(validation_list, use.names = TRUE, fill = TRUE)
ablation_dt <- rbindlist(ablation_list, use.names = TRUE, fill = TRUE)
top_dt <- rbindlist(top_list, use.names = TRUE, fill = TRUE)
shift_dt <- rbindlist(shift_list, use.names = TRUE, fill = TRUE)

fwrite(status_dt, file.path(OUT, "weighted_borda_ablation_status.tsv"), sep = "\t")
fwrite(validation_dt, file.path(OUT, "weighted_borda_validation_rank_correlations.tsv"), sep = "\t")
fwrite(ablation_dt, file.path(OUT, "weighted_borda_ppi_only_ablation_rank_correlations.tsv"), sep = "\t")
fwrite(top_dt, file.path(OUT, "weighted_borda_topk_overlap.tsv"), sep = "\t")
fwrite(shift_dt, file.path(OUT, "weighted_borda_rank_shift_summary_per_patient.tsv"), sep = "\t")

cat("\nStatus summary:\n")
print(status_dt[, .N, by = status])

if (nrow(status_dt[status == "failed"]) > 0) {
  cat("\nFailed rows:\n")
  print(status_dt[status == "failed"])
}

if (nrow(ablation_dt) == 0) {
  stop("No successful ablation results")
}

validation_summary <- validation_dt[, .(
  n_patients = .N,
  median_spearman = median(spearman, na.rm = TRUE),
  min_spearman = min(spearman, na.rm = TRUE),
  q25_spearman = quantile(spearman, 0.25, na.rm = TRUE),
  q75_spearman = quantile(spearman, 0.75, na.rm = TRUE),
  median_kendall = median(kendall, na.rm = TRUE)
), by = comparison]

ablation_summary <- ablation_dt[, .(
  n_patients = .N,
  median_spearman = median(spearman, na.rm = TRUE),
  min_spearman = min(spearman, na.rm = TRUE),
  q25_spearman = quantile(spearman, 0.25, na.rm = TRUE),
  q75_spearman = quantile(spearman, 0.75, na.rm = TRUE),
  median_kendall = median(kendall, na.rm = TRUE)
), by = comparison]

top_summary <- top_dt[, .(
  n_patients = .N,
  median_top_overlap_n = median(top_overlap_n, na.rm = TRUE),
  q25_top_overlap_n = quantile(top_overlap_n, 0.25, na.rm = TRUE),
  q75_top_overlap_n = quantile(top_overlap_n, 0.75, na.rm = TRUE),
  median_top_overlap_fraction = median(top_overlap_fraction, na.rm = TRUE),
  q25_top_overlap_fraction = quantile(top_overlap_fraction, 0.25, na.rm = TRUE),
  q75_top_overlap_fraction = quantile(top_overlap_fraction, 0.75, na.rm = TRUE),
  median_top_jaccard = median(top_jaccard, na.rm = TRUE),
  median_a_top_membrane_n = median(a_top_membrane_n, na.rm = TRUE),
  median_b_top_membrane_n = median(b_top_membrane_n, na.rm = TRUE)
), by = .(comparison, top_n)][order(comparison, top_n)]

shift_summary <- shift_dt[, .(
  n_patients = .N,
  median_rank_gain_all = median(median_rank_gain_all, na.rm = TRUE),
  median_rank_gain_membrane = median(median_rank_gain_membrane, na.rm = TRUE),
  median_rank_gain_nonmembrane = median(median_rank_gain_nonmembrane, na.rm = TRUE),
  median_mean_rank_gain_membrane = median(mean_rank_gain_membrane, na.rm = TRUE),
  median_mean_rank_gain_nonmembrane = median(mean_rank_gain_nonmembrane, na.rm = TRUE)
)]

display <- merge(
  ablation_summary,
  top_summary[comparison == "ablation_final_membrane_vs_ppi_only"],
  by = "comparison",
  all.x = TRUE
)

display <- display[, .(
  comparison,
  top_n,
  n_patients = n_patients.x,
  `Median Spearman rho` = round(median_spearman, 3),
  `Spearman IQR` = paste0(round(q25_spearman, 3), "-", round(q75_spearman, 3)),
  `Minimum Spearman rho` = round(min_spearman, 3),
  `Median Kendall tau` = round(median_kendall, 3),
  `Median top overlap n` = round(median_top_overlap_n, 1),
  `Top overlap n IQR` = paste0(round(q25_top_overlap_n, 1), "-", round(q75_top_overlap_n, 1)),
  `Median top overlap %` = round(100 * median_top_overlap_fraction, 1),
  `Top overlap % IQR` = paste0(round(100 * q25_top_overlap_fraction, 1), "-", round(100 * q75_top_overlap_fraction, 1)),
  `Median top Jaccard` = round(median_top_jaccard, 3),
  `Original top membrane n median` = round(median_a_top_membrane_n, 1),
  `PPI-only top membrane n median` = round(median_b_top_membrane_n, 1)
)][order(top_n)]

fwrite(validation_summary, file.path(OUT, "weighted_borda_validation_summary.tsv"), sep = "\t")
fwrite(ablation_summary, file.path(OUT, "weighted_borda_ppi_only_ablation_rank_summary.tsv"), sep = "\t")
fwrite(top_summary, file.path(OUT, "weighted_borda_topk_summary.tsv"), sep = "\t")
fwrite(shift_summary, file.path(OUT, "weighted_borda_rank_shift_summary_overall.tsv"), sep = "\t")
fwrite(display, file.path(OUT, "weighted_borda_ppi_only_ablation_summary_display.tsv"), sep = "\t")

cat("\nValidation summary:\n")
print(validation_summary)

cat("\nAblation summary display:\n")
print(display)

cat("\nRank shift summary:\n")
print(shift_summary)

msg("Done. Output: ", OUT)
