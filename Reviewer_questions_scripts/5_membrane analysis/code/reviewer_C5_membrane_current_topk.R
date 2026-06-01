#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

FINAL_DIR <- file.path(
  ROOT,
  "TCGA_melanoma/epitopes_prioritisation_complete_rescue_with_vaf/final_epitopes"
)

MEMBRANE_RDS <- file.path(
  ROOT,
  "TCGA_melanoma/sample_specific_networks/gene_lists/membrane_ensembl.rds"
)

OUT <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_N <- c(10, 20, 50, 100)
TOP_PROP <- c(0.01, 0.05, 0.10, 0.20)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

standardize_ensembl <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+$", "", x)
}

standardize_symbol <- function(x) {
  x <- as.character(x)
  trimws(toupper(x))
}

sample_from_file <- function(f) {
  x <- basename(f)
  x <- sub("\\.gz$", "", x)
  x <- sub("\\.tsv$", "", x)
  x <- sub("\\.csv$", "", x)
  x <- sub("_epitopes_final$", "", x)
  x
}

read_any <- function(f) {
  if (grepl("\\.tsv(\\.gz)?$", f)) {
    fread(f)
  } else if (grepl("\\.csv(\\.gz)?$", f)) {
    fread(f)
  } else {
    NULL
  }
}

message2("Reading membrane list: ", MEMBRANE_RDS)

membrane_raw <- readRDS(MEMBRANE_RDS)
membrane_raw <- unique(as.character(membrane_raw))

looks_ensembl <- mean(grepl("^ENSG", membrane_raw)) > 0.5

if (looks_ensembl) {
  message2("Membrane list appears to use Ensembl IDs. Mapping to HGNC symbols.")

  membrane_ens <- unique(standardize_ensembl(membrane_raw))

  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = membrane_ens,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )

  map <- as.data.table(map)
  map <- map[!is.na(SYMBOL) & nzchar(SYMBOL)]
  membrane_symbols <- unique(standardize_symbol(map$SYMBOL))

  fwrite(
    map,
    file.path(OUT, "membrane_ensembl_to_symbol_mapping.tsv"),
    sep = "\t"
  )

} else {
  message2("Membrane list does not look like Ensembl IDs. Treating as symbols.")
  membrane_symbols <- unique(standardize_symbol(membrane_raw))
}

message2("Membrane symbols available after mapping: ", length(membrane_symbols))

if (length(membrane_symbols) == 0) {
  stop("No membrane symbols available after mapping. Check membrane list.")
}

files <- list.files(FINAL_DIR, pattern = "\\.(tsv|csv)(\\.gz)?$", full.names = TRUE)
files <- sort(files)

if (length(files) == 0) {
  stop("No final epitope TSV/CSV files found in: ", FINAL_DIR)
}

message2("Final epitope files found: ", length(files))

per_patient <- list()
debug_cols <- list()
annotated_examples <- list()

for (f in files) {
  dt <- read_any(f)
  if (is.null(dt) || nrow(dt) == 0) next

  sample_id <- sample_from_file(f)

  gene_col <- "Gene.Name"
  rank_col <- "Borda_Rank"

  if (!(gene_col %in% names(dt))) {
    stop("Missing Gene.Name in: ", f)
  }

  if (!(rank_col %in% names(dt))) {
    stop("Missing Borda_Rank in: ", f)
  }

  debug_cols[[length(debug_cols) + 1]] <- data.table(
    file = f,
    sample = sample_id,
    gene_col = gene_col,
    rank_col = rank_col,
    n_rows = nrow(dt)
  )

  dt[, source_gene_symbol := standardize_symbol(get(gene_col))]
  dt[, is_membrane_gene := source_gene_symbol %in% membrane_symbols]

  dt[, Borda_Rank_numeric := as.numeric(get(rank_col))]
  setorder(dt, Borda_Rank_numeric)

  dt[, epitope_rank_current := seq_len(.N)]

  total_epitopes <- nrow(dt)
  total_membrane <- sum(dt$is_membrane_gene, na.rm = TRUE)
  background_fraction <- total_membrane / total_epitopes

  # Keep example top-ranked membrane hits for inspection.
  annotated_examples[[length(annotated_examples) + 1]] <- dt[
    is_membrane_gene == TRUE,
    .(
      sample = sample_id,
      Gene.Name,
      HLA.Allele,
      Peptide.Length,
      MT.Epitope.Seq,
      Median.MT.IC50.Score,
      Borda_Score,
      Borda_Rank,
      epitope_rank_current
    )
  ][order(epitope_rank_current)][1:min(.N, 50)]

  rows <- list()

  for (k in TOP_N) {
    top <- dt[seq_len(min(k, .N))]
    mem_n <- sum(top$is_membrane_gene, na.rm = TRUE)

    rows[[length(rows) + 1]] <- data.table(
      sample = sample_id,
      cutoff_type = "top_n",
      cutoff = k,
      total_epitopes = total_epitopes,
      total_membrane_epitopes = total_membrane,
      background_membrane_fraction = background_fraction,
      top_epitopes = nrow(top),
      top_membrane_epitopes = mem_n,
      top_membrane_fraction = mem_n / nrow(top),
      fold_enrichment_vs_background = ifelse(
        background_fraction > 0,
        (mem_n / nrow(top)) / background_fraction,
        NA_real_
      )
    )
  }

  for (p in TOP_PROP) {
    k <- max(1, ceiling(p * total_epitopes))
    top <- dt[seq_len(min(k, .N))]
    mem_n <- sum(top$is_membrane_gene, na.rm = TRUE)

    rows[[length(rows) + 1]] <- data.table(
      sample = sample_id,
      cutoff_type = "top_fraction",
      cutoff = p,
      total_epitopes = total_epitopes,
      total_membrane_epitopes = total_membrane,
      background_membrane_fraction = background_fraction,
      top_epitopes = nrow(top),
      top_membrane_epitopes = mem_n,
      top_membrane_fraction = mem_n / nrow(top),
      fold_enrichment_vs_background = ifelse(
        background_fraction > 0,
        (mem_n / nrow(top)) / background_fraction,
        NA_real_
      )
    )
  }

  per_patient[[length(per_patient) + 1]] <- rbindlist(rows)
}

debug_dt <- rbindlist(debug_cols, use.names = TRUE, fill = TRUE)
per_patient_dt <- rbindlist(per_patient, use.names = TRUE, fill = TRUE)
examples_dt <- rbindlist(annotated_examples, use.names = TRUE, fill = TRUE)

summary_dt <- per_patient_dt[, .(
  n_patients = .N,
  median_total_epitopes = median(total_epitopes),
  median_background_membrane_fraction = median(background_membrane_fraction, na.rm = TRUE),

  median_top_membrane_epitopes = median(top_membrane_epitopes, na.rm = TRUE),
  q25_top_membrane_epitopes = quantile(top_membrane_epitopes, 0.25, na.rm = TRUE),
  q75_top_membrane_epitopes = quantile(top_membrane_epitopes, 0.75, na.rm = TRUE),
  min_top_membrane_epitopes = min(top_membrane_epitopes, na.rm = TRUE),
  max_top_membrane_epitopes = max(top_membrane_epitopes, na.rm = TRUE),

  median_top_membrane_fraction = median(top_membrane_fraction, na.rm = TRUE),
  q25_top_membrane_fraction = quantile(top_membrane_fraction, 0.25, na.rm = TRUE),
  q75_top_membrane_fraction = quantile(top_membrane_fraction, 0.75, na.rm = TRUE),

  median_fold_enrichment_vs_background = median(fold_enrichment_vs_background, na.rm = TRUE),
  q25_fold_enrichment_vs_background = quantile(fold_enrichment_vs_background, 0.25, na.rm = TRUE),
  q75_fold_enrichment_vs_background = quantile(fold_enrichment_vs_background, 0.75, na.rm = TRUE)
), by = .(cutoff_type, cutoff)][order(cutoff_type, cutoff)]

display_dt <- summary_dt[, .(
  cutoff_type,
  cutoff,
  n_patients,
  `Background membrane %` = round(100 * median_background_membrane_fraction, 1),
  `Top membrane n median` = round(median_top_membrane_epitopes, 1),
  `Top membrane n IQR` = paste0(
    round(q25_top_membrane_epitopes, 1),
    "-",
    round(q75_top_membrane_epitopes, 1)
  ),
  `Top membrane % median` = round(100 * median_top_membrane_fraction, 1),
  `Top membrane % IQR` = paste0(
    round(100 * q25_top_membrane_fraction, 1),
    "-",
    round(100 * q75_top_membrane_fraction, 1)
  ),
  `Fold enrichment median` = round(median_fold_enrichment_vs_background, 2),
  `Fold enrichment IQR` = paste0(
    round(q25_fold_enrichment_vs_background, 2),
    "-",
    round(q75_fold_enrichment_vs_background, 2)
  )
)]

fwrite(debug_dt, file.path(OUT, "current_final_epitope_column_detection.tsv"), sep = "\t")
fwrite(per_patient_dt, file.path(OUT, "current_ranking_membrane_topk_per_patient.tsv"), sep = "\t")
fwrite(summary_dt, file.path(OUT, "current_ranking_membrane_topk_summary_numeric.tsv"), sep = "\t")
fwrite(display_dt, file.path(OUT, "current_ranking_membrane_topk_summary_display.tsv"), sep = "\t")
fwrite(examples_dt, file.path(OUT, "current_top_ranked_membrane_epitope_examples.tsv"), sep = "\t")

message2("Done. Output: ", OUT)

cat("\nDetected columns:\n")
print(debug_dt)

cat("\nDisplay summary:\n")
print(display_dt)

cat("\nTop membrane examples:\n")
print(head(examples_dt, 20))
