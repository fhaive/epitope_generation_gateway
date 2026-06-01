#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

BASE <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"

ORIGINAL_FILTERED_DIR <- file.path(
  BASE,
  "Sample_Specific_Networks_PPI_filtered",
  "filtered_networks_rds"
)

OUT <- file.path(BASE, "reviewer_B_edge_jaccard_consistency")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_PROPS <- c(0.01, 0.05, 0.10, 0.20)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

make_edge_key <- function(a, b) {
  paste(pmin(a, b), pmax(a, b), sep = "\r")
}

extract_edges_from_graph <- function(graph_file) {
  g <- readRDS(graph_file)

  if (!inherits(g, "igraph")) {
    stop("Expected igraph object in: ", graph_file)
  }

  if (ecount(g) == 0) {
    return(data.table(
      edge_key = character(),
      from = character(),
      to = character(),
      weight = numeric(),
      abs_weight = numeric()
    ))
  }

  ed <- as.data.table(as_data_frame(g, what = "edges"))

  if (!all(c("from", "to") %in% names(ed))) {
    stop("Could not find from/to columns in graph edges: ", graph_file)
  }

  if (!("weight" %in% names(ed))) {
    ed[, weight := 1.0]
  }

  ed <- ed[, .(
    from = as.character(from),
    to = as.character(to),
    weight = as.numeric(weight)
  )]

  ed <- ed[!is.na(from) & !is.na(to) & from != to]
  ed[, edge_key := make_edge_key(from, to)]
  ed[, abs_weight := abs(weight)]

  setorder(ed, edge_key, -abs_weight)
  ed <- ed[!duplicated(edge_key)]

  ed[, .(edge_key, from, to, weight, abs_weight)]
}

get_top_edge_keys <- function(edge_dt, top_prop) {
  if (nrow(edge_dt) == 0) {
    return(character(0))
  }

  x <- edge_dt[!is.na(abs_weight)]

  if (nrow(x) == 0) {
    return(character(0))
  }

  k <- max(1L, ceiling(top_prop * nrow(x)))
  setorder(x, -abs_weight, edge_key)

  x$edge_key[seq_len(min(k, nrow(x)))]
}

jaccard_stats <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))

  data.table(
    top_n_a = length(a),
    top_n_b = length(b),
    intersection_n = inter,
    union_n = union,
    jaccard = ifelse(union > 0, inter / union, NA_real_),
    overlap_fraction_a = ifelse(length(a) > 0, inter / length(a), NA_real_),
    overlap_fraction_b = ifelse(length(b) > 0, inter / length(b), NA_real_),
    overlap_fraction_min_top = ifelse(
      min(length(a), length(b)) > 0,
      inter / min(length(a), length(b)),
      NA_real_
    )
  )
}

process_same_patient_one <- function(row, original_edges) {
  target_sample <- as.character(row$target_sample)
  perturbation_id <- as.character(row$perturbation_id)
  graph_file <- as.character(row$final_graph_file)

  if (!file.exists(graph_file)) {
    stop("Missing jackknife final graph: ", graph_file)
  }

  if (is.null(original_edges[[target_sample]])) {
    stop("No original edge table for target sample: ", target_sample)
  }

  orig_ed <- original_edges[[target_sample]]
  jack_ed <- extract_edges_from_graph(graph_file)

  out <- vector("list", length(TOP_PROPS))

  for (j in seq_along(TOP_PROPS)) {
    p <- TOP_PROPS[j]

    orig_top <- get_top_edge_keys(orig_ed, p)
    jack_top <- get_top_edge_keys(jack_ed, p)

    js <- jaccard_stats(orig_top, jack_top)

    base_row <- data.table(
      analysis = "same_patient_original_vs_jackknife",
      target_index = as.integer(row$target_index),
      target_sample = target_sample,
      perturbation_id = perturbation_id,
      omitted_sample = as.character(row$omitted_sample),
      top_prop = p,
      original_edge_count = nrow(orig_ed),
      jackknife_edge_count = nrow(jack_ed)
    )

    out[[j]] <- cbind(base_row, js)
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

process_between_patient_one <- function(sample_a, sample_b, original_edges) {
  ed_a <- original_edges[[sample_a]]
  ed_b <- original_edges[[sample_b]]

  out <- vector("list", length(TOP_PROPS))

  for (j in seq_along(TOP_PROPS)) {
    p <- TOP_PROPS[j]

    top_a <- get_top_edge_keys(ed_a, p)
    top_b <- get_top_edge_keys(ed_b, p)

    js <- jaccard_stats(top_a, top_b)

    base_row <- data.table(
      analysis = "between_patient_original_pFINs",
      sample_a = sample_a,
      sample_b = sample_b,
      top_prop = p,
      edge_count_a = nrow(ed_a),
      edge_count_b = nrow(ed_b)
    )

    out[[j]] <- cbind(base_row, js)
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

summarize_jaccard <- function(dt, group_cols) {
  if (nrow(dt) == 0) {
    stop("Cannot summarize empty Jaccard table")
  }

  dt[, .(
    n_comparisons = .N,
    median_jaccard = median(jaccard, na.rm = TRUE),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    min_jaccard = min(jaccard, na.rm = TRUE),
    q025_jaccard = quantile(jaccard, 0.025, na.rm = TRUE),
    q25_jaccard = quantile(jaccard, 0.25, na.rm = TRUE),
    q75_jaccard = quantile(jaccard, 0.75, na.rm = TRUE),
    q975_jaccard = quantile(jaccard, 0.975, na.rm = TRUE),
    max_jaccard = max(jaccard, na.rm = TRUE),
    median_overlap_min_top = median(overlap_fraction_min_top, na.rm = TRUE),
    q025_overlap_min_top = quantile(overlap_fraction_min_top, 0.025, na.rm = TRUE),
    min_overlap_min_top = min(overlap_fraction_min_top, na.rm = TRUE),
    median_intersection_n = median(intersection_n, na.rm = TRUE),
    median_union_n = median(union_n, na.rm = TRUE),
    median_top_n_a = median(top_n_a, na.rm = TRUE),
    median_top_n_b = median(top_n_b, na.rm = TRUE)
  ), by = group_cols]
}

# ------------------------------------------------------------
# Resolve sample map
# ------------------------------------------------------------

message2("Resolving sample map")

section1_dirs <- list.dirs(BASE, recursive = FALSE, full.names = TRUE)
section1_dirs <- section1_dirs[grepl("^reviewer_A_section1_lioness_wci_target[0-9]+$", basename(section1_dirs))]

sample_map_list <- list()

for (d in section1_dirs) {
  man_file <- file.path(d, "section1_manifest.tsv")
  if (!file.exists(man_file)) next

  man <- fread(man_file)

  sample_map_list[[length(sample_map_list) + 1]] <- data.table(
    target_index = as.integer(unique(man$target_index)[1]),
    target_sample = as.character(unique(man$target_sample)[1])
  )
}

sample_map <- unique(rbindlist(sample_map_list, use.names = TRUE, fill = TRUE))
setorder(sample_map, target_index)

message2("Samples resolved: ", nrow(sample_map))

if (nrow(sample_map) != 30) {
  warning("Expected 30 samples, found: ", nrow(sample_map))
}

sample_map[, original_graph_file := file.path(
  ORIGINAL_FILTERED_DIR,
  paste0(target_sample, "_filtered.rds")
)]

missing_original <- sample_map[!file.exists(original_graph_file)]

if (nrow(missing_original) > 0) {
  fwrite(missing_original, file.path(OUT, "missing_original_filtered_graphs.tsv"), sep = "\t")
  stop("Missing original filtered graph files. See missing_original_filtered_graphs.tsv")
}

# ------------------------------------------------------------
# Read original pFIN edges
# ------------------------------------------------------------

message2("Reading original filtered pFIN graphs")

original_edges <- list()

for (i in seq_len(nrow(sample_map))) {
  idx <- sample_map$target_index[i]
  sample_id <- sample_map$target_sample[i]
  f <- sample_map$original_graph_file[i]

  message2("  original target ", idx, " | ", sample_id)

  original_edges[[sample_id]] <- extract_edges_from_graph(f)
}

original_qc <- rbindlist(lapply(names(original_edges), function(s) {
  ed <- original_edges[[s]]

  data.table(
    target_sample = s,
    n_edges = nrow(ed),
    median_abs_weight = median(ed$abs_weight, na.rm = TRUE),
    q25_abs_weight = quantile(ed$abs_weight, 0.25, na.rm = TRUE),
    q75_abs_weight = quantile(ed$abs_weight, 0.75, na.rm = TRUE)
  )
}), use.names = TRUE, fill = TRUE)

fwrite(original_qc, file.path(OUT, "original_pfin_edge_qc.tsv"), sep = "\t")

# ------------------------------------------------------------
# Read jackknife final pFIN manifest
# ------------------------------------------------------------

message2("Reading jackknife final pFIN manifests")

jack_dirs <- list.dirs(BASE, recursive = FALSE, full.names = TRUE)
jack_dirs <- jack_dirs[grepl("^reviewer_A_2B_membrane_added_target[0-9]+$", basename(jack_dirs))]

jack_manifest_files <- file.path(jack_dirs, "section2B_membrane_manifest.tsv")
jack_manifest_files <- jack_manifest_files[file.exists(jack_manifest_files)]

jack_all <- rbindlist(lapply(jack_manifest_files, function(f) {
  dt <- fread(f)
  dt[, source_folder := basename(dirname(f))]
  dt
}), use.names = TRUE, fill = TRUE)

message2("Jackknife final pFIN rows: ", nrow(jack_all))

if (nrow(jack_all) != 870) {
  warning("Expected 870 jackknife final pFIN rows, found: ", nrow(jack_all))
}

# ------------------------------------------------------------
# One-row diagnostic before full run
# ------------------------------------------------------------

message2("Running one-row diagnostic")

diag_row <- jack_all[1]
diag_result <- process_same_patient_one(diag_row, original_edges)

message2("Diagnostic result rows: ", nrow(diag_result))
print(diag_result)

# ------------------------------------------------------------
# B1: same-patient original-vs-jackknife edge Jaccard
# ------------------------------------------------------------

message2("Computing same-patient original-vs-jackknife edge Jaccard")

same_patient_results <- list()
same_patient_status <- vector("list", nrow(jack_all))

for (i in seq_len(nrow(jack_all))) {
  row <- jack_all[i]

  if ((i %% 25) == 0) {
    message2("  same-patient comparison ", i, " / ", nrow(jack_all))
  }

  task_start <- Sys.time()

  res <- tryCatch({
    result_dt <- process_same_patient_one(row, original_edges)

    list(
      result = result_dt,
      status = data.table(
        target_index = as.integer(row$target_index),
        target_sample = as.character(row$target_sample),
        perturbation_id = as.character(row$perturbation_id),
        status = "success",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = NA_character_
      )
    )
  }, error = function(e) {
    list(
      result = data.table(),
      status = data.table(
        target_index = as.integer(row$target_index),
        target_sample = as.character(row$target_sample),
        perturbation_id = as.character(row$perturbation_id),
        status = "failed",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = conditionMessage(e)
      )
    )
  })

  if (nrow(res$result) > 0) {
    same_patient_results[[length(same_patient_results) + 1]] <- res$result
  }

  same_patient_status[[i]] <- res$status

  if ((i %% 100) == 0) {
    tmp_dt <- rbindlist(same_patient_results, use.names = TRUE, fill = TRUE)
    tmp_status <- rbindlist(same_patient_status[seq_len(i)], use.names = TRUE, fill = TRUE)

    if (nrow(tmp_dt) > 0) {
      fwrite(tmp_dt, file.path(OUT, "same_patient_original_vs_jackknife_edge_jaccard_all.partial.tsv"), sep = "\t")
    }

    fwrite(tmp_status, file.path(OUT, "same_patient_edge_jaccard_status.partial.tsv"), sep = "\t")
  }
}

same_patient_dt <- rbindlist(same_patient_results, use.names = TRUE, fill = TRUE)
same_patient_status_dt <- rbindlist(same_patient_status, use.names = TRUE, fill = TRUE)

fwrite(same_patient_dt, file.path(OUT, "same_patient_original_vs_jackknife_edge_jaccard_all.tsv"), sep = "\t")
fwrite(same_patient_status_dt, file.path(OUT, "same_patient_edge_jaccard_status.tsv"), sep = "\t")

message2("Same-patient result rows: ", nrow(same_patient_dt))
message2("Same-patient status summary:")
print(same_patient_status_dt[, .N, by = status])

if (nrow(same_patient_dt) == 0) {
  stop("Same-patient result table is empty; cannot continue")
}

# ------------------------------------------------------------
# B2: between-patient original pFIN edge Jaccard
# ------------------------------------------------------------

message2("Computing between-patient original pFIN edge Jaccard")

samples <- sample_map$target_sample
pairs <- combn(samples, 2, simplify = FALSE)

between_results <- vector("list", length(pairs))

for (i in seq_along(pairs)) {
  pp <- pairs[[i]]

  if ((i %% 25) == 0) {
    message2("  between-patient comparison ", i, " / ", length(pairs))
  }

  between_results[[i]] <- process_between_patient_one(pp[1], pp[2], original_edges)
}

between_dt <- rbindlist(between_results, use.names = TRUE, fill = TRUE)

fwrite(between_dt, file.path(OUT, "between_patient_original_pfin_edge_jaccard_all.tsv"), sep = "\t")

message2("Between-patient result rows: ", nrow(between_dt))

if (nrow(between_dt) == 0) {
  stop("Between-patient result table is empty; cannot continue")
}

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

message2("Summarizing")

same_summary <- summarize_jaccard(same_patient_dt, c("analysis", "top_prop"))
between_summary <- summarize_jaccard(between_dt, c("analysis", "top_prop"))

summary_all <- rbindlist(list(same_summary, between_summary), use.names = TRUE, fill = TRUE)
setorder(summary_all, analysis, top_prop)

fwrite(same_summary, file.path(OUT, "same_patient_edge_jaccard_summary.tsv"), sep = "\t")
fwrite(between_summary, file.path(OUT, "between_patient_edge_jaccard_summary.tsv"), sep = "\t")
fwrite(summary_all, file.path(OUT, "edge_jaccard_summary_all.tsv"), sep = "\t")

display <- summary_all[, .(
  Analysis = fifelse(
    analysis == "same_patient_original_vs_jackknife",
    "Same patient: original vs jackknife",
    "Between patients: original pFINs"
  ),
  `Top edge fraction` = paste0(round(100 * top_prop), "%"),
  `N comparisons` = n_comparisons,
  `Median Jaccard` = round(median_jaccard, 3),
  `2.5% Jaccard` = round(q025_jaccard, 3),
  `Minimum Jaccard` = round(min_jaccard, 3),
  `75% Jaccard` = round(q75_jaccard, 3),
  `97.5% Jaccard` = round(q975_jaccard, 3),
  `Median overlap / smaller top set` = round(100 * median_overlap_min_top, 1),
  `2.5% overlap / smaller top set` = round(100 * q025_overlap_min_top, 1),
  `Minimum overlap / smaller top set` = round(100 * min_overlap_min_top, 1),
  `Median top edges A` = round(median_top_n_a),
  `Median top edges B` = round(median_top_n_b)
)]

fwrite(display, file.path(OUT, "edge_jaccard_supervisor_summary.tsv"), sep = "\t")

md_file <- file.path(OUT, "edge_jaccard_supervisor_summary.md")
con <- file(md_file, open = "wt")

writeLines("| Analysis | Top edge fraction | N comparisons | Median Jaccard | 2.5% Jaccard | Minimum Jaccard | 75% Jaccard | 97.5% Jaccard | Median overlap / smaller top set | 2.5% overlap / smaller top set | Minimum overlap / smaller top set | Median top edges A | Median top edges B |", con)
writeLines("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|", con)

for (i in seq_len(nrow(display))) {
  row <- display[i]

  writeLines(paste0(
    "| ", row$Analysis,
    " | ", row$`Top edge fraction`,
    " | ", row$`N comparisons`,
    " | ", row$`Median Jaccard`,
    " | ", row$`2.5% Jaccard`,
    " | ", row$`Minimum Jaccard`,
    " | ", row$`75% Jaccard`,
    " | ", row$`97.5% Jaccard`,
    " | ", row$`Median overlap / smaller top set`, "%",
    " | ", row$`2.5% overlap / smaller top set`, "%",
    " | ", row$`Minimum overlap / smaller top set`, "%",
    " | ", row$`Median top edges A`,
    " | ", row$`Median top edges B`,
    " |"
  ), con)
}

close(con)

message2("Done. Output: ", OUT)

cat("\nEdge Jaccard supervisor summary:\n\n")
print(display)

cat("\nMarkdown table:\n\n")
cat(readLines(md_file), sep = "\n")
cat("\n")

cat("\nSame-patient status summary:\n")
print(same_patient_status_dt[, .N, by = status])
