#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

BASE <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"

MODE <- Sys.getenv("MODE", unset = "target")
TARGET_INDEX <- as.integer(Sys.getenv("TARGET_INDEX", unset = "1"))

METRICS_ENV <- Sys.getenv("METRICS_TO_TEST", unset = "degree,strength,betweenness,lci")
METRICS_TO_TEST <- trimws(strsplit(METRICS_ENV, ",", fixed = TRUE)[[1]])
METRICS_TO_TEST <- METRICS_TO_TEST[nzchar(METRICS_TO_TEST)]

METRIC_LABEL <- gsub("[^A-Za-z0-9_.-]+", "_", paste(METRICS_TO_TEST, collapse = "_"))

OUT_BASE_NAME <- Sys.getenv("OUT_BASE_NAME", unset = "reviewer_A_2C_metrics")

SECTION2B_DIR <- file.path(BASE, paste0("reviewer_A_2B_membrane_added_target", TARGET_INDEX))
SECTION2B_MANIFEST <- file.path(SECTION2B_DIR, "section2B_membrane_manifest.tsv")

OUT <- file.path(BASE, paste0(OUT_BASE_NAME, "_", METRIC_LABEL, "_target", TARGET_INDEX))

BASELINE_DEGREE_DIR      <- file.path(BASE, "Network_Metrics_Degree")
BASELINE_STRENGTH_DIR    <- file.path(BASE, "Network_Metrics_Strength")
BASELINE_BETWEENNESS_DIR <- file.path(BASE, "Network_Metrics_Betweenness")
BASELINE_LCI_DIR         <- file.path(BASE, "Network_Metrics_LargestComponentImpact")

TOP_PROPS <- c(0.01, 0.05, 0.10, 0.20)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

sanitize <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

read_baseline_metric <- function(metric, sample_id) {
  if (metric == "degree") {
    f <- file.path(BASELINE_DEGREE_DIR, paste0(sample_id, "_degree.tsv"))
    value_col <- "degree"
  } else if (metric == "strength") {
    f <- file.path(BASELINE_STRENGTH_DIR, paste0(sample_id, "_strength.tsv"))
    value_col <- "strength"
  } else if (metric == "betweenness") {
    f <- file.path(BASELINE_BETWEENNESS_DIR, paste0(sample_id, "_betweenness.tsv"))
    value_col <- "betweenness"
  } else if (metric == "lci") {
    f <- file.path(BASELINE_LCI_DIR, paste0(sample_id, "_impact.tsv"))
    value_col <- "largest_component_impact"
  } else {
    stop("Unknown metric: ", metric)
  }

  if (!file.exists(f)) {
    stop("Missing baseline metric file: ", f)
  }

  dt <- fread(f)

  if (!("gene" %in% names(dt))) {
    stop("Missing gene column in: ", f)
  }

  if (!(value_col %in% names(dt))) {
    stop("Missing value column ", value_col, " in: ", f)
  }

  dt[, .(gene = as.character(gene), value = as.numeric(get(value_col)))]
}

top_genes <- function(dt, prop) {
  x <- copy(dt)
  x <- x[!is.na(value)]

  if (nrow(x) == 0) {
    return(character(0))
  }

  k <- max(1, ceiling(prop * nrow(x)))
  setorder(x, -value, gene)

  x$gene[seq_len(min(k, nrow(x)))]
}

compare_metric <- function(baseline_dt, replicate_dt, metric, target_sample,
                           perturbation_id, cohort_size, background_size) {
  m <- merge(
    baseline_dt,
    replicate_dt,
    by = "gene",
    suffixes = c("_baseline", "_replicate")
  )

  m <- m[!is.na(value_baseline) & !is.na(value_replicate)]

  spearman <- if (nrow(m) >= 3) {
    suppressWarnings(cor(m$value_baseline, m$value_replicate, method = "spearman"))
  } else {
    NA_real_
  }

  kendall <- if (nrow(m) >= 3) {
    suppressWarnings(cor(m$value_baseline, m$value_replicate, method = "kendall"))
  } else {
    NA_real_
  }

  cor_dt <- data.table(
    mode = "jackknife",
    target_sample = target_sample,
    perturbation_id = perturbation_id,
    cohort_size = cohort_size,
    background_size = background_size,
    metric = metric,
    n_common_genes = nrow(m),
    spearman_rank_correlation = spearman,
    kendall_rank_correlation = kendall
  )

  overlap_list <- list()

  for (p in TOP_PROPS) {
    base_top <- top_genes(copy(baseline_dt), p)
    rep_top <- top_genes(copy(replicate_dt), p)

    inter <- length(intersect(base_top, rep_top))
    union <- length(union(base_top, rep_top))

    overlap_list[[as.character(p)]] <- data.table(
      mode = "jackknife",
      target_sample = target_sample,
      perturbation_id = perturbation_id,
      cohort_size = cohort_size,
      background_size = background_size,
      metric = metric,
      top_prop = p,
      baseline_top_n = length(base_top),
      replicate_top_n = length(rep_top),
      intersection_n = inter,
      overlap_fraction_of_baseline_top = ifelse(length(base_top) > 0, inter / length(base_top), NA_real_),
      jaccard = ifelse(union > 0, inter / union, NA_real_)
    )
  }

  list(
    cor = cor_dt,
    overlap = rbindlist(overlap_list, use.names = TRUE, fill = TRUE)
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

  E(graph)$weight <- abs(E(graph)$weight)

  comp <- components(graph)
  largest_cc_idx <- which(comp$membership == which.max(comp$csize))
  baseline_cc_size <- length(largest_cc_idx)

  baseline_cc <- induced_subgraph(graph, largest_cc_idx)
  baseline_weight <- sum(E(baseline_cc)$weight)

  for (v in seq_len(n)) {
    g_temp <- delete_vertices(graph, v)

    if (vcount(g_temp) > 0 && ecount(g_temp) > 0) {
      comp_temp <- components(g_temp)
      new_cc_size <- max(comp_temp$csize)
      new_cc_idx <- which(comp_temp$membership == which.max(comp_temp$csize))
      new_cc <- induced_subgraph(g_temp, new_cc_idx)
      new_weight <- sum(E(new_cc)$weight)

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

compute_metrics <- function(g, metrics) {
  out <- list()
  gnames <- V(g)$name

  w <- E(g)$weight
  if (is.null(w)) {
    w <- rep(1, ecount(g))
  }

  if ("degree" %in% metrics) {
    out$degree <- data.table(
      gene = gnames,
      degree = as.numeric(degree(g, mode = "all"))
    )
  }

  if ("strength" %in% metrics) {
    out$strength <- data.table(
      gene = gnames,
      strength = as.numeric(strength(g, mode = "all", weights = abs(w)))
    )
  }

  if ("betweenness" %in% metrics) {
    w_len <- 1 / pmax(abs(w), .Machine$double.eps)

    out$betweenness <- data.table(
      gene = gnames,
      betweenness = as.numeric(betweenness(g, directed = FALSE, weights = w_len, normalized = TRUE))
    )
  }

  if ("lci" %in% metrics) {
    impact <- largest_component_impact(g)

    out$lci <- data.table(
      gene = names(impact),
      largest_component_impact = as.numeric(impact)
    )
  }

  out
}

metric_to_value_dt <- function(metric, dt) {
  if (metric == "degree") {
    dt[, .(gene = as.character(gene), value = as.numeric(degree))]
  } else if (metric == "strength") {
    dt[, .(gene = as.character(gene), value = as.numeric(strength))]
  } else if (metric == "betweenness") {
    dt[, .(gene = as.character(gene), value = as.numeric(betweenness))]
  } else if (metric == "lci") {
    dt[, .(gene = as.character(gene), value = as.numeric(largest_component_impact))]
  } else {
    stop("Unknown metric: ", metric)
  }
}

if (MODE == "target") {
  if (!file.exists(SECTION2B_MANIFEST)) {
    stop("Missing Section 2B manifest: ", SECTION2B_MANIFEST)
  }

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  for (m in METRICS_TO_TEST) {
    dir.create(file.path(OUT, "metrics", m), recursive = TRUE, showWarnings = FALSE)
  }

  message2("Mode: target")
  message2("Target index: ", TARGET_INDEX)
  message2("Metrics: ", paste(METRICS_TO_TEST, collapse = ", "))
  message2("Output: ", OUT)

  manifest <- fread(SECTION2B_MANIFEST)

  baseline_cache <- list()
  cor_results <- list()
  overlap_results <- list()
  status_list <- list()
  metric_manifest <- list()

  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i]

    perturbation_id <- row$perturbation_id
    safe_id <- sanitize(perturbation_id)
    target_sample <- row$target_sample

    message2("Target ", TARGET_INDEX, " | ", perturbation_id)

    task_start <- Sys.time()

    status_dt <- tryCatch({
      if (!file.exists(row$final_graph_file)) {
        stop("Missing final graph: ", row$final_graph_file)
      }

      g <- readRDS(row$final_graph_file)

      metric_list <- compute_metrics(g, METRICS_TO_TEST)

      for (metric in names(metric_list)) {
        if (metric == "degree") {
          metric_file <- file.path(OUT, "metrics", metric, paste0(safe_id, "_degree.tsv"))
        } else if (metric == "strength") {
          metric_file <- file.path(OUT, "metrics", metric, paste0(safe_id, "_strength.tsv"))
        } else if (metric == "betweenness") {
          metric_file <- file.path(OUT, "metrics", metric, paste0(safe_id, "_betweenness.tsv"))
        } else if (metric == "lci") {
          metric_file <- file.path(OUT, "metrics", metric, paste0(safe_id, "_impact.tsv"))
        }

        fwrite(metric_list[[metric]], metric_file, sep = "\t")

        cache_key <- paste(target_sample, metric, sep = "__")

        if (is.null(baseline_cache[[cache_key]])) {
          baseline_cache[[cache_key]] <- read_baseline_metric(metric, target_sample)
        }

        comp <- compare_metric(
          baseline_dt = baseline_cache[[cache_key]],
          replicate_dt = metric_to_value_dt(metric, metric_list[[metric]]),
          metric = metric,
          target_sample = target_sample,
          perturbation_id = perturbation_id,
          cohort_size = row$cohort_size,
          background_size = row$background_size
        )

        cor_results[[length(cor_results) + 1]] <- comp$cor
        overlap_results[[length(overlap_results) + 1]] <- comp$overlap

        metric_manifest[[length(metric_manifest) + 1]] <- data.table(
          target_index = TARGET_INDEX,
          target_sample = target_sample,
          perturbation_id = perturbation_id,
          metric = metric,
          metric_file = metric_file,
          final_graph_file = row$final_graph_file
        )
      }

      rm(g, metric_list)
      gc(verbose = FALSE)

      data.table(
        target_index = TARGET_INDEX,
        perturbation_id = perturbation_id,
        status = "success",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = NA_character_
      )

    }, error = function(e) {
      data.table(
        target_index = TARGET_INDEX,
        perturbation_id = perturbation_id,
        status = "failed",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = conditionMessage(e)
      )
    })

    status_list[[length(status_list) + 1]] <- status_dt

    if (length(cor_results) > 0) {
      fwrite(rbindlist(cor_results, use.names = TRUE, fill = TRUE),
             file.path(OUT, "stability_rank_correlations.tsv"),
             sep = "\t")
    }

    if (length(overlap_results) > 0) {
      fwrite(rbindlist(overlap_results, use.names = TRUE, fill = TRUE),
             file.path(OUT, "stability_topk_overlap.tsv"),
             sep = "\t")
    }

    if (length(status_list) > 0) {
      fwrite(rbindlist(status_list, use.names = TRUE, fill = TRUE),
             file.path(OUT, "section2C_metric_status.tsv"),
             sep = "\t")
    }

    if (length(metric_manifest) > 0) {
      fwrite(rbindlist(metric_manifest, use.names = TRUE, fill = TRUE),
             file.path(OUT, "section2C_metric_manifest.tsv"),
             sep = "\t")
    }
  }

  message2("Finished target ", TARGET_INDEX)
  quit(save = "no")
}

if (MODE == "merge") {
  MERGE_OUT <- file.path(BASE, "reviewer_A_2C_metrics_merged")
  dir.create(MERGE_OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(MERGE_OUT, "plots"), recursive = TRUE, showWarnings = FALSE)

  dirs <- list.dirs(BASE, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^reviewer_A_2C_metrics_.*_target[0-9]+$", basename(dirs))]

  message2("Merging dirs: ", length(dirs))

  read_all <- function(filename) {
    files <- file.path(dirs, filename)
    files <- files[file.exists(files)]

    if (length(files) == 0) {
      return(data.table())
    }

    rbindlist(lapply(files, function(f) {
      dt <- fread(f)
      dt[, source_folder := basename(dirname(f))]
      dt
    }), use.names = TRUE, fill = TRUE)
  }

  cor_all <- read_all("stability_rank_correlations.tsv")
  topk_all <- read_all("stability_topk_overlap.tsv")
  status_all <- read_all("section2C_metric_status.tsv")
  metric_manifest_all <- read_all("section2C_metric_manifest.tsv")

  fwrite(cor_all, file.path(MERGE_OUT, "filtered_metric_rank_correlations_all.tsv"), sep = "\t")
  fwrite(topk_all, file.path(MERGE_OUT, "filtered_metric_topk_overlap_all.tsv"), sep = "\t")
  fwrite(status_all, file.path(MERGE_OUT, "filtered_metric_status_all.tsv"), sep = "\t")
  fwrite(metric_manifest_all, file.path(MERGE_OUT, "filtered_metric_manifest_all.tsv"), sep = "\t")

  global_rank_summary <- cor_all[, .(
    n_replicates = .N,
    median_spearman = median(spearman_rank_correlation, na.rm = TRUE),
    mean_spearman = mean(spearman_rank_correlation, na.rm = TRUE),
    min_spearman = min(spearman_rank_correlation, na.rm = TRUE),
    q025_spearman = quantile(spearman_rank_correlation, 0.025, na.rm = TRUE),
    q25_spearman = quantile(spearman_rank_correlation, 0.25, na.rm = TRUE),
    q75_spearman = quantile(spearman_rank_correlation, 0.75, na.rm = TRUE),
    q975_spearman = quantile(spearman_rank_correlation, 0.975, na.rm = TRUE),
    max_spearman = max(spearman_rank_correlation, na.rm = TRUE),
    median_kendall = median(kendall_rank_correlation, na.rm = TRUE),
    min_common_genes = min(n_common_genes, na.rm = TRUE),
    max_common_genes = max(n_common_genes, na.rm = TRUE)
  ), by = metric][order(metric)]

  global_topk_summary <- topk_all[, .(
    n_replicates = .N,
    median_overlap_fraction = median(overlap_fraction_of_baseline_top, na.rm = TRUE),
    mean_overlap_fraction = mean(overlap_fraction_of_baseline_top, na.rm = TRUE),
    min_overlap_fraction = min(overlap_fraction_of_baseline_top, na.rm = TRUE),
    q025_overlap_fraction = quantile(overlap_fraction_of_baseline_top, 0.025, na.rm = TRUE),
    q25_overlap_fraction = quantile(overlap_fraction_of_baseline_top, 0.25, na.rm = TRUE),
    q75_overlap_fraction = quantile(overlap_fraction_of_baseline_top, 0.75, na.rm = TRUE),
    q975_overlap_fraction = quantile(overlap_fraction_of_baseline_top, 0.975, na.rm = TRUE),
    max_overlap_fraction = max(overlap_fraction_of_baseline_top, na.rm = TRUE),
    median_jaccard = median(jaccard, na.rm = TRUE),
    min_jaccard = min(jaccard, na.rm = TRUE)
  ), by = .(metric, top_prop)][order(metric, top_prop)]

  per_target_rank_summary <- cor_all[, .(
    n_replicates = .N,
    median_spearman = median(spearman_rank_correlation, na.rm = TRUE),
    min_spearman = min(spearman_rank_correlation, na.rm = TRUE),
    q25_spearman = quantile(spearman_rank_correlation, 0.25, na.rm = TRUE),
    q75_spearman = quantile(spearman_rank_correlation, 0.75, na.rm = TRUE),
    max_spearman = max(spearman_rank_correlation, na.rm = TRUE),
    median_kendall = median(kendall_rank_correlation, na.rm = TRUE)
  ), by = .(target_sample, metric)][order(metric, median_spearman)]

  worst_rank_perturbations <- cor_all[order(spearman_rank_correlation)][1:min(.N, 100)]

  fwrite(global_rank_summary, file.path(MERGE_OUT, "global_filtered_metric_rank_summary.tsv"), sep = "\t")
  fwrite(global_topk_summary, file.path(MERGE_OUT, "global_filtered_metric_topk_summary.tsv"), sep = "\t")
  fwrite(per_target_rank_summary, file.path(MERGE_OUT, "per_target_filtered_metric_rank_summary.tsv"), sep = "\t")
  fwrite(worst_rank_perturbations, file.path(MERGE_OUT, "worst_100_filtered_metric_rank_perturbations.tsv"), sep = "\t")

  pdf(file.path(MERGE_OUT, "plots", "filtered_metric_spearman_boxplots.pdf"), width = 10, height = 6)
  boxplot(
    spearman_rank_correlation ~ metric,
    data = as.data.frame(cor_all),
    las = 2,
    ylab = "Spearman rank correlation",
    main = "Filtered-network metric stability across jackknife perturbations"
  )
  abline(h = 0.90, lty = 2)
  abline(h = 0.80, lty = 3)
  dev.off()

  for (p in sort(unique(topk_all$top_prop))) {
    dtp <- topk_all[top_prop == p]

    pdf(file.path(MERGE_OUT, "plots", paste0("filtered_metric_top_", round(100*p), "pct_overlap_boxplots.pdf")), width = 10, height = 6)
    boxplot(
      overlap_fraction_of_baseline_top ~ metric,
      data = as.data.frame(dtp),
      las = 2,
      ylab = "Fraction of original top genes recovered",
      main = paste0("Filtered-network top ", round(100*p), "% metric overlap")
    )
    abline(h = 0.90, lty = 2)
    abline(h = 0.80, lty = 3)
    dev.off()
  }

  message2("Merged output: ", MERGE_OUT)
  print(global_rank_summary)
  print(global_topk_summary)

  quit(save = "no")
}

stop("Unknown MODE: ", MODE)
