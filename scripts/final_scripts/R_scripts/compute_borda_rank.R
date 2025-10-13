#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

opt <- OptionParser()
opt <- add_option(opt, "--input_csv",      type="character", help="Per-sample annotated_network CSV")
opt <- add_option(opt, "--config",         type="character", help="YAML config with Borda settings")
opt <- add_option(opt, "--out_csv",        type="character", help="Output CSV with Borda_Score and Borda_Rank")
opt <- add_option(opt, "--score_sigfigs",  type="integer",  default = 4,
                  help="Significant figures for Borda_Score in output (default: 4)")
args <- parse_args(opt)

stopifnot(!is.null(args$input_csv), !is.null(args$config), !is.null(args$out_csv))

# ---------- load ----------
dt  <- fread(args$input_csv)
cfg <- yaml::read_yaml(args$config)

borda       <- cfg$borda %||% stop("Missing 'borda' block in YAML.")
cols_cfg    <- borda$columns %||% stop("Missing 'borda.columns' in YAML.")
ties_method <- borda$ties_method %||% "average"
na_policy   <- tolower(borda$na_policy %||% "ignore")
if (!na_policy %in% c("ignore","worst","drop"))
  stop("na_policy must be one of: ignore|worst|drop")

# ---------- helpers ----------
normalize_points <- function(ranks) {
  n <- length(ranks)
  if (n <= 1L) return(rep(1, n))
  (n - ranks) / (n - 1)
}

# If weights don’t sum to 1, renormalize
w <- vapply(cols_cfg, function(x) as.numeric(x$weight %||% NA_real_), numeric(1))
if (any(!is.finite(w))) stop("All columns must have numeric weights in YAML.")
w <- w / sum(w)

dirs <- vapply(cols_cfg, function(x) tolower(x$direction %||% "higher_better"), character(1))
if (any(!dirs %in% c("higher_better","lower_better")))
  stop("direction must be 'higher_better' or 'lower_better'.")

# Allow Network_* <-> Net_* aliasing
aliases <- list(
  "Network_Betweenness" = c("Network_Betweenness","Net_Betweenness"),
  "Network_Degree"      = c("Network_Degree","Net_Degree"),
  "Network_Impact"      = c("Network_Impact","Net_Impact"),
  "Network_Strength"    = c("Network_Strength","Net_Strength"),
  "Network_WCI"         = c("Network_WCI","Net_WCI")
)

resolve_col <- function(name, cols) {
  if (name %in% cols) return(name)
  alts <- aliases[[name]]
  if (!is.null(alts)) {
    hit <- alts[alts %in% cols]
    if (length(hit)) return(hit[1])
  }
  dual <- if (startsWith(name, "Net_")) sub("^Net_", "Network_", name) else sub("^Network_", "Net_", name)
  if (dual %in% cols) return(dual)
  NA_character_
}

N <- nrow(dt)
score    <- numeric(N)
w_avail  <- numeric(N)   # per-row accumulated weight (for na_policy='ignore')
row_keep <- rep(TRUE, N) # used only for na_policy='drop'
used_cols <- character(0)

feat_names <- names(cols_cfg)
for (i in seq_along(feat_names)) {
  feat <- feat_names[i]
  effective_col <- resolve_col(feat, names(dt))
  if (is.na(effective_col)) {
    message(sprintf("[borda] SKIP: column '%s' not found (or alias).", feat))
    next
  }

  x <- dt[[effective_col]]
  suppressWarnings(xn <- as.numeric(x))
  good <- is.finite(xn)

  # For 'drop', track rows that must be removed if any feature is NA
  if (na_policy == "drop") row_keep <- row_keep & good

  sc <- rep(NA_real_, N)
  if (sum(good) > 0L) {
    to_rank <- xn[good]
    rk <- if (dirs[i] == "lower_better") {
      rank(to_rank, ties.method = ties_method)
    } else {
      rank(-to_rank, ties.method = ties_method)
    }
    sc[good] <- normalize_points(rk)
  }

  if (na_policy == "worst") {
    sc[!good] <- 0     # penalize missing
    score   <- score   + w[i] * sc
    w_avail <- w_avail + w[i]       # full weight always counts
  } else if (na_policy == "ignore") {
    # contribute only where we have data; accumulate weight only for those rows
    idx <- which(!is.na(sc))
    if (length(idx)) {
      score[idx]   <- score[idx]   + w[i] * sc[idx]
      w_avail[idx] <- w_avail[idx] + w[i]
    }
  } else { # drop
    # still accumulate now; we will subset rows at the end
    sc[!good] <- NA_real_
    idx <- which(!is.na(sc))
    if (length(idx)) {
      score[idx]   <- score[idx]   + w[i] * sc[idx]
      w_avail[idx] <- w_avail[idx] + w[i]
    }
  }

  used_cols <- c(used_cols, feat)
}

if (na_policy == "ignore") {
  idx <- which(w_avail > 0)
  if (length(idx)) score[idx] <- score[idx] / w_avail[idx]  # renormalize per-row so no penalty
  if (length(idx) < N) score[w_avail == 0] <- 0            # no contributing features → 0
} else if (na_policy == "drop") {
  keep <- which(row_keep)
  if (length(keep) < N) {
    dropped <- N - length(keep)
    message(sprintf("[borda] Dropping %d rows due to missing values (na_policy='drop').", dropped))
  }
  dt     <- dt[row_keep]
  score  <- score[row_keep]
}

# Attach score & rank (rank computed on full-precision scores)
dt[, Borda_Score := score]
dt[, Borda_Rank  := rank(-Borda_Score, ties.method = ties_method)]

# Round score to requested significant figures for output
sig <- as.integer(args$score_sigfigs %||% 4)
dt[, Borda_Score := ifelse(is.na(Borda_Score), NA_real_, signif(Borda_Score, digits = sig))]

fwrite(dt, args$out_csv)
message(sprintf(
  "Wrote %s | na_policy=%s | ties=%s | sigfigs=%d | features used: %s",
  args$out_csv, na_policy, ties_method, sig, paste(used_cols, collapse=", ")
))
