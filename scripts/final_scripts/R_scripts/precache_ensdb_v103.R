#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(AnnotationHub)
  library(AnnotationDbi)
  library(ensembldb)
  library(optparse)
})

# --- CLI ---
opt <- OptionParser()
opt <- add_option(opt, "--sentinel", type="character", help="Path to sentinel file to write")
args <- parse_args(opt)
if (is.null(args$sentinel)) stop("--sentinel is required")

# --- cache path ---
cache <- Sys.getenv("ANNOTATIONHUB_CACHE", unset = "")
if (!nzchar(cache)) stop("ANNOTATIONHUB_CACHE is not set")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

sqlite_path <- file.path(cache, "annotationhub.sqlite3")
index_path  <- file.path(cache, "annotationhub.index.rds")

purge_cache <- function(aggressive = FALSE) {
  # remove only the two fragile files
  if (file.exists(sqlite_path)) try(unlink(sqlite_path, force = TRUE), silent = TRUE)
  if (file.exists(index_path))  try(unlink(index_path,  force = TRUE), silent = TRUE)
  if (aggressive) {
    # nuke everything except the sentinel you are about to write
    keep <- basename(args$sentinel)
    all  <- list.files(cache, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    rmme <- all[basename(all) != keep]
    if (length(rmme)) try(unlink(rmme, recursive = TRUE, force = TRUE), silent = TRUE)
  }
}

open_ah <- function() {
  # Use explicit cache; if this still fails as "Corrupt Cache", caller will purge and retry.
  AnnotationHub(cache = cache)
}

safe_open_ah <- function() {
  # 1st attempt
  ah <- tryCatch(open_ah(), error = identity)
  if (!inherits(ah, "error")) return(ah)

  if (grepl("Corrupt Cache", conditionMessage(ah), ignore.case = TRUE)) {
    # 2nd attempt after light purge
    purge_cache(FALSE)
    ah2 <- tryCatch(open_ah(), error = identity)
    if (!inherits(ah2, "error")) return(ah2)

    if (grepl("Corrupt Cache", conditionMessage(ah2), ignore.case = TRUE)) {
      # 3rd attempt after aggressive purge
      purge_cache(TRUE)
      ah3 <- open_ah()  # let error bubble up if still failing
      return(ah3)
    } else {
      stop(ah2)
    }
  } else {
    stop(ah)
  }
}

# --- do the work ---
ah <- safe_open_ah()

q <- query(ah, "EnsDb.Hsapiens.v103")
if (length(q) < 1L) stop("EnsDb.Hsapiens.v103 not found in AnnotationHub")

edb <- q[[1L]]  # triggers download & caching

# sanity check (TP53)
sym <- AnnotationDbi::mapIds(edb, keys = "ENSG00000141510", keytype = "GENEID", column = "SYMBOL")
if (is.na(sym)) stop("Mapping check failed for TP53")

# sentinel
cat(format(Sys.time()), file = args$sentinel)
