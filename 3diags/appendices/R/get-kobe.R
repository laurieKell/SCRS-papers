#!/usr/bin/env Rscript
#
# Extract SS_output for all discovered runs, then return Kobe as one data.frame.
#
# Usage:
#   Rscript appendices/R/get-kobe.R --workers 4
#   Rscript appendices/R/get-kobe.R --from-cache --out appendices/data/sma-2026/kobe-all.csv

args <- commandArgs(trailingOnly = TRUE)
collectOnly <- "--collect-only" %in% args
fromCache <- "--from-cache" %in% args
refresh <- "--refresh" %in% args
parallel <- !("--no-parallel" %in% args)
workers <- NULL
wi <- match("--workers", args)
if (!is.na(wi) && wi < length(args)) {
  workers <- as.integer(args[[wi + 1L]])
}

out <- NULL
i <- match("--out", args)
if (!is.na(i) && i < length(args)) {
  out <- args[[i + 1L]]
}

if (!requireNamespace("sstruth", quietly = TRUE)) {
  stop("Package 'sstruth' is required.", call. = FALSE)
}

appendixDir <- normalizePath("appendices", winslash = "/", mustWork = FALSE)
study <- sstruth::ssStudyRuns(
  studyId = "sma-2026",
  appendixDir = appendixDir,
  discover = TRUE,
  includeSensitivity = TRUE
)

message("[get-kobe] ", nrow(study$runs), " runs under ", study$base)

if (!fromCache) {
  sstruth::SS_outputs(
    study$runs,
    refresh = refresh,
    parallel = parallel,
    workers = workers
  )
  if (collectOnly) {
    message("[get-kobe] collect complete (ss_output.rds in each run folder)")
    quit(save = "no", status = 0)
  }
}

kobe <- sstruth::ssKobe(
  study$runs,
  cache = TRUE,
  parallel = parallel,
  workers = workers
)

if (!is.null(out) && nzchar(out)) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(kobe, out, row.names = FALSE)
  message("[get-kobe] wrote ", nrow(kobe), " rows -> ", out)
} else {
  print(utils::head(kobe, 10))
  message("[get-kobe] ", nrow(kobe), " rows x ", ncol(kobe), " cols")
}

invisible(kobe)
