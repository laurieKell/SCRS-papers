#!/usr/bin/env Rscript

# Verify SS3 runs, collected summaries, figures, and appendix content.
#
# Usage:
#   Rscript appendices/R/verify-pipeline.R --appendix-dir appendices

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  study_id <- NULL
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--appendix-dir" && i < length(args)) {
      appendix_dir <- normalizePath(args[[i + 1L]], winslash = "/", mustWork = TRUE)
      i <- i + 2L
    } else if (args[[i]] == "--study" && i < length(args)) {
      study_id <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  if (is.null(appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  source(file.path(appendix_dir, "R", "longrun-utils.R"))
  source(file.path(appendix_dir, "R", "pipeline-utils.R"))
  cfg <- longrun_cfg(appendix_dir, study_id = study_id)
  out <- verify_pipeline(cfg, verbose = TRUE)
  root <- file.path(appendix_dir, "data", "pipeline")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, file.path(root, "verification.rds"))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      out$checks,
      file.path(root, "verification.json"),
      auto_unbox = TRUE,
      pretty = TRUE
    )
  }
  if (!isTRUE(out$ok)) {
    quit(status = 1L)
  }
}

if (sys.nframe() == 0L) {
  main()
}
