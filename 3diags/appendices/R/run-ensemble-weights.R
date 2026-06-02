#!/usr/bin/env Rscript

# Compute scenario ensemble weights from CV / hindcast skill (CVindexmase, CVindexrho).
# Uses retrospective one-step errors and B/BMSY horizon ratios; LOO CPUE when collected.
# Output: appendices/data/ensemble-weights/
#
# Usage:
#   Rscript appendices/R/run-ensemble-weights.R --appendix-dir appendices --collect
#
# Prerequisites: run-retrospectives.R --collect (or existing retroBacktestData inputs).
# Optional: run-crossvalidation.R --collect for CPUE-based weights.
# Set equal_weights = FALSE in config-mvln.R (or DIAG_MVLN_EQUAL_WEIGHTS=false) to apply
# these weights in the MVLN ensemble table at render time.

main <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- NULL
  i <- 1L
  while (i <= length(raw_args)) {
    if (raw_args[[i]] == "--appendix-dir" && i < length(raw_args)) {
      appendix_dir <- normalizePath(raw_args[[i + 1L]], winslash = "/", mustWork = TRUE)
      break
    }
    i <- i + 1L
  }
  if (is.null(appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  source(file.path(appendix_dir, "R", "longrun-utils.R"))
  source(file.path(appendix_dir, "R", "appendix-utils.R"))
  source(file.path(appendix_dir, "R", "ensemble-cv-weights.R"))

  args <- parse_longrun_args(raw_args, defaults = list(setup_on_all = FALSE))
  if (!isTRUE(args$collect)) {
    stop("Specify --collect (or --all).", call. = FALSE)
  }

  cfg <- longrun_cfg(appendix_dir, study_id = args$study_id)
  root <- longrun_job_root(cfg, "ensemble-weights")
  write_longrun_status(root, "collecting", list())

  out <- collect_ensemble_cv_weights(cfg)
  if (is.null(out) || is.null(out$weights) || !NROW(out$weights)) {
    write_longrun_status(root, "no_results", list())
    stop(
      "No CV metrics available. Run run-retrospectives.R --collect first.",
      call. = FALSE
    )
  }

  write_longrun_status(
    root,
    "complete",
    list(n_models = NROW(out$weights), combine = out$combine_method)
  )
  message("Ensemble CV weights written to ", root)
  if (!is.null(out$weights_combined_mase)) {
    message("CVindexmase: ", paste(round(out$weights_combined_mase, 4), collapse = ", "))
  }
  if (!is.null(out$weights_combined_rho)) {
    message("CVindexrho: ", paste(round(out$weights_combined_rho, 4), collapse = ", "))
  }
  if (!is.null(out$weight_comparison) && NROW(out$weight_comparison)) {
    message("Weight comparison table written (equal / ICES / CVindexmase / CVindexrho).")
  }
}

if (sys.nframe() == 0L) {
  main()
}
