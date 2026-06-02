#!/usr/bin/env Rscript

# Leave-one-out (LOO) CPUE cross-validation for the reference scenario.
# Each run omits one CPUE fleet from data.ss and re-estimates SS3.
# Runs live under appendices/data/crossvalidation/<scenario>/fleet-<id>/.
#
# Usage:
#   Rscript appendices/R/run-crossvalidation.R --appendix-dir appendices --scenario Scenario-1 --setup
#   Rscript appendices/R/run-crossvalidation.R --appendix-dir appendices --scenario Scenario-1 --run
#   Rscript appendices/R/run-crossvalidation.R --appendix-dir appendices --scenario Scenario-1 --collect
#   Rscript appendices/R/run-crossvalidation.R --appendix-dir appendices --scenario Scenario-1 --all
#
# Environment: DIAG_SS_BASE, SS3_EXE, DIAG_CV_FLEETS (comma-separated fleet ids)

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

  args <- parse_longrun_args(raw_args, defaults = list(setup_on_all = TRUE))
  args$appendix_dir <- appendix_dir
  if (is.null(args$scenario)) {
    args$scenario <- "Scenario-1"
  }
  if (!any(c(args$setup, args$run, args$collect))) {
    stop("Specify --setup, --run, --collect, or --all.", call. = FALSE)
  }

  cfg <- longrun_cfg(appendix_dir, study_id = args$study_id)
  cfg$reference_scenario <- args$scenario
  lr <- longrun_defaults_from_env(longrun_defaults())
  src_dir <- reference_path(cfg)
  if (!dir.exists(src_dir)) {
    stop("Reference scenario not found: ", src_dir, call. = FALSE)
  }

  root <- longrun_job_root(cfg, "crossvalidation", args$scenario)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  starter <- if (requireNamespace("r4ss", quietly = TRUE)) {
    r4ss::SS_readstarter(file.path(src_dir, "starter.ss"), verbose = FALSE)
  } else {
    list(datfile = "data.ss")
  }
  data_path <- file.path(src_dir, starter$datfile %||% "data.ss")
  fleets <- cpue_fleets_from_data(data_path)
  if (is.null(fleets) || !NROW(fleets)) {
    stop("No CPUE fleets found in ", data_path, call. = FALSE)
  }
  if (!is.null(lr$crossval_fleets) && length(lr$crossval_fleets)) {
    fleets <- fleets[fleets$fleet_id %in% lr$crossval_fleets, , drop = FALSE]
  }

  grid <- fleets
  grid$dir <- file.path(root, paste0("fleet-", grid$fleet_id))
  utils::write.csv(grid[, c("fleet_id", "fleet_name", "dir")], file.path(root, "grid.csv"), row.names = FALSE)

  if (isTRUE(args$setup) || isTRUE(args$run)) {
    write_longrun_status(root, "setup", list(n_fleets = nrow(grid)))
    files_dir <- file.path(root, "files")
    copy_ss_inputs(src_dir, files_dir)
    for (i in seq_len(nrow(grid))) {
      copy_ss_inputs(files_dir, grid$dir[i])
      patch_data_ss_loo_cpue(file.path(grid$dir[i], basename(data_path)), grid$fleet_id[i])
    }
    message("Setup complete: ", nrow(grid), " LOO scenarios under ", root)
  }

  if (isTRUE(args$run)) {
    write_longrun_status(root, "running", list(n_fleets = nrow(grid)))
    run_one <- function(i) {
      run_ss3_dir(
        grid$dir[i],
        ss_exe = lr$ss_exe,
        skip_finished = lr$skip_finished,
        show_in_console = lr$show_in_console
      )
    }
    results <- parallel_lapply(
      seq_len(nrow(grid)),
      run_one,
      parallel = args$parallel,
      workers = args$workers %||% lr$crossval_workers %||% lr$workers
    )
    n_ok <- sum(vapply(results, function(x) isTRUE(x$ok), logical(1)))
    write_longrun_status(root, "ran", list(n_ok = n_ok, n_total = nrow(grid)))
    message("LOO cross-validation runs: ", n_ok, "/", nrow(grid), " succeeded.")
  }

  if (isTRUE(args$collect)) {
    write_longrun_status(root, "collecting", list())
    out <- collect_crossvalidation(cfg, scenario = args$scenario)
    if (is.null(out) || !NROW(out$results)) {
      write_longrun_status(root, "collect_failed", list())
      stop("No cross-validation results collected under ", root, call. = FALSE)
    }
    write_longrun_status(
      root,
      "complete",
      list(n_results = NROW(out$results), n_fleets = length(unique(out$results$fleet_id)))
    )
    message("Collected ", NROW(out$results), " LOO results to ", root)
  }
}

if (sys.nframe() == 0L) {
  main()
}
