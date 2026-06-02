#!/usr/bin/env Rscript

# Retrospective peels for robustness diagnostics (step 1 of retro → hindcast workflow).
# Full SS3 re-fit per peel; outputs write in-place under <ss_base>/<scenario>/retrospectives/.
# Collected summaries (one-step errors, peel trajectories) land in appendices/data/retrospectives/.
# After this completes, run run-catch-backtests.R for forecast-only hindcasts (no re-fit).
#
# Usage:
#   Rscript appendices/R/run-retrospectives.R --appendix-dir appendices --run
#   Rscript appendices/R/run-retrospectives.R --appendix-dir appendices --collect
#   Rscript appendices/R/run-retrospectives.R --appendix-dir appendices --all
#   Rscript appendices/R/run-retrospectives.R --appendix-dir appendices --collect --import-from
#
# Environment:
#   DIAG_SS_BASE, SS3_EXE, DIAG_RETRO_YEARS, DIAG_RETRO_SCENARIOS

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- NULL
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--appendix-dir" && i < length(args)) {
      appendix_dir <- normalizePath(args[[i + 1L]], winslash = "/", mustWork = TRUE)
      break
    }
    i <- i + 1L
  }
  if (is.null(appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  source(file.path(appendix_dir, "R", "longrun-utils.R"))
  source(file.path(appendix_dir, "R", "appendix-utils.R"))

  args <- parse_longrun_args(
    commandArgs(trailingOnly = TRUE),
    defaults = list(setup_on_all = FALSE)
  )
  args$appendix_dir <- appendix_dir

  if (!any(c(args$run, args$collect))) {
    stop("Specify --run, --collect, or --all.", call. = FALSE)
  }

  cfg <- longrun_cfg(appendix_dir, study_id = args$study_id)
  load_appendix_helpers(cfg)
  lr <- longrun_defaults_from_env(longrun_defaults())
  scenarios <- args$scenarios %||% cfg$retro_scenarios
  root <- longrun_job_root(cfg, "retrospectives")

  if (!is.null(args$import_from)) {
    write_longrun_status(root, "importing", list(from = args$import_from))
    collect_retrospectives(cfg, scenarios = scenarios, import_from = args$import_from)
    write_longrun_status(root, "complete", list(imported_from = args$import_from))
    message("Imported retrospective collection from ", args$import_from)
    return(invisible(NULL))
  }

  if (isTRUE(args$run)) {
    dirs <- scenario_dirs(cfg, scenarios)
    bt <- longrun_params_from_cfg(cfg)
    bt$scenarioDirs <- scenarios
    bt$retroYears <- lr$retro_years
    bt$ssExe <- lr$ss_exe
    bt$retroVerbose <- lr$show_in_console
    bt$retroWorkers <- args$workers %||% lr$retro_workers

    write_longrun_status(root, "running", list(scenarios = scenarios, n_dirs = length(dirs)))
    message("Running retrospectives for: ", paste(scenarios, collapse = ", "))
    message("Peels: ", paste(bt$retroYears, collapse = ", "))

    backend_run_retros(cfg, dirs, bt = bt)
    write_longrun_status(root, "ran", list(scenarios = scenarios))
    message("Retrospective SS3 runs complete under each scenario's retrospectives/ folder.")
  }

  if (isTRUE(args$collect)) {
    write_longrun_status(root, "collecting", list(scenarios = scenarios))
    out <- collect_retrospectives(cfg, scenarios = scenarios)
    if (is.null(out)) {
      write_longrun_status(root, "collect_failed", list())
      stop("No retrospective data collected.", call. = FALSE)
    }
    write_longrun_status(
      root,
      "complete",
      list(
        n_scenarios = length(scenarios),
        n_one_step = NROW(out$data$oneStep),
        n_ratio = NROW(out$data$ratioHorizon)
      )
    )
    message("Collected retrospectives to ", root)
  }
}

if (sys.nframe() == 0L) {
  main()
}
