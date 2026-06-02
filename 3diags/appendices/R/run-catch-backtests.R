#!/usr/bin/env Rscript

# Unified back-test entry point (forecast-only SS3 passes on existing retro peels).
#
# Workflow (two steps):
#   1. run-retrospectives.R --run     full SS3 re-fit per peel (expensive)
#   2. run-catch-backtests.R --run    projection pass only — no re-estimation
#
# backtest_mode (config-longruns.R, DIAG_BACKTEST_MODE, or --backtest-mode):
#   appendix   — four hindcast scenarios per peel (default; SMA paper appendices)
#   management — reported-catch hindcast to current year, then F=0 and F=FMSY
#                for management_years under retrospectives/<peel>/backtest/
#
# Usage:
#   Rscript appendices/R/run-catch-backtests.R --appendix-dir appendices --run
#   Rscript appendices/R/run-catch-backtests.R --appendix-dir appendices --collect
#   Rscript appendices/R/run-catch-backtests.R --appendix-dir appendices --all
#   Rscript appendices/R/run-catch-backtests.R --appendix-dir appendices --run --backtest-mode management
#
# Prerequisites: converged retro peels with ss3.par (run-retrospectives.R --run).
# Environment: DIAG_SS_BASE, SS3_EXE, DIAG_RETRO_SCENARIOS, DIAG_BACKTEST_MODE

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

  args <- parse_longrun_args(raw_args, defaults = list(setup_on_all = FALSE))
  args$appendix_dir <- appendix_dir
  if (!any(c(args$run, args$collect))) {
    stop("Specify --run, --collect, or --all.", call. = FALSE)
  }

  cfg <- longrun_cfg(appendix_dir, study_id = args$study_id)
  load_appendix_helpers(cfg)
  lr <- longrun_defaults_from_env(longrun_defaults())
  if (!is.null(args$backtest_mode)) {
    lr$backtest_mode <- args$backtest_mode
  }
  mode <- backtest_mode(lr)
  fsc <- backtest_forecast_scenarios(lr)
  scenarios <- args$scenarios %||% cfg$retro_scenarios
  root <- longrun_job_root(cfg, "catch-backtests")

  if (isTRUE(args$run)) {
    write_longrun_status(
      root,
      "running",
      list(scenarios = scenarios, backtest_mode = mode)
    )
    for (sc in scenarios) {
      retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
      if (!dir.exists(retro_root)) {
        stop(
          "Retrospectives not found for ", sc, ": ", retro_root,
          "\nRun run-retrospectives.R --run first.",
          call. = FALSE
        )
      }
      message("Back-test projections for ", sc, " (", retro_root, ")")
      message("  Mode: ", mode)
      if (identical(mode, "management")) {
        message(
          "  Two-stage: reported catch to current year, then F=0 and F=FMSY for ",
          lr$management_years %||% 100L, " yr under ",
          lr$management_out_subdir %||% "backtest", "/"
        )
      } else {
        message(
          "  Forecast-only pass (no re-estimation",
          if (isTRUE(lr$forecast_nohess)) "; SS3 -nohess" else "",
          "); scenarios: ",
          paste(fsc, collapse = ", ")
        )
      }
      if (isTRUE(lr$parallel_forecast_peels)) {
        nw <- args$workers %||% lr$forecast_peel_workers
        message(
          "  Parallel peels: ",
          if (is.null(nw)) "auto (cores - 1)" else nw,
          " workers"
        )
      }
      backend_run_backtests(
        cfg,
        retro_root,
        lr = lr,
        workers = args$workers %||% lr$forecast_peel_workers
      )
    }
    write_longrun_status(
      root,
      "ran",
      list(scenarios = scenarios, backtest_mode = mode)
    )
    message("Back-test SS3 runs complete (mode = ", mode, ").")
  }

  if (isTRUE(args$collect)) {
    write_longrun_status(
      root,
      "collecting",
      list(scenarios = scenarios, backtest_mode = mode)
    )
    out <- collect_catch_backtests(cfg, scenarios = scenarios, lr = lr)
    if (is.null(out) || !NROW(out$summary)) {
      write_longrun_status(
        root,
        "no_results",
        list(scenarios = scenarios, backtest_mode = mode)
      )
      message(
        "No back-test forecast outputs found (mode = ", mode, "). ",
        "Run --run after retrospectives exist, then --collect again."
      )
      return(invisible(NULL))
    }
    write_longrun_status(
      root,
      "complete",
      list(
        n_rows = NROW(out$summary),
        scenarios = scenarios,
        backtest_mode = mode
      )
    )
    message(
      "Collected ", NROW(out$summary), " forecast rows (mode = ", mode, ") to ", root
    )
  }
}

if (sys.nframe() == 0L) {
  main()
}
