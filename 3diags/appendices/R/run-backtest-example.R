#!/usr/bin/env Rscript

# Management back-test example (tail-cut retro -> reported catch -> F=0 / F=FMSY).
#
# Generalise by changing study, scenario, peel years, or paths via CLI / env:
#
#   DIAG_STUDY=sma-2026
#   DIAG_SS_BASE=P:/rfmo/iccat/2026/shark/Analysis/ss3
#   DIAG_BACKTEST_SCENARIO=Scenario_2/Scn02_base
#   DIAG_BACKTEST_PEEL_YEARS=0,1,5,10,20
#   DIAG_MANAGEMENT_YEARS=100
#   SS3_EXE=ss3
#
# Usage (from 3diags/):
#   Rscript appendices/R/run-backtest-example.R --appendix-dir appendices --collect
#   Rscript appendices/R/run-backtest-example.R --appendix-dir appendices --run-retro --run-forecast --collect
#   Rscript appendices/R/run-backtest-example.R --appendix-dir appendices --demo --run-forecast --collect
#
# --demo uses sma-natl Scenario-2 under sstruth test data (retros on disk; no P: drive).

main <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- normalizePath(
    Sys.getenv("DIAG_APPENDIX_DIR", unset = file.path(getwd(), "appendices")),
    winslash = "/",
    mustWork = FALSE
  )
  study_id <- Sys.getenv("DIAG_STUDY", unset = "sma-2026")
  scenario_id <- NULL
  peel_years <- NULL
  run_retro <- FALSE
  run_forecast <- FALSE
  collect <- FALSE
  demo <- FALSE
  workers <- NULL

  i <- 1L
  while (i <= length(raw_args)) {
    key <- raw_args[[i]]
    if (key == "--appendix-dir" && i < length(raw_args)) {
      appendix_dir <- normalizePath(raw_args[[i + 1L]], winslash = "/", mustWork = TRUE)
      i <- i + 2L
    } else if (key == "--study" && i < length(raw_args)) {
      study_id <- raw_args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--scenario" && i < length(raw_args)) {
      scenario_id <- raw_args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--peel-years" && i < length(raw_args)) {
      peel_years <- raw_args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--workers" && i < length(raw_args)) {
      workers <- as.integer(raw_args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--run-retro") {
      run_retro <- TRUE
      i <- i + 1L
    } else if (key == "--run-forecast") {
      run_forecast <- TRUE
      i <- i + 1L
    } else if (key == "--collect") {
      collect <- TRUE
      i <- i + 1L
    } else if (key == "--all") {
      run_retro <- TRUE
      run_forecast <- TRUE
      collect <- TRUE
      i <- i + 1L
    } else if (key == "--demo") {
      demo <- TRUE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }

  if (!dir.exists(appendix_dir)) {
    stop("--appendix-dir not found: ", appendix_dir, call. = FALSE)
  }
  if (!any(c(run_retro, run_forecast, collect))) {
    stop("Specify at least one of --run-retro, --run-forecast, --collect, or --all.", call. = FALSE)
  }

  source(file.path(appendix_dir, "R", "backtest-example.R"), local = TRUE)

  cfg <- backtest_example_config(
    appendix_dir = appendix_dir,
    study_id = study_id,
    scenario_id = scenario_id %||% Sys.getenv("DIAG_BACKTEST_SCENARIO", unset = "")
  )
  if (!is.null(peel_years)) {
    cfg$peel_years <- parse_peel_years(peel_years)
  }
  if (isTRUE(demo)) {
    cfg$study_id <- "sma-natl-demo"
    cfg$ss_base <- "C:/active/flr/sstruth/data/ss3/sma-natl"
    cfg$scenario_id <- "Scenario-2"
    cfg$run_dir <- resolve_scenario_dir(cfg$ss_base, cfg$scenario_id)
    cfg$scenario_label <- "Scenario-2 (demo)"
    cfg$peel_years <- c(0, 1, 5, 10)
    cfg$ss_exe <- file.path(cfg$ss_base, "ss3_win.exe")
    message("[backtest] --demo: using ", cfg$run_dir)
  }

  message("[backtest] scenario: ", cfg$scenario_label)
  message("[backtest] run_dir:  ", cfg$run_dir)
  message("[backtest] peels:    ", paste(cfg$peel_years, collapse = ", "))
  message("[backtest] SS3 exe:  ", cfg$ss_exe)

  out <- run_backtest_example(
    cfg,
    run_retro = run_retro,
    run_forecast = run_forecast,
    collect = collect,
    workers = workers
  )

  if (isTRUE(collect) && !is.null(out$root)) {
    message("[backtest] collected to ", out$root)
    if (!is.null(out$forecasts)) {
      message("[backtest] forecast rows: ", NROW(out$forecasts$trajectories))
    }
  }
  invisible(out)
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x

if (sys.nframe() == 0L) {
  main()
}
