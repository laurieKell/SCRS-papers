# Defaults for long-running appendix jobs (retrospectives, catch back-tests, LOO cross-validation).
# Override via environment variables or render params (launch_* flags).
#
# backtest_mode:
#   appendix   — four hindcast scenarios per retro peel (SMA paper appendices)
#   management — reported-catch hindcast to current year, then F=0 / F=FMSY projections

longrun_defaults <- function() {
  ss_exe <- Sys.getenv("SS3_EXE", unset = "ss3")
  backtest_mode <- Sys.getenv("DIAG_BACKTEST_MODE", unset = "appendix")
  list(
    backtest_mode = backtest_mode,
    retro_years = 0:-10,
    retro_subdir = "retrospectives",
    forecast_horizons = 1:3,
    forecast_scenarios = c("catchReported", "fMeanFull", "fZero", "fMsy"),
    management_years = 100L,
    management_out_subdir = "backtest",
    management_forecast_scenarios = c("catchReported", "fZero", "fMsy"),
    ss_exe = ss_exe,
    ss_version = "3.30",
    retro_workers = NULL,
    forecast_peel_workers = NULL,
    crossval_workers = NULL,
    crossval_fleets = NULL,
    skip_finished = TRUE,
    parallel = TRUE,
    workers = NULL,
    show_in_console = FALSE,
    skip_existing_forecasts = TRUE,
    parallel_forecast_peels = TRUE,
    forecast_nohess = TRUE
  )
}

longrun_data_root <- function(appendix_dir, job, scenario = NULL, data_dir = NULL) {
  base <- data_dir %||% file.path(appendix_dir, "data")
  if (is.null(scenario)) {
    file.path(base, job)
  } else {
    file.path(base, job, scenario)
  }
}

longrun_params_from_cfg <- function(params) {
  defaults <- longrun_defaults()
  if (file.exists(file.path(params$appendix_dir, "config-longruns.R"))) {
    source(file.path(params$appendix_dir, "config-longruns.R"), local = TRUE)
  }
  list(
    ssBase = params$ss_base,
    scenarioDirs = params$retro_scenarios,
    ssVersion = params$ss_version %||% defaults$ss_version,
    runRetro = FALSE,
    retroYears = defaults$retro_years,
    retroSubdir = defaults$retro_subdir,
    forecastHorizons = defaults$forecast_horizons,
    parallelCollect = isTRUE(defaults$parallel),
    retroWorkers = defaults$retro_workers,
    ssExe = defaults$ss_exe,
    retroVerbose = FALSE,
    forecastScenarios = defaults$forecast_scenarios,
    backtestMode = defaults$backtest_mode,
    managementYears = defaults$management_years,
    managementOutSubdir = defaults$management_out_subdir,
    skipExistingForecastRuns = isTRUE(defaults$skip_existing_forecasts),
    parallelForecastPeels = isTRUE(defaults$parallel_forecast_peels),
    forecastPeelWorkers = defaults$forecast_peel_workers,
    forecastNohess = isTRUE(defaults$forecast_nohess)
  )
}
