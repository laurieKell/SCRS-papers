# Assessment backend dispatch (SS3 today; JABBA and others later).
# Templates and collectors call backend_*() rather than r4ss/sstruth directly.

backend_registry <- function() {
  c("ss3", "jabba")
}

backend_module_path <- function(appendix_dir, backend) {
  file.path(
    normalizePath(appendix_dir, winslash = "/", mustWork = TRUE),
    "R",
    paste0("backend-", backend, ".R")
  )
}

backend_fn <- function(cfg, name) {
  fn <- backend_for(cfg)$env[[name]]
  if (is.null(fn) || !is.function(fn)) {
    stop(
      "Backend '", cfg$backend %||% "ss3", "' has no function '", name, "'.",
      call. = FALSE
    )
  }
  fn
}

load_backend <- function(cfg) {
  backend <- cfg$backend %||% "ss3"
  if (!backend %in% backend_registry()) {
    stop(
      "Unknown backend '", backend, "'. Registered: ",
      paste(backend_registry(), collapse = ", "),
      call. = FALSE
    )
  }
  path <- backend_module_path(cfg$appendix_dir, backend)
  if (!file.exists(path)) {
    stop("Backend module not found: ", path, call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  source(path, local = env)
  list(name = backend, env = env)
}

backend_for <- function(cfg) {
  if (!is.null(cfg$.backend)) {
    return(cfg$.backend)
  }
  cfg$.backend <- load_backend(cfg)
  cfg$.backend
}

backend_capabilities <- function(cfg) {
  backend_fn(cfg, "backend_capabilities")(cfg)
}

backend_read_output <- function(cfg, path, ...) {
  backend_fn(cfg, "backend_read_output")(path, ...)
}

backend_scenario_paths <- function(cfg, scenarios = NULL) {
  backend_fn(cfg, "backend_scenario_paths")(cfg, scenarios = scenarios)
}

backend_reference_path <- function(cfg) {
  backend_fn(cfg, "backend_reference_path")(cfg)
}

backend_backtest_params <- function(cfg) {
  backend_fn(cfg, "backend_backtest_params")(cfg)
}

backend_retro_peel_tags <- function(cfg, retro_years) {
  backend_fn(cfg, "backend_retro_peel_tags")(retro_years)
}

backend_retro_dirs_available <- function(cfg) {
  backend_fn(cfg, "backend_retro_dirs_available")(cfg)
}

backend_run_retros <- function(cfg, dirs, bt = NULL, ...) {
  backend_fn(cfg, "backend_run_retros")(cfg, dirs, bt = bt, ...)
}

backend_run_forecast_peels <- function(cfg, retro_root, lr = NULL, workers = NULL, ...) {
  backend_fn(cfg, "backend_run_forecast_peels")(
    cfg, retro_root, lr = lr, workers = workers, ...
  )
}

backend_run_backtests <- function(cfg, retro_root, lr = NULL, workers = NULL, ...) {
  backend_fn(cfg, "backend_run_backtests")(
    cfg, retro_root, lr = lr, workers = workers, ...
  )
}

backend_retro_backtest_data <- function(cfg, bt = NULL) {
  backend_fn(cfg, "backend_retro_backtest_data")(cfg, bt = bt)
}

backend_collect_retrospectives <- function(cfg, scenarios = NULL, import_from = NULL) {
  backend_fn(cfg, "backend_collect_retrospectives")(
    cfg, scenarios = scenarios, import_from = import_from
  )
}

backend_load_internal <- function(cfg, path) {
  backend_fn(cfg, "backend_load_internal")(cfg, path)
}

backend_load_external <- function(cfg, path) {
  backend_fn(cfg, "backend_load_external")(cfg, path)
}

backend_load_robustness <- function(cfg, retro = NULL) {
  backend_fn(cfg, "backend_load_robustness")(cfg, retro = retro)
}

backend_load_case_study <- function(cfg) {
  backend_fn(cfg, "backend_load_case_study")(cfg)
}

backend_load_pe <- function(cfg, paths = NULL, for_plots = FALSE, parallel = TRUE, workers = NULL) {
  backend_fn(cfg, "backend_load_pe")(
    cfg,
    paths = paths,
    for_plots = for_plots,
    parallel = parallel,
    workers = workers
  )
}

backend_load_fls <- function(cfg, paths = NULL) {
  backend_fn(cfg, "backend_load_fls")(cfg, paths = paths)
}

backend_read_forecast_output <- function(cfg, path, strict = TRUE) {
  backend_fn(cfg, "backend_read_forecast_output")(path, strict = strict)
}

backend_read_profile_output <- function(cfg, run_dir) {
  backend_fn(cfg, "backend_read_profile_output")(run_dir)
}

backend_summarise_run <- function(cfg, path, rep = NULL) {
  backend_fn(cfg, "backend_summarise_run")(path, rep = rep)
}

backend_diag_ss <- function(cfg, path) {
  backend_fn(cfg, "backend_diag_ss")(path)
}

backend_plot_cpue_residuals <- function(cfg, diag_df) {
  backend_fn(cfg, "backend_plot_cpue_residuals")(diag_df)
}

backend_plot_cpue_std_residuals <- function(cfg, diag_df) {
  backend_fn(cfg, "backend_plot_cpue_std_residuals")(diag_df)
}

backend_plot_sp_ssb <- function(cfg, sp_df) {
  backend_fn(cfg, "backend_plot_sp_ssb")(sp_df)
}

backend_plot_retro_ratios <- function(cfg, retro_data, fore_col, ref_col, ylab, title) {
  backend_fn(cfg, "backend_plot_retro_ratios")(
    retro_data, fore_col, ref_col, ylab, title
  )
}

backend_collect_catch_backtests <- function(cfg, scenarios = NULL, lr = NULL) {
  backend_fn(cfg, "backend_collect_catch_backtests")(cfg, scenarios = scenarios, lr = lr)
}

# Facade names for templates (cfg-first API).
assessment_read_output <- function(cfg, path, covar = FALSE, ...) {
  backend_read_output(cfg, path, covar = covar, ...)
}

assessment_reference_report <- function(cfg) {
  backend_read_output(cfg, reference_path(cfg))
}

load_assessment_internal <- function(cfg, path = reference_path(cfg)) {
  backend_load_internal(cfg, path)
}

load_assessment_external <- function(cfg, path = reference_path(cfg)) {
  backend_load_external(cfg, path)
}

load_assessment_case_study <- function(cfg) {
  backend_load_case_study(cfg)
}
