# JABBA / biomass-dynamic backend (stub for future implementation).
#
# Target interface mirrors backend-ss3.R:
#   backend_read_output, backend_scenario_paths, backend_run_retros, etc.
# Collectors under appendices/data/ should consume normalized tables so templates
# stay backend-agnostic once JABBA draws and hindcasts are wired.

backend_capabilities <- function(cfg) {
  list(
    read_output = FALSE,
    retrospectives = TRUE,
    catch_backtests = TRUE,
    likelihood_profiles = FALSE,
    mvln_hessian = FALSE,
    loo_crossvalidation = FALSE,
    length_composition = FALSE,
    ss3diags_runs_test = FALSE
  )
}

jabba_not_implemented <- function(what) {
  stop(
    "JABBA backend '", what, "' is not implemented yet. ",
    "Add backend-jabba.R methods and normalized collectors first.",
    call. = FALSE
  )
}

backend_read_output <- function(path, ...) jabba_not_implemented("read_output")
backend_scenario_paths <- function(cfg, scenarios = NULL) jabba_not_implemented("scenario_paths")
backend_reference_path <- function(cfg) jabba_not_implemented("reference_path")
backend_backtest_params <- function(cfg) jabba_not_implemented("backtest_params")
backend_retro_peel_tags <- function(retro_years) jabba_not_implemented("retro_peel_tags")
backend_retro_dirs_available <- function(cfg) jabba_not_implemented("retro_dirs_available")
backend_run_retros <- function(cfg, ...) jabba_not_implemented("run_retros")
backend_run_forecast_peels <- function(cfg, ...) jabba_not_implemented("run_forecast_peels")
backend_run_backtests <- function(cfg, ...) jabba_not_implemented("run_backtests")
backend_retro_backtest_data <- function(cfg, ...) jabba_not_implemented("retro_backtest_data")
backend_collect_retrospectives <- function(cfg, ...) jabba_not_implemented("collect_retrospectives")
backend_load_internal <- function(cfg, path) jabba_not_implemented("load_internal")
backend_load_external <- function(cfg, path) jabba_not_implemented("load_external")
backend_load_robustness <- function(cfg, retro = NULL) jabba_not_implemented("load_robustness")
backend_load_case_study <- function(cfg) jabba_not_implemented("load_case_study")
backend_load_pe <- function(cfg, paths = NULL, for_plots = FALSE) {
  jabba_not_implemented("load_pe")
}
backend_load_fls <- function(cfg, paths = NULL) jabba_not_implemented("load_fls")
backend_read_forecast_output <- function(path, strict = TRUE) jabba_not_implemented("read_forecast_output")
backend_read_profile_output <- function(run_dir) jabba_not_implemented("read_profile_output")
backend_summarise_run <- function(path, rep = NULL) jabba_not_implemented("summarise_run")
backend_diag_ss <- function(path) jabba_not_implemented("diag_ss")
backend_plot_cpue_residuals <- function(diag_df) jabba_not_implemented("plot_cpue_residuals")
backend_plot_cpue_std_residuals <- function(diag_df) jabba_not_implemented("plot_cpue_std_residuals")
backend_plot_sp_ssb <- function(sp_df) jabba_not_implemented("plot_sp_ssb")
backend_plot_retro_ratios <- function(retro_data, fore_col, ref_col, ylab, title) {
  jabba_not_implemented("plot_retro_ratios")
}
backend_collect_catch_backtests <- function(cfg, scenarios = NULL) {
  jabba_not_implemented("collect_catch_backtests")
}
