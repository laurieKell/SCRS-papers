# Case-study defaults for 3diags appendix templates.
# Study bundles live in appendices/studies/<study_id>.R (see study-config.R).
# Override at render time, e.g. render_appendices(run_code = TRUE, study_id = "sma-2026").

case_study_defaults <- function(appendix_dir = NULL, study_id = NULL, overrides = list()) {
  if (is.null(appendix_dir)) {
    appendix_dir <- normalizePath("appendices", winslash = "/", mustWork = FALSE)
    if (!dir.exists(appendix_dir)) {
      appendix_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    }
  }
  study_path <- file.path(appendix_dir, "R", "study-config.R")
  if (!file.exists(study_path)) {
    stop("Cannot find study-config.R at ", study_path, call. = FALSE)
  }
  source(study_path, local = TRUE)
  build_case_study_config(
    appendix_dir = appendix_dir,
    study_id = study_id,
    overrides = overrides
  )
}

ssRunDir <- function(path) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::ssRunDir(path))
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(path)) {
    return(path)
  }
  if (file.exists(file.path(path, "Report.sso"))) {
    return(path)
  }
  inner <- file.path(path, basename(path))
  if (dir.exists(inner) && file.exists(file.path(inner, "Report.sso"))) {
    return(normalizePath(inner, winslash = "/", mustWork = TRUE))
  }
  path
}

scenario_paths <- function(params) {
  raw <- file.path(params$ss_base, params$scenarios)
  if (!identical(params$backend %||% "ss3", "ss3")) {
    return(raw)
  }
  vapply(raw, ssRunDir, character(1), USE.NAMES = FALSE)
}

reference_path <- function(params) {
  ssRunDir(file.path(params$ss_base, params$reference_scenario))
}

profile_grid_from_cfg <- function(cfg) {
  profile_cfg <- file.path(cfg$appendix_dir, "config-profile-grid.R")
  if (file.exists(profile_cfg)) {
    source(profile_cfg, local = TRUE)
    defaults <- profile_grid_defaults()
  } else {
    defaults <- list(
      M = c(0.5, 0.75, 1, 1.25, 1.5),
      beta = c(0.5, 1, 1.5),
      beta_label = "SR_surv_Beta"
    )
  }
  grid <- cfg$jobs$likelihood_profiles$grid
  if (is.null(grid)) {
    return(defaults)
  }
  utils::modifyList(defaults, grid)
}

production_plot_opts <- function(cfg, which = c("scenarios", "profile")) {
  diags::productionPlotOpts(cfg, which = which)
}

kobe_plot_opts <- function(cfg, which = c("scenarios", "profile")) {
  diags::kobePlotOpts(cfg, which = which)
}
