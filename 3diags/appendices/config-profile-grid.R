# Defaults for M × SR_surv_Beta likelihood-profile grids (SMA natl case study).
# Used by appendices/R/run-likelihood-profiles.R and internal appendix plots.

profile_grid_defaults <- function() {
  list(
    M = c(0.5, 0.75, 1, 1.25, 1.5),
    beta = c(0.5, 1, 1.5),
    ctlfile = "control.ss",
    beta_label = "SR_surv_Beta",
    ss_exe = Sys.getenv("SS3_EXE", unset = "ss3"),
    parallel = TRUE,
    workers = NULL,
    skip_finished = TRUE,
    show_in_console = FALSE
  )
}

profile_grid_scenarios <- function(M = NULL, beta = NULL) {
  defaults <- profile_grid_defaults()
  M <- M %||% defaults$M
  beta <- beta %||% defaults$beta
  scen <- expand.grid(M = as.numeric(M), beta = as.numeric(beta), stringsAsFactors = FALSE)
  scen$i <- seq_len(nrow(scen))
  scen
}

profile_data_root <- function(appendix_dir, scenario_dir, data_dir = NULL) {
  base <- data_dir %||% file.path(appendix_dir, "data")
  normalizePath(
    file.path(base, "likelihood-profiles", basename(scenario_dir)),
    winslash = "/",
    mustWork = FALSE
  )
}
