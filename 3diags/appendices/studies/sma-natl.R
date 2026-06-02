# Shortfin mako national SS3 case study — 2025 SMA local copy.
# For SCRS 2026 WG runs use sma-2026.R. For the journal paper use flrpapers/3strands (sma-2025).
study_spec <- function() {
  list(
    study_id = "sma-natl",
    label = "Shortfin mako (national SS3)",
    backend = "ss3",
    assessment_id = "sma-ss3",
    model = list(
      class = "SS3",
      version = "3.30"
    ),
    paths = list(
      assessment_base = "C:/active/flr/sstruth/data/ss3/sma-natl",
      data_dir = NULL,
      pe_dir = "../pe-ss3-sma"
    ),
    scenarios = list(
      ids = c("Scenario-1", "Scenario-2", "Scenario-3"),
      reference = "Scenario-1",
      retrospective = c("Scenario-2", "Scenario-3"),
      register = data.frame(
        id = c("Scenario-1", "Scenario-2", "Scenario-3"),
        label = c("Reference", "Alternative 1", "Alternative 2"),
        purpose = c(
          "Reference assessment",
          "Alternative structural assumption 1",
          "Alternative structural assumption 2"
        ),
        stringsAsFactors = FALSE
      )
    ),
    jobs = list(
      likelihood_profiles = list(
        enabled = TRUE,
        scenario = "Scenario-1",
        import_from = "C:/active/scrs-papers/ss-ll/data/ICCAT_2017_Model_Run_3_SSv330.23.02/00_Model_run_files",
        grid = list(
          M = c(0.5, 0.75, 1, 1.25, 1.5),
          beta = c(0.5, 1, 1.5),
          beta_label = "SR_surv_Beta"
        )
      ),
      crossvalidation = list(
        enabled = FALSE,
        scenario = "Scenario-1"
      )
    ),
    capabilities = list(
      likelihood_profiles = TRUE,
      length_composition = TRUE,
      mvln_hessian = TRUE,
      retrospectives = TRUE,
      catch_backtests = TRUE,
      loo_crossvalidation = TRUE,
      ss3diags_runs_test = TRUE,
      flr_selmat = TRUE
    ),
    paper_figures = c(
      "sma-fig-ts-1.png",
      "sma-fig-pe-phase-1.png",
      "selmat.png",
      "sma-fig-res-time-1.png",
      "sma-fig-res-acf-1.png",
      "4panel.png"
    ),
    plots = list(
      production_functions = list(
        profile = list(
          ylim = c(0, 1500),
          show_trajectories = TRUE,
          trajectory_yield_scale = 0.25
        )
      ),
      kobe = list(
        scenarios = list(xlim = c(0, 2.5), ylim = c(0, 2.5), show_trajectories = TRUE),
        profile = list(xlim = c(0, 2.5), ylim = c(0, 2.5), show_trajectories = TRUE)
      )
    )
  )
}
