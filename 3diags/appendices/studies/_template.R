# Copy to studies/<your-stock-id>.R and edit paths, scenarios, and capabilities.
# Select at runtime: DIAG_STUDY=<your-stock-id> or --study <your-stock-id>

study_spec <- function() {
  list(
    study_id = "your-stock-id",
    label = "Your stock (SS3)",
    backend = "ss3",
    assessment_id = "your-stock-ss3",
    model = list(
      class = "SS3",
      version = "3.30"
    ),
    paths = list(
      assessment_base = "C:/path/to/ss3/scenario-parent",
      data_dir = NULL,
      pe_dir = NULL
    ),
    scenarios = list(
      ids = c("Scenario-1", "Scenario-2"),
      reference = "Scenario-1",
      retrospective = c("Scenario-2"),
      register = data.frame(
        id = c("Scenario-1", "Scenario-2"),
        label = c("Reference", "Alternative"),
        purpose = c("Reference assessment", "Robustness scenario"),
        stringsAsFactors = FALSE
      )
    ),
    jobs = list(
      likelihood_profiles = list(
        enabled = FALSE,
        scenario = "Scenario-1",
        import_from = NULL,
        grid = list(
          M = c(0.75, 1, 1.25),
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
      likelihood_profiles = FALSE,
      length_composition = TRUE,
      mvln_hessian = TRUE,
      retrospectives = TRUE,
      catch_backtests = TRUE,
      loo_crossvalidation = FALSE,
      ss3diags_runs_test = TRUE,
      flr_selmat = TRUE
    ),
    paper_figures = character(0)
  )
}
