# Shortfin mako — SCRS paper scenarios (2026 ICCAT shark WG, SCRS/2026/143).
#
# Canonical SS3 runs for this repo (main paper + appendices):
#   P:/rfmo/iccat/2026/shark/Analysis/ss3/
#     Scenario_1/Sc1_Base, Scenario_2/Scn02_base, Scenario_3/Sc3_Base
# Sensitivity runs (checklist Sens columns; see content/sma/scenario-runs.tex):
#   Scenario_1: Sc1_SensA_InitF0p05, Sc1_SensB_InitF0p005, Sc1_SensC_oldM_newH,
#               Sc1_SensD_includeCondAgeAtLength
#   Scenario_2: Scn02_SensA_InitF_0p05, Scn02_SensB_InitF_0p005, Scn02_SensC_oldM_newH
#   Scenario_3: Sc3_SensA_InitF0p05, Sc3_SensB_InitF0p005, Sc3_SensC_oldM_newH,
#               Sc3_SensD_excludeCondAgeAtLength
# (Report.sso is one level deeper, e.g. Sc1_Base/Sc1_Base — resolved automatically.)
# Override only for testing: DIAG_SS_BASE or --ss-base.
# Collected outputs: appendices/data/sma-2026/

study_spec <- function() {
  list(
    study_id = "sma-2026",
    label = "Shortfin mako (2026 ICCAT shark WG)",
    backend = "ss3",
    assessment_id = "sma-ss3-2026",
    wg_document = "SCRS/2026/143",
    model = list(
      class = "SS3",
      version = "3.30"
    ),
    paths = list(
      assessment_base = "P:/rfmo/iccat/2026/shark/Analysis/ss3",
      data_dir = "data/sma-2026",
      pe_dir = "../pe-ss3-sma"
    ),
    scenarios = list(
      ids = c(
        "Scenario_1/Sc1_Base",
        "Scenario_2/Scn02_base",
        "Scenario_3/Sc3_Base"
      ),
      reference = "Scenario_1/Sc1_Base",
      retrospective = c("Scenario_2/Scn02_base", "Scenario_3/Sc3_Base"),
      register = data.frame(
        id = c(
          "Scenario_1/Sc1_Base",
          "Scenario_2/Scn02_base",
          "Scenario_3/Sc3_Base"
        ),
        label = c("Scenario 1", "Scenario 2", "Scenario 3"),
        purpose = c(
          "Reference assessment (2026 shark WG)",
          "Alternative structural assumption 1 (2026 shark WG)",
          "Alternative structural assumption 2 (2026 shark WG)"
        ),
        stringsAsFactors = FALSE
      ),
      sensitivity_register = data.frame(
        scenario = c(
          rep("Scenario 1", 4L),
          rep("Scenario 2", 3L),
          rep("Scenario 3", 4L)
        ),
        id = c(
          "Scenario_1/Sc1_SensA_InitF0p05",
          "Scenario_1/Sc1_SensB_InitF0p005",
          "Scenario_1/Sc1_SensC_oldM_newH",
          "Scenario_1/Sc1_SensD_includeCondAgeAtLength",
          "Scenario_2/Scn02_SensA_InitF_0p05",
          "Scenario_2/Scn02_SensB_InitF_0p005",
          "Scenario_2/Scn02_SensC_oldM_newH",
          "Scenario_3/Sc3_SensA_InitF0p05",
          "Scenario_3/Sc3_SensB_InitF0p005",
          "Scenario_3/Sc3_SensC_oldM_newH",
          "Scenario_3/Sc3_SensD_excludeCondAgeAtLength"
        ),
        label = c(
          "Scenario 1 sens A (Init F 0.05)",
          "Scenario 1 sens B (Init F 0.005)",
          "Scenario 1 sens C (old M, new H)",
          "Scenario 1 sens D (include CAAL growth)",
          "Scenario 2 sens A (Init F 0.05)",
          "Scenario 2 sens B (Init F 0.005)",
          "Scenario 2 sens C (old M, new H)",
          "Scenario 3 sens A (Init F 0.05)",
          "Scenario 3 sens B (Init F 0.005)",
          "Scenario 3 sens C (old M, new H)",
          "Scenario 3 sens D (exclude CAAL growth)"
        ),
        purpose = c(
          rep("Alternative initial F prior", 2L),
          "Alternative M and steepness (modified PW)",
          "CAAL growth sensitivity (include CAAL)",
          "Alternative initial F prior",
          "Alternative initial F prior",
          "Alternative M and steepness (modified PW)",
          rep("Alternative initial F prior", 2L),
          "Alternative M and steepness (modified PW)",
          "CAAL growth sensitivity (exclude CAAL)"
        ),
        stringsAsFactors = FALSE
      )
    ),
    jobs = list(
      likelihood_profiles = list(
        enabled = TRUE,
        scenario = "Scenario_1/Sc1_Base",
        import_from = NULL,
        grid = list(
          M = c(0.5, 0.75, 1, 1.25, 1.5),
          beta = c(0.5, 1, 1.5),
          beta_label = "SR_surv_Beta"
        )
      ),
      crossvalidation = list(
        enabled = TRUE,
        scenario = "Scenario_1/Sc1_Base"
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
        scenarios = list(
          width = 12,
          height = 4,
          show_trajectories = TRUE,
          trajectory_yield_scale = 0.25
        ),
        profile = list(
          ylim = c(0, 1500),
          show_trajectories = TRUE,
          trajectory_yield_scale = 0.25
        )
      ),
      kobe = list(
        scenarios = list(
          width = 12,
          height = 5,
          xlim = c(0, 2.5),
          ylim = c(0, 11),
          facet_scales = "fixed",
          show_trajectories = TRUE,
          show_equilibrium = TRUE,
          trajectory_linewidth = 0.9,
          trajectory_colour = "grey20"
        ),
        profile = list(
          xlim = c(0, 2.5),
          ylim = c(0, 2.5),
          show_trajectories = TRUE,
          show_equilibrium = TRUE
        )
      )
    )
  )
}
