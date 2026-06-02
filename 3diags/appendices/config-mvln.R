# MVLN ensemble defaults for robustness appendix (best assessment vs scenarios).

mvln_defaults <- function() {
  mc <- 500L
  mc_env <- Sys.getenv("DIAG_MVLN_MC", unset = "")
  if (nzchar(mc_env)) {
    mc <- as.integer(mc_env)
  }
  equal_weights <- TRUE
  eq_env <- Sys.getenv("DIAG_MVLN_EQUAL_WEIGHTS", unset = "")
  if (nzchar(eq_env)) {
    equal_weights <- tolower(eq_env) %in% c("1", "true", "yes")
  }
  weight_method <- Sys.getenv("DIAG_MVLN_WEIGHT_METHOD", unset = "cvindexmase")
  list(
    mc = max(100L, mc),
    seed = 1224L,
    terminal_years = NULL,
    equal_weights = equal_weights,
    weight_method = weight_method,
    weight_combine = "average",
    highlight_reference = TRUE
  )
}
