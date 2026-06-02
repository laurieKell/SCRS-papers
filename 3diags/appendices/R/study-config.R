# Study bundles: one file per stock under appendices/studies/<study_id>.R
# Generic loading lives in diags; this file sets SCRS/SMA defaults only.

default_study_id <- function() {
  Sys.getenv("DIAG_STUDY", unset = "sma-2026")
}

study_spec_path <- function(study_id, appendix_dir) {
  diags::studySpecPath(study_id, appendix_dir)
}

list_studies <- function(appendix_dir) {
  diags::listStudies(appendix_dir)
}

load_study_spec <- function(study_id = NULL, appendix_dir = NULL) {
  diags::loadStudySpec(
    study_id,
    appendix_dir,
    default_study_id = default_study_id()
  )
}

normalize_study_path <- function(path, project_dir, appendix_dir) {
  diags::normalizeStudyPath(path, project_dir, appendix_dir)
}

study_spec_to_config <- function(spec, appendix_dir, overrides = list()) {
  diags::studySpecToConfig(
    spec,
    appendix_dir,
    overrides = overrides,
    default_assessment_base = "P:/rfmo/iccat/2026/shark/Analysis/ss3"
  )
}

build_case_study_config <- function(appendix_dir = NULL, study_id = NULL, overrides = list()) {
  diags::buildCaseStudyConfig(
    appendix_dir = appendix_dir,
    study_id = study_id,
    overrides = overrides,
    default_study_id = default_study_id(),
    default_assessment_base = "P:/rfmo/iccat/2026/shark/Analysis/ss3"
  )
}

has_capability <- function(cfg, capability) {
  diags::hasCapability(cfg, capability)
}

require_capability <- function(cfg, capability, message = NULL) {
  diags::requireCapability(cfg, capability, message = message)
}
