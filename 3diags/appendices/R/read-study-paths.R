#!/usr/bin/env Rscript
# Print study paths for shell resume scripts: ss_base|retro1,retro2|n_retro_years
args <- commandArgs(trailingOnly = TRUE)
appendix_dir <- "appendices"
study_id <- Sys.getenv("DIAG_STUDY", unset = "sma-natl")
i <- 1L
while (i <= length(args)) {
  if (args[[i]] == "--appendix-dir" && i < length(args)) {
    appendix_dir <- args[[i + 1L]]
    i <- i + 2L
  } else if (args[[i]] == "--study" && i < length(args)) {
    study_id <- args[[i + 1L]]
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}
source(file.path(appendix_dir, "R", "study-config.R"))
source(file.path(appendix_dir, "config-case-study.R"))
source(file.path(appendix_dir, "config-longruns.R"))
cfg <- case_study_defaults(appendix_dir, study_id)
lr <- longrun_defaults()
cat(cfg$ss_base, "|", paste(cfg$retro_scenarios, collapse = ","), "|", length(lr$retro_years), sep = "")
