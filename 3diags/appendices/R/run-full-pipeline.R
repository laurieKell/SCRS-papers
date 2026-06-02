#!/usr/bin/env Rscript

# End-to-end 3diags analysis pipeline (retros → hindcasts → ensemble → profiles → render → verify).
#
# Usage:
#   Rscript appendices/R/run-full-pipeline.R --appendix-dir appendices --fresh --fresh-retros
#   Rscript appendices/R/run-full-pipeline.R --appendix-dir appendices --verify-only
#   Rscript appendices/R/run-full-pipeline.R --appendix-dir appendices --collect-render
#
# Flags:
#   --fresh           Clear collected data + hindcast forecast folders before run
#   --fresh-retros    Also delete retrospective peel directories (full SS3 re-fit)
#   --skip-retros     Skip retrospective SS3 runs (use existing peels)
#   --skip-hindcasts  Skip catch back-test SS3 runs
#   --include-cv      Run LOO CPUE cross-validation (Scenario-1; slow)
#   --import-profiles Import likelihood profile grid from ss-ll (default)
#   --run-profiles    Run local M×beta profile grid instead of import
#   --collect-render  Only collect, render appendices, verify (no SS3)
#   --verify-only     Verification checks only
#
# Environment: SS3_EXE, DIAG_SS_BASE, DIAG_SKIP_EXISTING_FORECASTS=false (forced when --fresh)

parse_pipeline_args <- function(args) {
  out <- list(
    appendix_dir = NULL,
    study_id = NULL,
    fresh = FALSE,
    fresh_retros = FALSE,
    skip_retros = FALSE,
    skip_hindcasts = FALSE,
    include_cv = FALSE,
    import_profiles = TRUE,
    run_profiles = FALSE,
    collect_render = FALSE,
    verify_only = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--appendix-dir" && i < length(args)) {
      out$appendix_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--study" && i < length(args)) {
      out$study_id <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--fresh") {
      out$fresh <- TRUE
      i <- i + 1L
    } else if (key == "--fresh-retros") {
      out$fresh_retros <- TRUE
      i <- i + 1L
    } else if (key == "--skip-retros") {
      out$skip_retros <- TRUE
      i <- i + 1L
    } else if (key == "--skip-hindcasts") {
      out$skip_hindcasts <- TRUE
      i <- i + 1L
    } else if (key == "--include-cv") {
      out$include_cv <- TRUE
      i <- i + 1L
    } else if (key == "--import-profiles") {
      out$import_profiles <- TRUE
      out$run_profiles <- FALSE
      i <- i + 1L
    } else if (key == "--run-profiles") {
      out$run_profiles <- TRUE
      out$import_profiles <- FALSE
      i <- i + 1L
    } else if (key == "--collect-render") {
      out$collect_render <- TRUE
      i <- i + 1L
    } else if (key == "--verify-only") {
      out$verify_only <- TRUE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}

run_rscript <- function(script, args, label, wd = getwd()) {
  message("\n=== ", label, " ===")
  cmd <- paste(c("Rscript", shQuote(script), args), collapse = " ")
  message("wd: ", wd)
  message(cmd)
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(wd)
  status <- system(cmd)
  if (!identical(status, 0L) && !identical(status, 0)) {
    stop(label, " failed (exit ", status, ")", call. = FALSE)
  }
  invisible(status)
}

main <- function() {
  args <- parse_pipeline_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  appendix_dir <- normalizePath(args$appendix_dir, winslash = "/", mustWork = TRUE)

  source(file.path(appendix_dir, "R", "longrun-utils.R"))
  source(file.path(appendix_dir, "R", "pipeline-utils.R"))
  source(file.path(appendix_dir, "R", "study-config.R"))
  cfg <- longrun_cfg(appendix_dir, study_id = args$study_id)
  project_dir <- cfg$project_dir
  study_flag <- if (!is.null(args$study_id)) {
    paste("--study", shQuote(args$study_id, type = "cmd"))
  } else {
    ""
  }
  ad <- paste(
    "--appendix-dir", shQuote(appendix_dir, type = "cmd"),
    study_flag
  )
  prof_scenario <- cfg$profile_scenario %||% cfg$reference_scenario
  cv_scenario <- cfg$cv_scenario %||% cfg$reference_scenario

  if (isTRUE(args$verify_only)) {
    verify_pipeline(cfg, verbose = TRUE)
    return(invisible(NULL))
  }

  if (isTRUE(args$fresh) || isTRUE(args$fresh_retros)) {
    message("=== Fresh reset ===")
    if (isTRUE(args$fresh)) {
      rm <- clear_collected_data(appendix_dir, data_dir = cfg$data_dir)
      nh <- clear_hindcast_forecast_dirs(cfg)
      message("Cleared collected data: ", paste(names(rm), collapse = ", "))
      message("Removed hindcast forecast dirs: ", nh)
      Sys.setenv(DIAG_SKIP_EXISTING_FORECASTS = "false")
    }
    if (isTRUE(args$fresh_retros)) {
      nr <- clear_retrospective_peels(cfg)
      message("Removed retrospective roots: ", nr)
    }
    write_pipeline_status(appendix_dir, "fresh", "done")
  }

  if (!isTRUE(args$collect_render)) {
    if (!isTRUE(args$skip_retros)) {
      write_pipeline_status(appendix_dir, "retrospectives", "running")
      run_rscript(
        file.path(appendix_dir, "R", "run-retrospectives.R"),
        paste(ad, "--all"),
        "Retrospectives (run + collect)",
        wd = project_dir
      )
      write_pipeline_status(appendix_dir, "retrospectives", "complete")
    } else {
      message("Skipping retrospective SS3 runs (--skip-retros)")
    }

    if (!isTRUE(args$skip_hindcasts)) {
      write_pipeline_status(appendix_dir, "catch-backtests", "running")
      run_rscript(
        file.path(appendix_dir, "R", "run-catch-backtests.R"),
        paste(ad, "--all"),
        "Catch back-tests (run + collect)",
        wd = project_dir
      )
      write_pipeline_status(appendix_dir, "catch-backtests", "complete")
    } else {
      message("Skipping catch back-tests (--skip-hindcasts)")
    }

    write_pipeline_status(appendix_dir, "ensemble-weights", "running")
    run_rscript(
      file.path(appendix_dir, "R", "run-ensemble-weights.R"),
      paste(ad, "--collect"),
      "Ensemble weights (collect)",
      wd = project_dir
    )
    write_pipeline_status(appendix_dir, "ensemble-weights", "complete")

    if (isTRUE(args$run_profiles)) {
      write_pipeline_status(appendix_dir, "likelihood-profiles", "running")
      run_rscript(
        file.path(appendix_dir, "R", "run-likelihood-profiles.R"),
        paste(ad, "--scenario", prof_scenario, "--all"),
        "Likelihood profiles (setup + run + collect)",
        wd = project_dir
      )
    } else if (has_capability(cfg, "likelihood_profiles")) {
      import_from <- cfg$profile_import_from
      if (!is.null(import_from) && nzchar(import_from) && dir.exists(import_from)) {
        write_pipeline_status(appendix_dir, "likelihood-profiles", "importing")
        run_rscript(
          file.path(appendix_dir, "R", "run-likelihood-profiles.R"),
          paste(
            ad,
            "--scenario", prof_scenario,
            "--import-from",
            shQuote(import_from, type = "cmd")
          ),
          "Likelihood profiles (import + collect)",
          wd = project_dir
        )
      } else {
        stop("Profile import path not found for study ", cfg$study_id, call. = FALSE)
      }
    } else {
      message("Skipping likelihood profiles (disabled for study ", cfg$study_id, ")")
    }
    write_pipeline_status(appendix_dir, "likelihood-profiles", "complete")

    if (isTRUE(args$include_cv)) {
      write_pipeline_status(appendix_dir, "crossvalidation", "running")
      run_rscript(
        file.path(appendix_dir, "R", "run-crossvalidation.R"),
        paste(ad, "--scenario", cv_scenario, "--all"),
        "LOO CPUE cross-validation",
        wd = project_dir
      )
      write_pipeline_status(appendix_dir, "crossvalidation", "complete")
    } else {
      message("Skipping LOO cross-validation (pass --include-cv to enable)")
    }
  }

  write_pipeline_status(appendix_dir, "render", "running")
  run_rscript(
    file.path(appendix_dir, "render_appendices.R"),
    paste("--run-code", study_flag),
    "Render appendices",
    wd = project_dir
  )
  write_pipeline_status(appendix_dir, "render", "complete")

  out <- verify_pipeline(cfg, verbose = TRUE)
  root <- file.path(appendix_dir, "data", "pipeline")
  saveRDS(out, file.path(root, "verification.rds"))
  write_pipeline_status(
    appendix_dir,
    "pipeline",
    if (isTRUE(out$ok)) "complete" else "failed",
    extra = out$checks
  )
  if (!isTRUE(out$ok)) {
    quit(status = 1L)
  }
}

if (sys.nframe() == 0L) {
  main()
}
