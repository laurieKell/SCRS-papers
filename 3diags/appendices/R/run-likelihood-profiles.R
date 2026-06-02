#!/usr/bin/env Rscript

# Setup and run M x SR_surv_Beta likelihood-profile SS3 grids for appendix diagnostics.
# SS3 I/O delegates to sstruth; study paths and status files stay in appendix-utils.
#
# Usage:
#   Rscript appendices/R/run-likelihood-profiles.R --appendix-dir appendices --scenario Scenario-1 --setup
#   Rscript appendices/R/run-likelihood-profiles.R --appendix-dir appendices --scenario Scenario-1 --run
#   Rscript appendices/R/run-likelihood-profiles.R --appendix-dir appendices --scenario Scenario-1 --collect
#   Rscript appendices/R/run-likelihood-profiles.R --appendix-dir appendices --scenario Scenario-1 --all
#
# Environment:
#   DIAG_SS_BASE      SS3 scenario parent directory
#   SS3_EXE           SS3 executable name/path (default ss3)
#   DIAG_PROFILE_M    Comma-separated M multipliers
#   DIAG_PROFILE_BETA Comma-separated SR_surv_Beta values

if (!exists("diag_strict", mode = "function", inherits = FALSE)) {
  for (p in c(
    file.path(getwd(), "appendices", "R", "strict-utils.R"),
    file.path(getwd(), "R", "strict-utils.R")
  )) {
    if (file.exists(p)) {
      source(p, local = FALSE)
      break
    }
  }
}

parse_profile_args <- function(args) {
  out <- list(
    appendix_dir = NULL,
    scenario = "Scenario-1",
    setup = FALSE,
    run = FALSE,
    collect = FALSE,
    all = FALSE,
    parallel = TRUE,
    workers = NULL,
    import_from = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--appendix-dir" && i < length(args)) {
      out$appendix_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--scenario" && i < length(args)) {
      out$scenario <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--workers" && i < length(args)) {
      out$workers <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--setup") {
      out$setup <- TRUE
      i <- i + 1L
    } else if (key == "--run") {
      out$run <- TRUE
      i <- i + 1L
    } else if (key == "--collect") {
      out$collect <- TRUE
      i <- i + 1L
    } else if (key == "--all") {
      out$all <- TRUE
      i <- i + 1L
    } else if (key == "--import-from" && i < length(args)) {
      out$import_from <- args[[i + 1L]]
      out$collect <- TRUE
      i <- i + 2L
    } else if (key == "--no-parallel") {
      out$parallel <- FALSE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  if (isTRUE(out$all)) {
    out$setup <- TRUE
    out$run <- TRUE
    out$collect <- TRUE
  }
  out
}

profile_grid_from_env <- function() {
  defaults <- profile_grid_defaults()
  if (nzchar(Sys.getenv("DIAG_PROFILE_M", unset = ""))) {
    defaults$M <- as.numeric(strsplit(Sys.getenv("DIAG_PROFILE_M"), ",", fixed = TRUE)[[1]])
  }
  if (nzchar(Sys.getenv("DIAG_PROFILE_BETA", unset = ""))) {
    defaults$beta <- as.numeric(strsplit(Sys.getenv("DIAG_PROFILE_BETA"), ",", fixed = TRUE)[[1]])
  }
  if (nzchar(Sys.getenv("SS3_EXE", unset = ""))) {
    defaults$ss_exe <- Sys.getenv("SS3_EXE")
  }
  defaults
}

main <- function() {
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required for likelihood profiles.", call. = FALSE)
  }

  args <- parse_profile_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  if (!any(c(args$setup, args$run, args$collect)) && is.null(args$import_from)) {
    stop("Specify --setup, --run, --collect, --import-from, or --all.", call. = FALSE)
  }

  appendix_dir <- normalizePath(args$appendix_dir, winslash = "/", mustWork = TRUE)
  source(file.path(appendix_dir, "R", "appendix-utils.R"))
  source(file.path(appendix_dir, "R", "longrun-utils.R"), local = FALSE)
  source(file.path(appendix_dir, "config-case-study.R"))
  source(file.path(appendix_dir, "config-profile-grid.R"))

  cfg <- case_study_defaults(appendix_dir = appendix_dir)
  cfg$reference_scenario <- args$scenario
  src_dir <- reference_path(cfg)
  if (!dir.exists(src_dir)) {
    stop("Reference scenario not found: ", src_dir, call. = FALSE)
  }

  root <- profile_root(cfg, args$scenario)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(args$import_from)) {
    write_profile_status(root, "importing", list(from = args$import_from))
    prof <- import_profile_grid(args$import_from, root)
    write_profile_status(
      root,
      "complete",
      list(
        imported_from = args$import_from,
        n_manifest = nrow(prof$manifest),
        n_component_rows = nrow(prof$components)
      )
    )
    message("Imported profile grid from ", args$import_from)
    return(invisible(NULL))
  }

  grid_cfg <- profile_grid_from_env()
  scen <- profile_grid_scenarios(M = grid_cfg$M, beta = grid_cfg$beta)
  scen$dir <- file.path(root, as.character(scen$M), as.character(scen$beta))
  utils::write.csv(scen[, c("i", "M", "beta", "dir")], file.path(root, "grid.csv"), row.names = FALSE)

  files_dir <- file.path(root, "files")
  if (isTRUE(args$setup) || isTRUE(args$run)) {
    write_profile_status(root, "setup", list(n_scenarios = nrow(scen)))
    sstruth::copySsInputs(src_dir, files_dir)
    for (i in seq_len(nrow(scen))) {
      sstruth::setupSs3ProfileScenario(
        src_dir = files_dir,
        dest_dir = scen$dir[i],
        M = scen$M[i],
        beta = scen$beta[i],
        ctlfile = grid_cfg$ctlfile
      )
    }
    message("Setup complete: ", nrow(scen), " scenarios under ", root)
  }

  if (isTRUE(args$run)) {
    write_profile_status(root, "running", list(n_scenarios = nrow(scen)))
    run_one <- function(i) {
      run_profile_scenario(
        scen$dir[i],
        ss_exe = grid_cfg$ss_exe,
        skip_finished = grid_cfg$skip_finished,
        show_in_console = grid_cfg$show_in_console
      )
    }
    results <- parallel_lapply(
      seq_len(nrow(scen)),
      run_one,
      parallel = args$parallel,
      workers = args$workers %||% grid_cfg$workers
    )
    n_ok <- sum(vapply(results, function(x) isTRUE(x$ok), logical(1)))
    write_profile_status(root, "ran", list(n_ok = n_ok, n_total = nrow(scen)))
    message("SS3 runs complete: ", n_ok, "/", nrow(scen), " succeeded.")
  }

  if (isTRUE(args$collect)) {
    prof <- collect_likelihood_profiles(cfg, args$scenario)
    if (is.null(prof)) {
      write_profile_status(root, "collect_failed", list())
      stop("No profile results collected under ", root, call. = FALSE)
    }
    write_profile_status(
      root,
      "complete",
      list(
        n_manifest = nrow(prof$manifest),
        n_component_rows = nrow(prof$components)
      )
    )
    message(
      "Collected ",
      nrow(prof$manifest),
      " runs; component rows: ",
      nrow(prof$components)
    )
  }
}

if (sys.nframe() == 0L) {
  main()
}
