# Fresh-reset helpers and post-run verification for the 3diags appendix pipeline.

forecast_scenario_dirs <- function(retro_root, peel_tags, forecast_scenarios, lr = NULL) {
  out <- character()
  for (peel in peel_tags) {
    for (fsc in forecast_scenarios) {
      if (is.null(lr)) {
        out <- c(out, file.path(retro_root, peel, fsc))
      } else if (exists("forecast_run_dir", mode = "function")) {
        out <- c(out, forecast_run_dir(retro_root, peel, fsc, lr))
      } else {
        out <- c(out, file.path(retro_root, peel, fsc))
      }
    }
  }
  out
}

retro_peel_tags <- function(retro_years, cfg = NULL) {
  if (!is.null(cfg)) {
    return(backend_retro_peel_tags(cfg, retro_years))
  }
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    vapply(retro_years, function(y) {
      if (y == 0L) {
        "retro0"
      } else {
        paste0("retro", y)
      }
    }, character(1))
  } else {
    sstruth::retroTag(as.integer(retro_years))
  }
}

clear_collected_data <- function(appendix_dir, jobs = c(
  "retrospectives",
  "catch-backtests",
  "ensemble-weights",
  "likelihood-profiles"
), data_dir = NULL) {
  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  data_root <- data_dir %||% file.path(appendix_dir, "data")
  removed <- list()
  for (job in jobs) {
    root <- file.path(data_root, job)
    if (!dir.exists(root)) {
      next
    }
    fl <- list.files(root, full.names = TRUE, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    fl <- fl[!grepl("[/\\\\]\\.", fl)]
    if (length(fl)) {
      unlink(fl, recursive = TRUE, force = TRUE)
      removed[[job]] <- length(fl)
    }
  }
  cv_root <- file.path(data_root, "crossvalidation")
  if (dir.exists(cv_root)) {
    for (sc_dir in list.dirs(cv_root, recursive = FALSE, full.names = TRUE)) {
      for (nm in c("results.csv", "results.rds", "status.json", "status.rds", "manifest.csv")) {
        f <- file.path(sc_dir, nm)
        if (file.exists(f)) {
          unlink(f)
          removed[[basename(sc_dir)]] <- (removed[[basename(sc_dir)]] %||% 0L) + 1L
        }
      }
    }
  }
  invisible(removed)
}

clear_hindcast_forecast_dirs <- function(cfg, scenarios = NULL, retro_years = NULL, forecast_scenarios = NULL, lr = NULL) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  scenarios <- scenarios %||% cfg$retro_scenarios
  retro_years <- retro_years %||% lr$retro_years
  forecast_scenarios <- forecast_scenarios %||% backtest_forecast_scenarios(lr)
  peel_tags <- retro_peel_tags(retro_years)
  n_removed <- 0L
  for (sc in scenarios) {
    retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
    if (!dir.exists(retro_root)) {
      next
    }
    dirs <- forecast_scenario_dirs(retro_root, peel_tags, forecast_scenarios, lr = lr)
    for (d in dirs) {
      if (dir.exists(d)) {
        unlink(d, recursive = TRUE, force = TRUE)
        n_removed <- n_removed + 1L
      }
    }
    work <- file.path(retro_root, ".ss3_projection_work")
    if (dir.exists(work)) {
      unlink(work, recursive = TRUE, force = TRUE)
    }
  }
  invisible(n_removed)
}

clear_retrospective_peels <- function(cfg, scenarios = NULL) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  lr <- longrun_defaults()
  scenarios <- scenarios %||% cfg$retro_scenarios
  n_removed <- 0L
  for (sc in scenarios) {
    retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
    if (!dir.exists(retro_root)) {
      next
    }
    unlink(retro_root, recursive = TRUE, force = TRUE)
    n_removed <- n_removed + 1L
  }
  invisible(n_removed)
}

count_forecast_reports <- function(cfg, scenarios = NULL, retro_years = NULL, forecast_scenarios = NULL, lr = NULL) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  scenarios <- scenarios %||% cfg$retro_scenarios
  retro_years <- retro_years %||% lr$retro_years
  forecast_scenarios <- forecast_scenarios %||% backtest_forecast_scenarios(lr)
  peel_tags <- retro_peel_tags(retro_years)
  rows <- list()
  for (sc in scenarios) {
    retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
    for (peel in peel_tags) {
      for (fsc in forecast_scenarios) {
        run_dir <- if (exists("forecast_run_dir", mode = "function")) {
          forecast_run_dir(retro_root, peel, fsc, lr)
        } else {
          file.path(retro_root, peel, fsc)
        }
        rows[[length(rows) + 1L]] <- data.frame(
          scenario = sc,
          peel = peel,
          forecast_scenario = fsc,
          has_report = file.exists(file.path(run_dir, "Report.sso")),
          has_par = any(file.exists(file.path(retro_root, peel, c("ss3.par", "ss.par")))),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

count_retro_peels <- function(cfg, scenarios = NULL, retro_years = NULL) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  lr <- longrun_defaults()
  scenarios <- scenarios %||% cfg$retro_scenarios
  retro_years <- retro_years %||% lr$retro_years
  peel_tags <- retro_peel_tags(retro_years)
  rows <- list()
  for (sc in scenarios) {
    retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
    for (peel in peel_tags) {
      peel_dir <- file.path(retro_root, peel)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = sc,
        peel = peel,
        exists = dir.exists(peel_dir),
        has_par = any(file.exists(file.path(peel_dir, c("ss3.par", "ss.par")))),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

verify_pipeline <- function(cfg, verbose = TRUE) {
  source(file.path(cfg$appendix_dir, "R", "appendix-utils.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- longrun_defaults_from_env(longrun_defaults())
  n_peels <- length(lr$retro_years)
  n_fcast <- length(backtest_forecast_scenarios(lr))
  expected_hindcast <- length(cfg$retro_scenarios) * n_peels * n_fcast

  retro_df <- count_retro_peels(cfg)
  fcast_df <- count_forecast_reports(cfg)

  retro <- load_collected_retrospectives(cfg)
  catch <- load_collected_catch_backtests(cfg)

  ens_path <- file.path(cfg$appendix_dir, "data", "ensemble-weights", "ensemble_weight_comparison.csv")
  prof_path <- file.path(
    cfg$appendix_dir, "data", "likelihood-profiles", cfg$reference_scenario, "manifest.csv"
  )
  ens <- if (file.exists(ens_path)) {
    utils::read.csv(ens_path, stringsAsFactors = FALSE)
  } else if (diag_strict(cfg)) {
    stop("Ensemble weights not collected: ", ens_path, call. = FALSE)
  } else {
    NULL
  }
  prof <- if (file.exists(prof_path)) {
    utils::read.csv(prof_path, stringsAsFactors = FALSE)
  } else if (diag_strict(cfg)) {
    stop("Likelihood profile manifest not found: ", prof_path, call. = FALSE)
  } else {
    NULL
  }

  paper_figs <- cfg$paper_figures
  if (is.null(paper_figs) || !length(paper_figs)) {
    paper_figs <- c(
      "sma-fig-ts-1.png",
      "sma-fig-pe-phase-1.png",
      "selmat.png",
      "sma-fig-res-time-1.png",
      "sma-fig-res-acf-1.png",
      "4panel.png"
    )
  }
  fig_dir <- file.path(cfg$project_dir, "figs")
  fig_ok <- vapply(paper_figs, function(f) file.exists(file.path(fig_dir, f)), logical(1))

  content_files <- c(
    file.path(cfg$appendix_dir, "content", "internal", "checks.tex"),
    file.path(cfg$appendix_dir, "content", "external", "output.tex"),
    file.path(cfg$appendix_dir, "content", "robustness", "output.tex")
  )
  content_ok <- file.exists(content_files)

  retro_one_step_n <- function(retro) {
    if (is.null(retro)) {
      return(0L)
    }
    if (!is.null(retro$data$oneStep)) {
      return(NROW(retro$data$oneStep))
    }
    if (!is.null(retro$oneStep)) {
      return(NROW(retro$oneStep))
    }
    0L
  }

  checks <- list(
    retro_peels_with_par = sum(retro_df$has_par, na.rm = TRUE),
    retro_peels_expected = nrow(retro_df),
    hindcast_reports = sum(fcast_df$has_report, na.rm = TRUE),
    hindcast_expected = expected_hindcast,
    retro_one_step_rows = retro_one_step_n(retro),
    catch_forecast_rows = if (!is.null(catch) && !is.null(catch$summary)) {
      NROW(catch$summary)
    } else {
      0L
    },
    ensemble_weights = !is.null(ens) && NROW(ens) > 0L,
    profile_runs = if (!is.null(prof)) NROW(prof) else 0L,
    paper_figures = sum(fig_ok),
    paper_figures_expected = length(paper_figs),
    appendix_content = sum(content_ok),
    appendix_content_expected = length(content_files)
  )

  ok <- isTRUE(checks$retro_peels_with_par >= checks$retro_peels_expected) &&
    isTRUE(checks$hindcast_reports >= checks$hindcast_expected) &&
    isTRUE(checks$catch_forecast_rows >= checks$hindcast_expected) &&
    isTRUE(checks$retro_one_step_rows > 0L) &&
    isTRUE(checks$ensemble_weights) &&
    isTRUE(checks$profile_runs > 0L) &&
    isTRUE(checks$paper_figures == checks$paper_figures_expected) &&
    isTRUE(checks$appendix_content == checks$appendix_content_expected)

  out <- list(ok = ok, checks = checks, retro = retro_df, hindcast = fcast_df)
  if (isTRUE(verbose)) {
    message("=== Pipeline verification ===")
    message(
      "Retros: ",
      checks$retro_peels_with_par, "/", checks$retro_peels_expected,
      " peels with par; one-step rows=", checks$retro_one_step_rows
    )
    message(
      "Hindcasts: ",
      checks$hindcast_reports, "/", checks$hindcast_expected,
      " Report.sso; collected rows=", checks$catch_forecast_rows
    )
    message(
      "Ensemble weights: ", if (checks$ensemble_weights) "OK" else "MISSING",
      "; profile runs: ", checks$profile_runs
    )
    message(
      "Paper figs: ", checks$paper_figures, "/", checks$paper_figures_expected,
      "; appendix content files: ", checks$appendix_content, "/", checks$appendix_content_expected
    )
    message("Overall: ", if (ok) "PASS" else "FAIL (see checks)")
  }
  invisible(out)
}

write_pipeline_status <- function(appendix_dir, step, status, extra = list()) {
  root <- file.path(appendix_dir, "data", "pipeline")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  payload <- c(
    list(step = step, status = status, updated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    extra
  )
  jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
  if (jsonlite_ok) {
    jsonlite::write_json(payload, file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  }
  saveRDS(payload, file.path(root, "status.rds"))
  invisible(payload)
}
