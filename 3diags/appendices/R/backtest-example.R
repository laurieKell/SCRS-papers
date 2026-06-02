# Reusable management back-test workflow (sstruth::ssBacktest).
# Source from run-backtest-example.R or ss-backtest-example.Rmd.

parse_peel_years <- function(x) {
  if (is.numeric(x)) {
    return(sort(unique(as.integer(x))))
  }
  as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}

#' Resolve one SS3 scenario directory from study spec + relative id
resolve_scenario_dir <- function(ss_base, scenario_id) {
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required.", call. = FALSE)
  }
  path <- file.path(ss_base, scenario_id)
  sstruth::ssRunDir(path)
}

#' Default back-test configuration (override fields as needed)
backtest_example_config <- function(
  appendix_dir,
  study_id = Sys.getenv("DIAG_STUDY", unset = "sma-2026"),
  scenario_id = Sys.getenv("DIAG_BACKTEST_SCENARIO", unset = ""),
  ss_base = Sys.getenv("DIAG_SS_BASE", unset = ""),
  peel_years = parse_peel_years(Sys.getenv("DIAG_BACKTEST_PEEL_YEARS", unset = "0,1,5,10,20")),
  management_years = as.integer(Sys.getenv("DIAG_MANAGEMENT_YEARS", unset = "100")),
  retro_subdir = "retrospectives",
  backtest_subdir = "backtest",
  forecast_horizons = 1:3,
  ss_exe = Sys.getenv("SS3_EXE", unset = "ss3"),
  ss_version = "3.30"
) {
  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  spec_path <- file.path(appendix_dir, "studies", paste0(study_id, ".R"))
  if (!file.exists(spec_path)) {
    stop("Study spec not found: ", spec_path, call. = FALSE)
  }
  source(spec_path, local = TRUE)
  spec <- study_spec()

  if (!nzchar(ss_base)) {
    ss_base <- spec$paths$assessment_base
  }
  if (!nzchar(scenario_id)) {
    scenario_id <- spec$scenarios$retrospective[1]
    if (is.null(scenario_id) || !length(scenario_id)) {
      scenario_id <- spec$scenarios$ids[1]
    }
  }

  run_dir <- resolve_scenario_dir(ss_base, scenario_id)
  label <- spec$scenarios$register$label[
    match(scenario_id, spec$scenarios$register$id)
  ]
  if (is.na(label)) {
    label <- scenario_id
  }

  list(
    appendix_dir = appendix_dir,
    study_id = study_id,
    scenario_id = scenario_id,
    scenario_label = label,
    ss_base = normalizePath(ss_base, winslash = "/", mustWork = FALSE),
    run_dir = run_dir,
    peel_years = peel_years,
    management_years = management_years,
    retro_subdir = retro_subdir,
    backtest_subdir = backtest_subdir,
    forecast_horizons = forecast_horizons,
    ss_exe = ss_exe,
    ss_version = ss_version,
    forecast_scenarios = c("catchReported", "fZero", "fMsy")
  )
}

backtest_output_root <- function(cfg, job = "management-backtest") {
  file.path(cfg$appendix_dir, "data", cfg$study_id, job, cfg$scenario_id)
}

ss_backtest_params <- function(cfg) {
  sstruth::ssBacktestParams(
    runDir = cfg$run_dir,
    peelYears = cfg$peel_years,
    retroSubdir = cfg$retro_subdir,
    forecastHorizons = cfg$forecast_horizons,
    runRetro = FALSE,
    ssExe = cfg$ss_exe,
    ssVersion = cfg$ss_version
  )
}

#' Run retrospective re-fits (expensive; needs SS3)
run_backtest_retros <- function(cfg, workers = NULL) {
  params <- ss_backtest_params(cfg)
  params$runRetro <- TRUE
  params$retroWorkers <- workers
  message("[backtest] retros: ", cfg$scenario_label, " (", cfg$run_dir, ")")
  message("[backtest] peel years: ", paste(cfg$peel_years, collapse = ", "))
  sstruth::runRetros(cfg$run_dir, params)
  invisible(params)
}

#' Run two-stage management projections (reported catch -> F=0 / F=FMSY)
run_backtest_forecasts <- function(cfg, workers = NULL, parallel_peels = FALSE) {
  retro_root <- file.path(cfg$run_dir, cfg$retro_subdir)
  if (!dir.exists(retro_root)) {
    stop(
      "Retrospectives not found: ", retro_root,
      "\nRun run_backtest_retros() first.",
      call. = FALSE
    )
  }
  message("[backtest] forecasts: ", cfg$scenario_label)
  message(
    "[backtest] hindcast (reported catch) then F=0 / F=FMSY for ",
    cfg$management_years, " yr"
  )
  sstruth::runBacktestProjectionsPeels(
    retroRoot = retro_root,
    retroYears = sstruth::retroYearsFromPeel(cfg$peel_years),
    outSubdir = cfg$backtest_subdir,
    managementYears = cfg$management_years,
    ssExe = cfg$ss_exe,
    ssVersion = cfg$ss_version,
    parallelPeels = isTRUE(parallel_peels),
    peelWorkers = workers,
    skipExistingForecastRuns = TRUE,
    nohess = TRUE
  )
}

#' Collect retrospective one-step / ratio tables
collect_backtest_retro <- function(cfg) {
  params <- ss_backtest_params(cfg)
  sstruth::retroBacktestData(params)
}

#' Resolve forecast run directory (management backtest/ or legacy peel/scenario/)
resolve_forecast_run_dir <- function(retro_root, peel, forecast_scenario, backtest_subdir = "backtest") {
  candidates <- unique(c(
    file.path(retro_root, peel, backtest_subdir, forecast_scenario),
    file.path(retro_root, peel, forecast_scenario)
  ))
  for (run_dir in candidates) {
    if (file.exists(file.path(run_dir, "Report.sso"))) {
      return(run_dir)
    }
  }
  NA_character_
}

#' Collect forecast trajectories from retro peel forecast folders
collect_management_forecasts <- function(cfg) {
  if (!requireNamespace("r4ss", quietly = TRUE)) {
    stop("Package 'r4ss' is required.", call. = FALSE)
  }
  retro_root <- file.path(cfg$run_dir, cfg$retro_subdir)
  if (!dir.exists(retro_root)) {
    return(NULL)
  }
  peel_tags <- sstruth::retroTag(sstruth::retroYearsFromPeel(cfg$peel_years))
  peel_tags <- peel_tags[dir.exists(file.path(retro_root, peel_tags))]

  rows <- list()
  manifest <- list()
  for (peel in peel_tags) {
    peel_yr <- sub("^retro", "", peel)
    for (fsc in cfg$forecast_scenarios) {
      run_dir <- resolve_forecast_run_dir(
        retro_root,
        peel,
        fsc,
        backtest_subdir = cfg$backtest_subdir
      )
      if (is.na(run_dir)) {
        next
      }
      rep <- tryCatch(
        r4ss::SS_output(run_dir, forecast = TRUE, verbose = FALSE, printstats = FALSE),
        error = function(e) NULL
      )
      if (is.null(rep) || is.null(rep$timeseries)) {
        next
      }
      ts <- rep$timeseries
      era_col <- if ("Era" %in% names(ts)) "Era" else NULL
      if (!is.null(era_col)) {
        ts <- ts[as.character(ts[[era_col]]) %in% c("TIME", "FORE"), , drop = FALSE]
      }
      if (!NROW(ts)) {
        next
      }
      spawn <- if ("SpawnBio" %in% names(ts)) ts$SpawnBio else ts$SSB
      catch <- if ("Dead_Catch" %in% names(ts)) ts$Dead_Catch else NA_real_
      fval <- if ("F" %in% names(ts)) ts$F else NA_real_
      df <- data.frame(
        scenario_id = cfg$scenario_id,
        scenario_label = cfg$scenario_label,
        peel = peel_yr,
        forecast_scenario = fsc,
        year = ts$Yr,
        era = if (!is.null(era_col)) as.character(ts[[era_col]]) else NA_character_,
        spawn_bio = spawn,
        dead_catch = catch,
        sum_f = fval,
        run_dir = run_dir,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- df
      manifest[[length(manifest) + 1L]] <- data.frame(
        peel = peel_yr,
        forecast_scenario = fsc,
        dir = run_dir,
        n_years = NROW(df),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(NULL)
  }
  list(
    trajectories = do.call(rbind, rows),
    manifest = do.call(rbind, manifest)
  )
}

#' Write collected tables under appendices/data/<study>/management-backtest/
collect_backtest_example <- function(cfg) {
  root <- backtest_output_root(cfg)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  retro <- tryCatch(collect_backtest_retro(cfg), error = function(e) {
    message("[backtest] retro collect: ", conditionMessage(e))
    NULL
  })

  fcast <- collect_management_forecasts(cfg)

  out <- list(root = root, config = cfg, retro = retro, forecasts = fcast)
  if (!is.null(retro)) {
    if (!is.null(retro$oneStep) && NROW(retro$oneStep)) {
      utils::write.csv(retro$oneStep, file.path(root, "one-step.csv"), row.names = FALSE)
    }
    if (!is.null(retro$ratioHorizon) && NROW(retro$ratioHorizon)) {
      utils::write.csv(retro$ratioHorizon, file.path(root, "ratio-horizon.csv"), row.names = FALSE)
    }
    saveRDS(retro, file.path(root, "retro-collected.rds"))
  }
  if (!is.null(fcast) && NROW(fcast$trajectories)) {
    utils::write.csv(
      fcast$trajectories,
      file.path(root, "forecast-trajectories.csv"),
      row.names = FALSE
    )
    utils::write.csv(fcast$manifest, file.path(root, "forecast-manifest.csv"), row.names = FALSE)
  }
  saveRDS(out, file.path(root, "collected.rds"))
  jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
  status <- list(
    scenario = cfg$scenario_id,
    peel_years = cfg$peel_years,
    has_retro = !is.null(retro),
    has_forecasts = !is.null(fcast) && NROW(fcast$trajectories) > 0,
    updated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  if (jsonlite_ok) {
    jsonlite::write_json(status, file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  }
  out
}

#' Full workflow wrapper
run_backtest_example <- function(
  cfg,
  run_retro = FALSE,
  run_forecast = FALSE,
  collect = TRUE,
  workers = NULL,
  parallel_peels = FALSE
) {
  if (isTRUE(run_retro)) {
    run_backtest_retros(cfg, workers = workers)
  }
  forecast_status <- NULL
  if (isTRUE(run_forecast)) {
    forecast_status <- run_backtest_forecasts(
      cfg,
      workers = workers,
      parallel_peels = parallel_peels
    )
  }
  out <- NULL
  if (isTRUE(collect)) {
    out <- collect_backtest_example(cfg)
    out$forecast_status <- forecast_status
  } else {
    out <- list(config = cfg, forecast_status = forecast_status)
  }
  invisible(out)
}
