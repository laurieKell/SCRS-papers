# SS3 assessment backend — all r4ss / sstruth I/O and collection for 3diags.
# Templates and collectors should call backend_*() via backend.R, not packages directly.

ss3_require_r4ss <- function() {
  if (!requireNamespace("r4ss", quietly = TRUE)) {
    stop("Package 'r4ss' is required for SS3 diagnostics.", call. = FALSE)
  }
  invisible(TRUE)
}

ss3_require_sstruth <- function() {
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required for SS3 diagnostics.", call. = FALSE)
  }
  invisible(TRUE)
}

backend_capabilities <- function(cfg) {
  caps <- cfg$capabilities %||% list()
  list(
    read_output = TRUE,
    retrospectives = isTRUE(caps$retrospectives %||% TRUE),
    catch_backtests = isTRUE(caps$catch_backtests %||% TRUE),
    likelihood_profiles = isTRUE(caps$likelihood_profiles %||% FALSE),
    mvln_hessian = isTRUE(caps$mvln_hessian %||% TRUE),
    loo_crossvalidation = isTRUE(caps$loo_crossvalidation %||% FALSE),
    length_composition = isTRUE(caps$length_composition %||% TRUE),
    ss3diags_runs_test = isTRUE(caps$ss3diags_runs_test %||% TRUE)
  )
}

backend_read_output <- function(path, covar = FALSE, forecast = FALSE, ...) {
  ss3_require_r4ss()
  r4ss::SS_output(
    dir = path,
    verbose = FALSE,
    printstats = FALSE,
    forecast = isTRUE(forecast),
    covar = isTRUE(covar),
    hidewarn = TRUE,
    NoCompOK = TRUE
  )
}

ssRunDir <- function(path) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::ssRunDir(path))
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(path)) {
    return(path)
  }
  if (file.exists(file.path(path, "Report.sso"))) {
    return(path)
  }
  inner <- file.path(path, basename(path))
  if (dir.exists(inner) && file.exists(file.path(inner, "Report.sso"))) {
    return(normalizePath(inner, winslash = "/", mustWork = TRUE))
  }
  path
}

backend_scenario_paths <- function(cfg, scenarios = NULL) {
  scenarios <- scenarios %||% cfg$scenarios
  raw <- file.path(cfg$ss_base, scenarios)
  vapply(raw, ssRunDir, character(1), USE.NAMES = FALSE)
}

backend_reference_path <- function(cfg) {
  ssRunDir(file.path(cfg$ss_base, cfg$reference_scenario))
}

backend_backtest_params <- function(cfg) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  longrun_params_from_cfg(cfg)
}

backend_retro_peel_tags <- function(retro_years) {
  ss3_require_sstruth()
  sstruth::retroTag(as.integer(retro_years))
}

backend_retro_dirs_available <- function(cfg) {
  bt <- backend_backtest_params(cfg)
  retro_paths <- file.path(bt$ssBase, bt$scenarioDirs, bt$retroSubdir)
  any(dir.exists(retro_paths))
}

backend_run_retros <- function(cfg, dirs, bt = NULL, ...) {
  ss3_require_sstruth()
  bt <- bt %||% backend_backtest_params(cfg)
  bt$runRetro <- FALSE
  dots <- list(...)
  for (nm in names(dots)) {
    bt[[nm]] <- dots[[nm]]
  }
  sstruth::runRetros(dirs, bt)
}

backend_run_forecast_peels <- function(cfg, retro_root, lr = NULL, workers = NULL, ...) {
  ss3_require_sstruth()
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  sstruth::runFlevelProjectionsPeels(
    retroRoot = retro_root,
    ssExe = lr$ss_exe,
    ssVersion = lr$ss_version,
    scenarios = lr$forecast_scenarios,
    retroYears = lr$retro_years,
    runLogical = TRUE,
    skipExistingForecastRuns = isTRUE(lr$skip_existing_forecasts),
    nohess = isTRUE(lr$forecast_nohess),
    parallelPeels = isTRUE(lr$parallel_forecast_peels),
    peelWorkers = workers %||% lr$forecast_peel_workers,
    ...
  )
}

backend_run_management_backtests <- function(cfg, retro_root, lr = NULL, workers = NULL, ...) {
  ss3_require_sstruth()
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  sstruth::runBacktestProjectionsPeels(
    retroRoot = retro_root,
    ssExe = lr$ss_exe,
    ssVersion = lr$ss_version,
    retroYears = lr$retro_years,
    outSubdir = lr$management_out_subdir %||% "backtest",
    managementYears = lr$management_years %||% 100L,
    runLogical = TRUE,
    skipExistingForecastRuns = isTRUE(lr$skip_existing_forecasts),
    nohess = isTRUE(lr$forecast_nohess),
    parallelPeels = isTRUE(lr$parallel_forecast_peels),
    peelWorkers = workers %||% lr$forecast_peel_workers,
    ...
  )
}

backend_run_backtests <- function(cfg, retro_root, lr = NULL, workers = NULL, ...) {
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  if (identical(backtest_mode(lr), "management")) {
    backend_run_management_backtests(cfg, retro_root, lr = lr, workers = workers, ...)
  } else {
    backend_run_forecast_peels(cfg, retro_root, lr = lr, workers = workers, ...)
  }
}

backend_retro_backtest_data <- function(cfg, bt = NULL) {
  ss3_require_sstruth()
  bt <- bt %||% backend_backtest_params(cfg)
  bt$runRetro <- FALSE
  sstruth::retroBacktestData(bt)
}

backend_collect_retrospectives <- function(cfg, scenarios = NULL, import_from = NULL) {
  ss3_require_sstruth()
  scenarios <- scenarios %||% cfg$retro_scenarios

  if (!is.null(import_from)) {
    import_from <- normalizePath(import_from, winslash = "/", mustWork = TRUE)
    src_rds <- file.path(import_from, "collected.rds")
    if (file.exists(src_rds)) {
      return(readRDS(src_rds))
    }
    if (dir.exists(import_from)) {
      bt <- backend_backtest_params(cfg)
      bt$ssBase <- import_from
      bt$scenarioDirs <- scenarios
      bt$runRetro <- FALSE
      return(sstruth::retroBacktestData(bt))
    }
    stop("import_from must be collected.rds parent or ss_base directory.", call. = FALSE)
  }

  bt <- backend_backtest_params(cfg)
  bt$scenarioDirs <- scenarios
  bt$runRetro <- FALSE
  backend_retro_backtest_data(cfg, bt)
}

backend_load_internal <- function(cfg, path) {
  ss3_require_sstruth()
  if (!dir.exists(path)) {
    stop("SS3 scenario path not found: ", path, call. = FALSE)
  }
  list(
    diag = sstruth::diagSS(path),
    ssu = sstruth::ssU(path),
    ss_len = sstruth::ssLen(path),
    run = sstruth::ssRunSummary(path, verbose = FALSE)
  )
}

backend_load_external <- function(cfg, path) {
  ss3_require_sstruth()
  if (!dir.exists(path)) {
    stop("SS3 scenario path not found: ", path, call. = FALSE)
  }
  rep <- backend_read_output(path)
  list(
    diag = sstruth::diagSS(path),
    sp = sstruth::spSsbFromSs3(rep),
    ss_pe = sstruth::ssPeDiagnose(path)
  )
}

backend_load_robustness <- function(cfg, retro = NULL) {
  ss3_require_sstruth()
  if (is.null(retro)) {
    if (backend_retro_dirs_available(cfg)) {
      retro <- backend_retro_backtest_data(cfg)
    } else {
      stop(
        "Retrospective directories not found under ", cfg$ss_base,
        "; run appendices/R/run-retrospectives.R --run then --collect.",
        call. = FALSE
      )
    }
  }
  paths <- backend_scenario_paths(cfg)
  paths <- paths[dir.exists(paths)]
  if (!length(paths)) {
    stop("No SS3 scenario directories found for robustness diagnostics.", call. = FALSE)
  }
  rows <- lapply(paths, function(p) {
    d <- sstruth::ssPeDiagnose(p)
    cbind(
      scenario = gsub("Scenario-", "Scenario ", basename(p)),
      as.data.frame(d, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  })
  list(retro = retro, ss_pe = do.call(rbind, rows))
}

backend_load_pe <- function(cfg, paths = NULL, for_plots = FALSE, parallel = TRUE, workers = NULL) {
  ss3_require_sstruth()
  paths <- paths %||% backend_scenario_paths(cfg)
  paths <- paths[dir.exists(paths)]
  if (!length(paths)) {
    return(NULL)
  }
  use_parallel <- isTRUE(parallel)
  if (isTRUE(for_plots) && identical(Sys.info()[["sysname"]], "Windows")) {
    use_parallel <- isTRUE(getOption("diag.parallel_curvess", FALSE))
  }
  runs <- data.frame(
    id = basename(paths),
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  pe <- sstruth::ssCurve(
    runs,
    col = "id",
    for_plots = isTRUE(for_plots),
    cache = TRUE,
    parallel = use_parallel,
    workers = workers
  )
  if (!is.null(pe) && "run" %in% names(pe$tseries) && !"id" %in% names(pe$tseries)) {
    for (nm in names(pe)) {
      if (is.data.frame(pe[[nm]]) && "run" %in% names(pe[[nm]])) {
        names(pe[[nm]])[names(pe[[nm]]) == "run"] <- "id"
      }
    }
  }
  pe
}

backend_load_fls <- function(cfg, paths = NULL) {
  if (!requireNamespace("FLCore", quietly = TRUE)) {
    stop("Package 'FLCore' is required to read SS3 scenarios.", call. = FALSE)
  }
  paths <- paths %||% backend_scenario_paths(cfg)
  paths <- paths[dir.exists(paths)]
  if (!length(paths)) {
    return(NULL)
  }
  read_fn <- if (requireNamespace("sstruth", quietly = TRUE)) {
    sstruth::readTruth
  } else if (requireNamespace("ss3om", quietly = TRUE)) {
    ss3om::readFLSss3
  } else {
    stop("Package 'sstruth' or 'ss3om' is required to read SS3 scenarios.", call. = FALSE)
  }
  stocks <- lapply(paths, read_fn)
  names(stocks) <- basename(paths)
  FLCore::FLStocks(stocks)
}

backend_load_case_study <- function(cfg) {
  paths <- backend_scenario_paths(cfg)
  paths_ok <- paths[dir.exists(paths)]
  list(
    paths = paths_ok,
    pe = backend_load_pe(cfg, paths_ok, for_plots = TRUE),
    fls = backend_load_fls(cfg, paths_ok)
  )
}

backend_read_forecast_output <- function(path, strict = TRUE) {
  ss3_require_r4ss()
  if (isTRUE(strict)) {
    backend_read_output(path, forecast = TRUE)
  } else {
    tryCatch(
      backend_read_output(path, forecast = TRUE),
      error = function(e) NULL
    )
  }
}

backend_read_profile_output <- function(run_dir) {
  ss3_require_sstruth()
  sstruth::readSs3Run(run_dir)
}

backend_diag_ss <- function(path) {
  ss3_require_sstruth()
  sstruth::diagSS(path)
}

backend_plot_cpue_residuals <- function(diag_df) {
  ss3_require_sstruth()
  sstruth::plotResiduals(diag_df)
}

backend_plot_cpue_std_residuals <- function(diag_df) {
  ss3_require_sstruth()
  sstruth::plotStdResiduals(diag_df)
}

backend_plot_sp_ssb <- function(sp_df) {
  ss3_require_sstruth()
  sstruth::plotSpSsb(sp_df, addSmooth = TRUE, showYears = FALSE)
}

backend_plot_retro_ratios <- function(retro_data, fore_col, ref_col, ylab, title) {
  ss3_require_sstruth()
  sstruth::plotRatios(
    retro_data$ratioHorizon,
    foreCol = fore_col,
    refCol = ref_col,
    ylab = ylab,
    title = title
  )
}

backend_summarise_run <- function(path, rep = NULL) {
  ss3_require_r4ss()
  rep <- rep %||% backend_read_output(path)
  info <- r4ss::SSsummarize(list(rep), verbose = FALSE)
  total_ll <- info$likelihoods[
    info$likelihoods$Label == "TOTAL",
    "model1",
    drop = TRUE
  ]
  data.frame(
    scenario = basename(path),
    n_params = info$npars,
    total_ll = round(total_ll, 3),
    max_grad = signif(info$maxgrad, 3),
    log_det_hessian = if (!is.null(rep$log_det_hessian)) {
      signif(rep$log_det_hessian, 4)
    } else {
      NA_real_
    },
    n_iterations = rep$N_iterations %||% NA_integer_,
    n_warnings = rep$Nwarnings %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

backend_collect_catch_backtests <- function(cfg, scenarios = NULL, lr = NULL) {
  source(file.path(cfg$appendix_dir, "config-longruns.R"), local = TRUE)
  source(file.path(cfg$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  lr <- lr %||% longrun_defaults_from_env(longrun_defaults())
  mode <- backtest_mode(lr)
  forecast_scenarios <- backtest_forecast_scenarios(lr)
  scenarios <- scenarios %||% cfg$retro_scenarios

  rows <- list()
  manifest <- list()
  for (sc in scenarios) {
    retro_root <- file.path(cfg$ss_base, sc, lr$retro_subdir)
    if (!dir.exists(retro_root)) {
      next
    }
    peel_tags <- backend_retro_peel_tags(lr$retro_years)
    peel_tags <- peel_tags[dir.exists(file.path(retro_root, peel_tags))]
    for (peel in peel_tags) {
      peel_yr <- sub("^retro", "", peel)
      for (fsc in forecast_scenarios) {
        run_dir <- forecast_run_dir(retro_root, peel, fsc, lr)
        if (!dir.exists(run_dir)) {
          next
        }
        rep <- backend_read_forecast_output(run_dir, strict = isTRUE(cfg$strict_errors %||% TRUE))
        row <- forecast_summary_from_rep(rep, sc, peel_yr, fsc, run_dir)
        if (!is.null(row)) {
          row$backtest_mode <- mode
          rows[[length(rows) + 1L]] <- row
          manifest[[length(manifest) + 1L]] <- data.frame(
            scenario = sc,
            peel = peel_yr,
            forecast_scenario = fsc,
            backtest_mode = mode,
            dir = run_dir,
            ok = !is.null(rep),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(rows)) {
    if (diag_strict(cfg)) {
      stop("Catch back-test collection produced no rows.", call. = FALSE)
    }
    return(NULL)
  }
  list(
    summary = do.call(rbind, rows),
    manifest = do.call(rbind, manifest)
  )
}
