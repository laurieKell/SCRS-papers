# Shared helpers for long-running appendix jobs (retros, catch back-tests, LOO CV).

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

parse_longrun_args <- function(args, defaults = list()) {
  out <- utils::modifyList(
    list(
      appendix_dir = NULL,
      study_id = NULL,
      scenario = NULL,
      scenarios = NULL,
      setup = FALSE,
      run = FALSE,
      collect = FALSE,
      all = FALSE,
      parallel = TRUE,
      workers = NULL,
      import_from = NULL,
      backtest_mode = NULL
    ),
    defaults
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
    } else if (key == "--scenario" && i < length(args)) {
      out$scenario <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--scenarios" && i < length(args)) {
      out$scenarios <- strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
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
      i <- i + 2L
    } else if (key == "--backtest-mode" && i < length(args)) {
      out$backtest_mode <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--no-parallel") {
      out$parallel <- FALSE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  if (isTRUE(out$all)) {
    out$setup <- isTRUE(out$setup) || identical(defaults$setup_on_all, TRUE)
    out$run <- TRUE
    out$collect <- TRUE
  }
  out
}

parse_int_csv <- function(x) {
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

longrun_cfg <- function(appendix_dir, study_id = NULL) {
  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  source(file.path(appendix_dir, "R", "study-config.R"), local = FALSE)
  source(file.path(appendix_dir, "config-case-study.R"), local = FALSE)
  source(file.path(appendix_dir, "config-longruns.R"), local = FALSE)
  case_study_defaults(appendix_dir = appendix_dir, study_id = study_id)
}

longrun_job_root <- function(cfg, job, scenario = NULL) {
  longrun_data_root(cfg$appendix_dir, job, scenario, data_dir = cfg$data_dir)
}

write_longrun_status <- function(root, status, extra = list()) {
  if (requireNamespace("diags", quietly = TRUE)) {
    return(diags::writeLongrunStatus(root, status, extra))
  }
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  payload <- c(list(status = status, updated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), extra)
  path_json <- file.path(root, "status.json")
  path_rds <- file.path(root, "status.rds")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(payload, path_json, auto_unbox = TRUE, pretty = TRUE)
  } else {
    saveRDS(payload, path_rds)
  }
  invisible(payload)
}

read_longrun_status <- function(root) {
  if (requireNamespace("diags", quietly = TRUE)) {
    st <- diags::readLongrunStatus(root)
    if (!is.null(st)) {
      return(st)
    }
  }
  path_json <- file.path(root, "status.json")
  path_rds <- file.path(root, "status.rds")
  if (file.exists(path_json)) {
    if (diag_strict()) {
      return(read_json_strict(path_json, "longrun status.json"))
    }
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      return(jsonlite::fromJSON(path_json, simplifyVector = TRUE))
    }
  }
  if (file.exists(path_rds)) {
    if (diag_strict()) {
      return(read_rds_strict(path_rds, "longrun status.rds"))
    }
    return(readRDS(path_rds))
  }
  NULL
}

longrun_log_path <- function(root) {
  file.path(root, "run.log")
}

longrun_defaults_from_env <- function(defaults) {
  if (nzchar(Sys.getenv("SS3_EXE", unset = ""))) {
    defaults$ss_exe <- Sys.getenv("SS3_EXE")
  }
  if (nzchar(Sys.getenv("DIAG_RETRO_YEARS", unset = ""))) {
    defaults$retro_years <- parse_int_csv(Sys.getenv("DIAG_RETRO_YEARS"))
  }
  if (nzchar(Sys.getenv("DIAG_RETRO_SCENARIOS", unset = ""))) {
    defaults$retro_scenarios <- strsplit(Sys.getenv("DIAG_RETRO_SCENARIOS"), ",", fixed = TRUE)[[1]]
  }
  if (nzchar(Sys.getenv("DIAG_CV_FLEETS", unset = ""))) {
    defaults$crossval_fleets <- parse_int_csv(Sys.getenv("DIAG_CV_FLEETS"))
  }
  if (nzchar(Sys.getenv("DIAG_CV_MC", unset = ""))) {
    defaults$crossval_mc <- as.integer(Sys.getenv("DIAG_CV_MC"))
  }
  if (nzchar(Sys.getenv("DIAG_SKIP_EXISTING_FORECASTS", unset = ""))) {
    val <- tolower(Sys.getenv("DIAG_SKIP_EXISTING_FORECASTS"))
    defaults$skip_existing_forecasts <- val %in% c("1", "true", "yes")
  }
  if (nzchar(Sys.getenv("DIAG_BACKTEST_MODE", unset = ""))) {
    defaults$backtest_mode <- Sys.getenv("DIAG_BACKTEST_MODE")
  }
  if (nzchar(Sys.getenv("DIAG_MANAGEMENT_YEARS", unset = ""))) {
    defaults$management_years <- as.integer(Sys.getenv("DIAG_MANAGEMENT_YEARS"))
  }
  defaults
}

backtest_mode <- function(lr) {
  mode <- lr$backtest_mode %||% "appendix"
  mode <- tolower(as.character(mode))
  if (!mode %in% c("appendix", "management")) {
    stop("backtest_mode must be 'appendix' or 'management', not '", mode, "'.", call. = FALSE)
  }
  mode
}

backtest_forecast_scenarios <- function(lr) {
  if (identical(backtest_mode(lr), "management")) {
    lr$management_forecast_scenarios %||% c("catchReported", "fZero", "fMsy")
  } else {
    lr$forecast_scenarios
  }
}

forecast_run_dir <- function(retro_root, peel_tag, forecast_scenario, lr) {
  if (identical(backtest_mode(lr), "management")) {
    file.path(
      retro_root,
      peel_tag,
      lr$management_out_subdir %||% "backtest",
      forecast_scenario
    )
  } else {
    file.path(retro_root, peel_tag, forecast_scenario)
  }
}

scenario_dirs <- function(cfg, scenarios) {
  paths <- file.path(cfg$ss_base, scenarios)
  ok <- dir.exists(paths)
  if (!all(ok)) {
    stop(
      "Missing scenario director(ies): ",
      paste(paths[!ok], collapse = ", "),
      call. = FALSE
    )
  }
  paths
}

parallel_lapply <- function(x, fun, parallel = TRUE, workers = NULL) {
  if (requireNamespace("diags", quietly = TRUE)) {
    return(diags::parallelLapply(x, fun, parallel = parallel, workers = workers))
  }
  if (isTRUE(parallel) && length(x) > 1L) {
    nw <- workers
    if (is.null(nw)) {
      nw <- max(1L, parallel::detectCores(logical = TRUE) - 2L)
    }
    nw <- min(as.integer(nw), length(x))
    if (nw > 1L && requireNamespace("future.apply", quietly = TRUE) &&
        requireNamespace("future", quietly = TRUE)) {
      old <- future::plan(future::multisession, workers = nw)
      on.exit(future::plan(old), add = TRUE)
      return(future.apply::future_lapply(x, fun))
    }
  }
  lapply(x, fun)
}

ss_input_files <- function() {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::ssInputFiles())
  }
  c("control.ss", "data.ss", "starter.ss", "forecast.ss")
}

copy_ss_inputs <- function(from_dir, to_dir) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::copySsInputs(from_dir, to_dir))
  }
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in ss_input_files()) {
    src <- file.path(from_dir, f)
    if (file.exists(src)) {
      file.copy(src, file.path(to_dir, f), overwrite = TRUE)
    }
  }
  invisible(to_dir)
}

run_ss3_dir <- function(dir, ss_exe, skip_finished = TRUE, show_in_console = FALSE) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    res <- sstruth::runSs3Dir(
      dir,
      ss_exe = ss_exe,
      skip_finished = skip_finished,
      show_in_console = show_in_console,
      save_output = TRUE
    )
    if (isTRUE(diag_strict()) && !isTRUE(res$ok) && !isTRUE(res$skipped)) {
      stop("SS3 run failed in ", dir, call. = FALSE)
    }
    return(res)
  }
  if (isTRUE(skip_finished) && file.exists(file.path(dir, "Report.sso"))) {
    return(list(ok = TRUE, skipped = TRUE, dir = dir))
  }
  ok <- tryCatch({
    r4ss::run(
      dir,
      exe = ss_exe,
      skipfinished = isTRUE(skip_finished),
      show_in_console = isTRUE(show_in_console)
    )
    TRUE
  }, error = function(e) {
    if (diag_strict()) {
      stop("SS3 run failed in ", dir, ": ", conditionMessage(e), call. = FALSE)
    }
    message("SS3 run failed in ", dir, ": ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) {
    ss_rep <- if (diag_strict()) {
      r4ss::SS_output(
        dir,
        verbose = FALSE,
        printstats = FALSE,
        covar = FALSE,
        hidewarn = TRUE,
        NoCompOK = TRUE
      )
    } else {
      tryCatch(
        r4ss::SS_output(
          dir,
          verbose = FALSE,
          printstats = FALSE,
          covar = FALSE,
          hidewarn = TRUE,
          NoCompOK = TRUE
        ),
        error = function(e) NULL
      )
    }
    if (!is.null(ss_rep)) {
      save(ss_rep, file = file.path(dir, "ss.RData"))
    } else if (diag_strict()) {
      stop("SS3 output read failed after run in ", dir, call. = FALSE)
    }
  } else if (diag_strict()) {
    stop("SS3 run failed in ", dir, call. = FALSE)
  }
  list(ok = isTRUE(ok), skipped = FALSE, dir = dir)
}

cpue_fleets_from_data <- function(data_path) {
  if (!file.exists(data_path)) {
    return(NULL)
  }
  lines <- readLines(data_path, warn = FALSE)
  hdr <- grep("^#_type[[:space:]]+surveytiming", lines)
  if (!length(hdr)) {
    return(NULL)
  }
  start <- hdr[[1]] + 1L
  end <- grep("^#_CPUE_and_surveyabundance_observations", lines)[[1]] - 1L
  if (is.na(end) || end < start) {
    return(NULL)
  }
  block <- lines[seq(start, end)]
  block <- block[nzchar(trimws(block)) & !grepl("^#", trimws(block))]
  rows <- lapply(block, function(ln) {
    comment <- sub("^.*(#_.*)$", "\\1", ln)
    if (!grepl("^#_", comment)) {
      return(NULL)
    }
    fleet_id <- suppressWarnings(as.integer(sub("^#_", "", comment)))
    name_part <- sub("[[:space:]]#_.*$", "", trimws(ln))
    parts <- strsplit(trimws(name_part), "[[:space:]]+", perl = TRUE)[[1]]
    if (length(parts) < 1L || parts[[1]] != "3") {
      return(NULL)
    }
    fleet_name <- parts[[length(parts)]]
    data.frame(
      fleet_id = fleet_id,
      fleet_name = fleet_name,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$fleet_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

patch_data_ss_loo_cpue <- function(data_path, fleet_id) {
  lines <- readLines(data_path, warn = FALSE)
  cpue_hdr <- grep("^#_CPUE_data", lines)
  if (!length(cpue_hdr)) {
    stop("No #_CPUE_data section in ", data_path, call. = FALSE)
  }
  start <- cpue_hdr[[1]] + 2L
  end_candidates <- grep("^#_", lines)
  end_candidates <- end_candidates[end_candidates > start]
  end <- if (length(end_candidates)) end_candidates[[1]] - 1L else length(lines)
  keep_lines <- character(0)
  for (i in seq(start, end)) {
    ln <- lines[[i]]
    if (!nzchar(trimws(ln)) || grepl("^#", trimws(ln))) {
      keep_lines <- c(keep_lines, ln)
      next
    }
    parts <- strsplit(trimws(ln), "[[:space:]]+", perl = TRUE)[[1]]
    if (length(parts) < 3L) {
      keep_lines <- c(keep_lines, ln)
      next
    }
    idx <- suppressWarnings(as.integer(parts[[3]]))
    if (is.finite(idx) && idx == as.integer(fleet_id)) {
      next
    }
    keep_lines <- c(keep_lines, ln)
  }
  lines <- c(lines[seq_len(start - 1L)], keep_lines, lines[seq(end + 1L, length(lines))])
  writeLines(lines, data_path, useBytes = FALSE)
  invisible(data_path)
}

terminal_bratio_f <- function(rep, year = NULL) {
  if (is.null(rep)) {
    return(list(bratio = NA_real_, f = NA_real_, fmsy = NA_real_, ffmsy = NA_real_))
  }
  dq <- rep[["derived_quants"]]
  if (is.null(year)) {
    year <- rep$endyr
  }
  br <- NA_real_
  f <- NA_real_
  fmsy <- NA_real_
  if (!is.null(dq) && NROW(dq)) {
    br_lab <- paste0("Bratio_", year)
    f_lab <- paste0("F_", year)
    hit <- dq[["Label"]] == br_lab
    if (any(hit)) {
      br <- suppressWarnings(as.numeric(dq[hit, "Value"][[1]]))
    }
    hit <- dq[["Label"]] == f_lab
    if (any(hit)) {
      f <- suppressWarnings(as.numeric(dq[hit, "Value"][[1]]))
    }
    hit <- dq[["Label"]] == "annF_MSY"
    if (any(hit)) {
      fmsy <- suppressWarnings(as.numeric(dq[hit, "Value"][[1]]))
    }
  }
  ffmsy <- if (is.finite(f) && is.finite(fmsy) && fmsy > 0) f / fmsy else NA_real_
  list(bratio = br, f = f, fmsy = fmsy, ffmsy = ffmsy, year = year)
}

forecast_summary_from_rep <- function(rep, scenario, peel, forecast_scenario, run_dir) {
  if (is.null(rep)) {
    return(NULL)
  }
  ts <- rep$timeseries
  if (is.null(ts) || !NROW(ts)) {
    return(NULL)
  }
  fore <- ts[ts$Era == "FORE", , drop = FALSE]
  if (!NROW(fore)) {
    return(NULL)
  }
  yr <- min(fore$Yr, na.rm = TRUE)
  ix <- which(fore$Yr == yr)[1]
  spawn <- if ("SpawnBio" %in% names(fore)) fore$SpawnBio[[ix]] else NA_real_
  catch <- if ("Dead_Catch" %in% names(fore)) fore$Dead_Catch[[ix]] else NA_real_
  fval <- if ("F" %in% names(fore)) fore$F[[ix]] else NA_real_
  data.frame(
    scenario = scenario,
    peel = peel,
    forecast_scenario = forecast_scenario,
    forecast_year = yr,
    spawn_bio = spawn,
    dead_catch = catch,
    sum_f = fval,
    run_dir = run_dir,
    ok = TRUE,
    stringsAsFactors = FALSE
  )
}

launch_longrun_job <- function(script, appendix_dir, args, root, label) {
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  log_path <- longrun_log_path(root)
  cmd_args <- c("--vanilla", script, "--appendix-dir", appendix_dir, args)
  write_longrun_status(root, "queued", list(job = label, log = log_path))
  system2(
    normalizePath(Sys.which("Rscript"), winslash = "/", mustWork = TRUE),
    cmd_args,
    stdout = log_path,
    stderr = log_path,
    wait = FALSE
  )
  message(label, " job queued in background; log: ", log_path)
  invisible(log_path)
}
