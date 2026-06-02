#!/usr/bin/env Rscript
#
# CLI wrapper for sstruth SS_output collect / bind workflow.
# Study-specific paths come from appendices/studies/<study_id>.R.
#
# Phase 1 — cache SS_output in each run folder:
#   Rscript appendices/R/extract-ss-output.R \
#     --study sma-2026 --include-sensitivity --collect --workers 4
#
# Phase 2 — bind Kobe from cache:
#   Rscript appendices/R/extract-ss-output.R \
#     --study sma-2026 --include-sensitivity --from-cache \
#     --slot Kobe --out appendices/data/sma-2026/kobe-all.csv
#
# Environment: DIAG_SS_BASE overrides study assessment_base path.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

requireSstruth <- function() {
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop(
      "Package 'sstruth' is required. ",
      "Install/load dev version from C:/active/flr/sstruth.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

parseExtractArgs <- function(args) {
  out <- list(
    appendixDir = NULL,
    studyId = Sys.getenv("DIAG_STUDY", unset = "sma-2026"),
    ssBase = NULL,
    scenarios = NULL,
    includeSensitivity = FALSE,
    discover = NULL,
    slot = NULL,
    out = NULL,
    saveReps = NULL,
    loadReps = NULL,
    covar = FALSE,
    forecast = FALSE,
    listSlots = FALSE,
    scan = NULL,
    collect = FALSE,
    fromCache = FALSE,
    refresh = FALSE,
    parallel = TRUE,
    workers = NULL,
    scenarioCol = "scenario"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--appendix-dir", "--study", "--ss-base", "--scenarios",
                   "--slot", "--element", "--out", "--save-reps", "--load-reps",
                   "--scenario-col", "--scan", "--workers") &&
      i < length(args)) {
      val <- args[[i + 1L]]
      if (key == "--appendix-dir") out$appendixDir <- val
      else if (key == "--study") out$studyId <- val
      else if (key == "--ss-base") out$ssBase <- val
      else if (key == "--scenarios") out$scenarios <- val
      else if (key %in% c("--slot", "--element")) out$slot <- val
      else if (key == "--out") out$out <- val
      else if (key == "--save-reps") out$saveReps <- val
      else if (key == "--load-reps") out$loadReps <- val
      else if (key == "--scenario-col") out$scenarioCol <- val
      else if (key == "--scan") out$scan <- val
      else if (key == "--workers") out$workers <- val
      i <- i + 2L
    } else if (key == "--include-sensitivity") {
      out$includeSensitivity <- TRUE
      i <- i + 1L
    } else if (key == "--discover") {
      out$discover <- TRUE
      i <- i + 1L
    } else if (key == "--no-discover") {
      out$discover <- FALSE
      i <- i + 1L
    } else if (key %in% c("--covar", "--forecast", "--list-slots", "--list-elements",
                          "--collect", "--from-cache", "--refresh", "--parallel")) {
      if (key %in% c("--list-slots", "--list-elements")) out$listSlots <- TRUE
      else if (key == "--from-cache") out$fromCache <- TRUE
      else if (key == "--collect") out$collect <- TRUE
      else if (key == "--refresh") out$refresh <- TRUE
      else if (key == "--parallel") out$parallel <- TRUE
      else if (key == "--covar") out$covar <- TRUE
      else if (key == "--forecast") out$forecast <- TRUE
      i <- i + 1L
    } else if (key == "--no-parallel") {
      out$parallel <- FALSE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  if (!is.null(out$scenarios) && nzchar(out$scenarios)) {
    out$scenarios <- trimws(strsplit(out$scenarios, ",", fixed = TRUE)[[1]])
  }
  if (!is.null(out$workers) && nzchar(out$workers)) {
    out$workers <- as.integer(out$workers)
  }
  out
}

extractHelp <- function() {
  cat(
    "Usage: Rscript appendices/R/extract-ss-output.R [options]\n\n",
    "Core functions live in sstruth (ssRuns, SS_outputs, ssBindRuns, ...).\n\n",
    "Options:\n",
    "  --study ID                 Study bundle (default: sma-2026)\n",
    "  --appendix-dir DIR         Appendices root (default: appendices/)\n",
    "  --ss-base PATH             SS3 scenario parent directory\n",
    "  --scenarios A,B,C          Comma-separated scenario ids (under ss-base)\n",
    "  --include-sensitivity      Discover all Report.sso runs under ss-base\n",
    "  --discover                 Scan ss-base for Report.sso\n",
    "  --no-discover              Use study spec + sensitivity_register\n",
    "  --collect                  Phase 1: cache ss_output.rds per run folder\n",
    "  --from-cache               Phase 2: load cache only\n",
    "  --slot NAME                SS_output slot, e.g. Kobe, timeseries (--element alias)\n",
    "  --parallel / --no-parallel Parallel SS_output read (default: parallel)\n",
    "  --workers N                Parallel workers\n",
    "  --out FILE                 Write CSV output\n",
    "  --scan DIR                 List Report.sso run directories\n",
    sep = ""
  )
}

resolveAppendixDir <- function(path = NULL) {
  if (!is.null(path) && nzchar(path) && dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  env <- Sys.getenv("DIAG_APPENDIX_CORE", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  for (candidate in c(
    file.path(getwd(), "appendices"),
    normalizePath("appendices", winslash = "/", mustWork = FALSE)
  )) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot find appendices/ directory.", call. = FALSE)
}

loadStudyPaths <- function(
  appendixDir,
  studyId,
  includeSensitivity = FALSE,
  discover = NULL
) {
  requireSstruth()
  out <- sstruth::ssStudyRuns(
    studyId = studyId,
    appendixDir = appendixDir,
    discover = discover,
    includeSensitivity = includeSensitivity
  )
  list(
    ssBase = out$base,
    scenarios = out$scenarios,
    runs = out$runs
  )
}

resolveStudyRuns <- function(parsed, appendixDir) {
  study <- NULL
  runs <- NULL
  scenarios <- parsed$scenarios

  if (!is.null(parsed$ssBase) && nzchar(parsed$ssBase)) {
    ssBase <- normalizePath(parsed$ssBase, winslash = "/", mustWork = FALSE)
    if (is.null(scenarios) || !length(scenarios)) {
      study <- loadStudyPaths(
        appendixDir,
        parsed$studyId,
        includeSensitivity = isTRUE(parsed$includeSensitivity),
        discover = parsed$discover
      )
      scenarios <- study$scenarios
      runs <- study$runs
    }
  } else {
    study <- loadStudyPaths(
      appendixDir,
      parsed$studyId,
      includeSensitivity = isTRUE(parsed$includeSensitivity),
      discover = parsed$discover
    )
    ssBase <- study$ssBase
    scenarios <- scenarios %||% study$scenarios
    runs <- study$runs
  }

  if (is.null(ssBase)) {
    ssBase <- study$ssBase %||% ""
  }

  if (isTRUE(parsed$discover) && (is.null(runs) || !NROW(runs)) &&
      nzchar(ssBase) && dir.exists(ssBase)) {
    runs <- sstruth::ssRuns(ssBase)
    scenarios <- runs$id
  }

  if ((is.null(runs) || !NROW(runs)) && length(scenarios)) {
    runs <- sstruth::ssRuns(ssBase, scenarios)
  }

  if (is.null(scenarios) || !length(scenarios)) {
    stop("No scenarios resolved; use --scenarios, --study, or --scan.", call. = FALSE)
  }

  list(study = study, ssBase = ssBase, scenarios = scenarios, runs = runs)
}

collectStudy <- function(
  appendixDir = "appendices",
  studyId = "sma-2026",
  includeSensitivity = FALSE,
  discover = NULL,
  manifest = NULL,
  ...
) {
  requireSstruth()
  appendixDir <- resolveAppendixDir(appendixDir)
  study <- loadStudyPaths(
    appendixDir,
    studyId,
    includeSensitivity = includeSensitivity,
    discover = discover
  )
  status <- sstruth::SS_outputs(study$runs, ...)
  if (!is.null(manifest) && nzchar(manifest)) {
    dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(status, manifest, row.names = FALSE)
    message("[collect] manifest -> ", manifest)
  }
  structure(list(study = study, status = status), class = "ssCollectStudy")
}

extractStudy <- function(
  appendixDir = "appendices",
  studyId = "sma-2026",
  slot,
  includeSensitivity = FALSE,
  discover = NULL,
  fromCache = FALSE,
  ...
) {
  requireSstruth()
  appendixDir <- resolveAppendixDir(appendixDir)
  study <- loadStudyPaths(
    appendixDir,
    studyId,
    includeSensitivity = includeSensitivity,
    discover = discover
  )
  df <- sstruth::ssExtract(
    base = study$ssBase,
    ids = study$scenarios,
    slot = slot,
    runs = study$runs,
    cache = fromCache,
    ...
  )
  structure(list(study = study, data = df), class = "ssExtractStudy")
}

main <- function() {
  requireSstruth()
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args) {
    extractHelp()
    return(invisible(NULL))
  }

  parsed <- parseExtractArgs(args)

  if (!is.null(parsed$scan) && nzchar(parsed$scan)) {
    df <- sstruth::ssRuns(parsed$scan)
    if (!is.null(parsed$out) && nzchar(parsed$out)) {
      dir.create(dirname(parsed$out), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(df, parsed$out, row.names = FALSE)
      message("[scan] found ", nrow(df), " runs -> ", parsed$out)
    } else if (nrow(df)) {
      writeLines(df$path)
    }
    return(invisible(df))
  }

  appendixDir <- resolveAppendixDir(parsed$appendixDir)
  ctx <- resolveStudyRuns(parsed, appendixDir)

  if (isTRUE(parsed$collect)) {
    status <- sstruth::SS_outputs(
      ctx$runs,
      covar = parsed$covar,
      forecast = parsed$forecast,
      refresh = parsed$refresh,
      parallel = parsed$parallel,
      workers = parsed$workers
    )
    if (!is.null(parsed$out) && nzchar(parsed$out)) {
      dir.create(dirname(parsed$out), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(status, parsed$out, row.names = FALSE)
      message("[collect] manifest -> ", parsed$out)
    } else {
      print(status)
    }
    return(invisible(list(study = ctx$study, runs = ctx$runs, status = status)))
  }

  if (isTRUE(parsed$listSlots)) {
    reps <- sstruth::ssLoad(
      ctx$runs,
      cache = parsed$fromCache,
      covar = parsed$covar,
      forecast = parsed$forecast,
      refresh = parsed$refresh,
      parallel = parsed$parallel,
      workers = parsed$workers
    )
    print(sstruth::ssSlots(reps[[1]]))
    return(invisible(reps))
  }

  if (is.null(parsed$slot) || !nzchar(parsed$slot)) {
    stop(
      "Provide --slot (e.g. Kobe) or use --collect / --list-slots.",
      call. = FALSE
    )
  }

  df <- sstruth::ssExtract(
    base = ctx$ssBase,
    ids = ctx$scenarios,
    slot = parsed$slot,
    runs = ctx$runs,
    cache = parsed$fromCache,
    covar = parsed$covar,
    forecast = parsed$forecast,
    refresh = parsed$refresh,
    parallel = parsed$parallel,
    workers = parsed$workers,
    col = parsed$scenarioCol %||% "scenario"
  )

  if (!is.null(parsed$out) && nzchar(parsed$out)) {
    dir.create(dirname(parsed$out), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, parsed$out, row.names = FALSE)
    message("[extract] wrote ", nrow(df), " rows -> ", parsed$out)
  } else {
    print(df)
  }

  invisible(list(data = df, runs = ctx$runs))
}

if (sys.nframe() == 0L) {
  main()
}
