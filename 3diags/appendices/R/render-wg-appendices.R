#!/usr/bin/env Rscript

# Render SCRS/2026/143 working-group appendices from 2026 SS3 runs.
#
# Usage:
#   Rscript appendices/R/render-wg-appendices.R \
#     --appendix-dir appendices \
#     --study sma-2026 \
#     --ss-base /path/to/new/scenario/parent \
#     --run-code --compile
#
# Environment:
#   DIAG_SS_BASE   SS3 scenario parent (required if paths not in study spec)
#   DIAG_STUDY     Default sma-2026
#   SS3_EXE        SS3 executable

parse_wg_args <- function(args) {
  out <- list(
    appendix_dir = NULL,
    study_id = "sma-2026",
    ss_base = NULL,
    run_code = FALSE,
    compile = FALSE,
    tasks = NULL
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
    } else if (key == "--ss-base" && i < length(args)) {
      out$ss_base <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--run-code") {
      out$run_code <- TRUE
      i <- i + 1L
    } else if (key == "--compile") {
      out$compile <- TRUE
      i <- i + 1L
    } else if (key == "--tasks" && i < length(args)) {
      out$tasks <- strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
      out$tasks <- trimws(out$tasks)
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

compile_wg_pdf <- function(project_dir) {
  tex <- file.path(project_dir, "SCRS_2026_143_appendices.tex")
  if (!file.exists(tex)) {
    stop("Cannot find ", tex, call. = FALSE)
  }
  owd <- setwd(project_dir)
  on.exit(setwd(owd), add = TRUE)
  for (cmd in c(
    "pdflatex -interaction=nonstopmode SCRS_2026_143_appendices.tex",
    "bibtex SCRS_2026_143_appendices",
    "pdflatex -interaction=nonstopmode SCRS_2026_143_appendices.tex",
    "pdflatex -interaction=nonstopmode SCRS_2026_143_appendices.tex"
  )) {
    message("[compile] ", cmd)
    status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (status != 0L) {
      stop("LaTeX failed: ", cmd, call. = FALSE)
    }
  }
  pdf <- file.path(project_dir, "SCRS_2026_143_appendices.pdf")
  message("[compile] wrote ", normalizePath(pdf, winslash = "/"))
  invisible(pdf)
}

main <- function() {
  args <- parse_wg_args(commandArgs(trailingOnly = TRUE))
  appendix_dir <- args$appendix_dir
  if (is.null(appendix_dir)) {
    for (candidate in c(
      file.path(getwd(), "appendices"),
      normalizePath("appendices", winslash = "/", mustWork = FALSE)
    )) {
      if (dir.exists(candidate)) {
        appendix_dir <- candidate
        break
      }
    }
  }
  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  project_dir <- normalizePath(file.path(appendix_dir, ".."), winslash = "/", mustWork = TRUE)

  if (!is.null(args$ss_base) && nzchar(args$ss_base)) {
    Sys.setenv(DIAG_SS_BASE = normalizePath(args$ss_base, winslash = "/", mustWork = FALSE))
  }
  if (!nzchar(Sys.getenv("DIAG_SS_BASE", unset = ""))) {
    cfg <- case_study_defaults(appendix_dir = appendix_dir, study_id = args$study_id)
    if (dir.exists(cfg$ss_base)) {
      Sys.setenv(DIAG_SS_BASE = cfg$ss_base)
    }
  }
  if (!nzchar(Sys.getenv("DIAG_SS_BASE", unset = ""))) {
    stop(
      "Set --ss-base or DIAG_SS_BASE to the 2026 shark WG SS3 scenario parent ",
      "(expected P:/rfmo/iccat/2026/shark/Analysis/ss3).",
      call. = FALSE
    )
  }

  Sys.setenv(DIAG_STUDY = args$study_id)

  if (!is.null(args$tasks) && length(args$tasks)) {
    diag_script <- file.path(appendix_dir, "R", "run-diagnostics.R")
    if (!file.exists(diag_script)) {
      stop("Cannot find ", diag_script, call. = FALSE)
    }
    task_arg <- paste(args$tasks, collapse = ",")
    cmd <- c(
      diag_script,
      "--appendix-dir", appendix_dir,
      "--study", args$study_id,
      "--ss-base", Sys.getenv("DIAG_SS_BASE"),
      "--tasks", task_arg
    )
    message("[wg] running diagnostics: ", task_arg)
    status <- system2(normalizePath(Sys.which("Rscript"), winslash = "/"), cmd)
    if (status != 0L) {
      stop("Diagnostic tasks failed.", call. = FALSE)
    }
  }

  render_script <- file.path(appendix_dir, "render_appendices.R")
  if (!file.exists(render_script)) {
    stop("Cannot find ", render_script, call. = FALSE)
  }
  message("[wg] rendering appendices (study=", args$study_id, ")")
  source(render_script, local = TRUE)
  render_appendices(
    run_code = isTRUE(args$run_code),
    appendix_dir = appendix_dir,
    study_id = args$study_id
  )

  if (isTRUE(args$compile)) {
    compile_wg_pdf(project_dir)
  }
}

if (sys.nframe() == 0L) {
  main()
}
