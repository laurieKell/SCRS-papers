#!/usr/bin/env Rscript

render_appendix <- function(
  strand = c("internal", "external", "robustness"),
  run_code = FALSE,
  launch_profiles = FALSE,
  appendix_dir = NULL,
  study_id = NULL
) {
  strand <- match.arg(strand)
  template <- paste0(strand, "-diagnostics.Rmd")

  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE)[1])

  if (is.null(appendix_dir)) {
    appendix_dir <- if (!is.na(file_arg) && nzchar(file_arg)) {
      dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
    } else {
      normalizePath("appendices", winslash = "/", mustWork = FALSE)
    }
  }

  source(file.path(appendix_dir, "R", "appendix-utils.R"), local = TRUE)
  render_params <- merge_params(
    overrides = list(
      run_code = isTRUE(run_code),
      launch_profiles = isTRUE(launch_profiles)
    ),
    appendix_dir = appendix_dir,
    study_id = study_id
  )

  input <- file.path(appendix_dir, "templates", template)
  output_dir <- file.path(appendix_dir, "content", strand)
  output_file <- "generated.tex"

  if (!file.exists(input)) {
    stop("Cannot find template: ", input, call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render(
    input = input,
    output_format = rmarkdown::latex_fragment(),
    params = render_params,
    output_dir = output_dir,
    output_file = output_file,
    quiet = FALSE
  )
  unlink(file.path(output_dir, output_file))

  invisible(render_params)
}

render_appendices <- function(run_code = FALSE, launch_profiles = FALSE, appendix_dir = NULL, study_id = NULL) {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE)[1])

  if (is.null(appendix_dir)) {
    appendix_dir <- if (!is.na(file_arg) && nzchar(file_arg)) {
      dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
    } else {
      normalizePath("appendices", winslash = "/", mustWork = FALSE)
    }
  }

  source(file.path(appendix_dir, "R", "appendix-utils.R"), local = TRUE)
  render_params <- merge_params(
    overrides = list(
      run_code = isTRUE(run_code),
      launch_profiles = isTRUE(launch_profiles)
    ),
    appendix_dir = appendix_dir,
    study_id = study_id
  )
  templates <- c(
    "internal-diagnostics.Rmd",
    "external-diagnostics.Rmd",
    "robustness-diagnostics.Rmd"
  )

  for (template in templates) {
    render_appendix(
      strand = sub("-diagnostics.Rmd$", "", template),
      run_code = run_code,
      launch_profiles = launch_profiles,
      appendix_dir = appendix_dir,
      study_id = study_id
    )
  }

  invisible(render_params)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_code <- "--run-code" %in% args
  study_id <- NULL
  i <- which(args == "--study")
  if (length(i)) {
    study_id <- args[[i[[1L]] + 1L]]
  }
  render_appendices(run_code = run_code, study_id = study_id)
}
