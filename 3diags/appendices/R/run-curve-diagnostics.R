#!/usr/bin/env Rscript

# Render SS3 curve diagnostic panel (Kobe, production function, SP, PE).
#
# Usage (from 3diags/):
#   Rscript appendices/R/run-curve-diagnostics.R
#   Rscript appendices/R/run-curve-diagnostics.R --appendix-dir appendices --study sma-2026
#
# Default SS3 base: P:/rfmo/iccat/2026/shark/Analysis/ss3
# (via appendices/studies/sma-2026.R). Override with DIAG_SS_BASE.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- normalizePath(
    Sys.getenv("DIAG_APPENDIX_DIR", unset = file.path(getwd(), "appendices")),
    winslash = "/",
    mustWork = FALSE
  )
  study_id <- Sys.getenv("DIAG_STUDY", unset = "sma-2026")

  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--appendix-dir" && i < length(args)) {
      appendix_dir <- normalizePath(args[[i + 1L]], winslash = "/", mustWork = TRUE)
      i <- i + 2L
    } else if (args[[i]] == "--study" && i < length(args)) {
      study_id <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }

  if (!dir.exists(appendix_dir)) {
    stop("appendix_dir not found: ", appendix_dir, call. = FALSE)
  }

  Sys.setenv(DIAG_APPENDIX_DIR = appendix_dir)
  Sys.setenv(DIAG_STUDY = study_id)
  Sys.setenv(DIAG_PROJECT_DIR = dirname(appendix_dir))

  rmd <- file.path(appendix_dir, "Rmd", "ss-curve-diagnostics.Rmd")
  if (!file.exists(rmd)) {
    stop("Rmd not found: ", rmd, call. = FALSE)
  }

  message("Rendering ", rmd)
  message("Study: ", study_id)
  sstruth_root <- Sys.getenv("SSTRUTH_ROOT", unset = "C:/active/flr/sstruth")
  if (dir.exists(sstruth_root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(sstruth_root, quiet = TRUE)
  }
  Sys.setenv(SSTRUTH_PARALLEL = "false")

  out <- rmarkdown::render(
    rmd,
    output_dir = file.path(appendix_dir, "Rmd"),
    quiet = FALSE
  )
  message("Done: ", out)
  invisible(out)
}

if (sys.nframe() == 0L) {
  main()
}
