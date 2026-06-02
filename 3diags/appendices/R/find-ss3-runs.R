#!/usr/bin/env Rscript
#
# List all SS3 run directories (containing Report.sso) under a folder.
#
# Usage:
#   Rscript appendices/R/find-ss3-runs.R "P:/rfmo/iccat/2026/shark/Analysis/ss3"
#   Rscript appendices/R/find-ss3-runs.R "P:/..." --out runs.csv
#
# In R:
#   sstruth::ssRuns("P:/.../ss3")

parseFindArgs <- function(args) {
  out <- list(root = NULL, out = NULL)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--out" && i < length(args)) {
      out$out <- args[[i + 1L]]
      i <- i + 2L
    } else if (is.null(out$root) && !startsWith(args[[i]], "-")) {
      out$root <- args[[i]]
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args || !length(args)) {
    cat(
      "Usage: Rscript appendices/R/find-ss3-runs.R ROOT [options]\n\n",
      "Uses sstruth::ssRuns().\n\n",
      "Options:\n",
      "  --out FILE    Write id,path CSV (default: print absolute paths)\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required.", call. = FALSE)
  }

  parsed <- parseFindArgs(args)
  if (is.null(parsed$root) || !nzchar(parsed$root)) {
    stop("Provide ROOT folder to scan.", call. = FALSE)
  }

  df <- sstruth::ssRuns(parsed$root)
  if (!is.null(parsed$out) && nzchar(parsed$out)) {
    dir.create(dirname(parsed$out), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, parsed$out, row.names = FALSE)
    message("[scan] found ", nrow(df), " runs -> ", parsed$out)
  } else if (nrow(df)) {
    writeLines(df$path)
  } else {
    message("[scan] no Report.sso found under ", parsed$root)
  }
  invisible(df)
}

if (sys.nframe() == 0L) {
  main()
}
