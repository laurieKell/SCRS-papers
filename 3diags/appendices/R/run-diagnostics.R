#!/usr/bin/env Rscript
#
# Run modular SS3 diagnostic tasks in parallel where possible.
#
# First pass (local profile re-fit + key figures):
#   Rscript appendices/R/run-diagnostics.R \
#     --appendix-dir appendices \
#     --study sma-2026 \
#     --ss-base "P:/rfmo/iccat/2026/shark/Analysis/ss3" \
#     --scenarios Scenario_1/Sc1_Base,Scenario_2/Scn02_base,Scenario_3/Sc3_Base \
#     --profile-scenario Scenario_1/Sc1_Base \
#     --tasks first_pass \
#     --workers 6
#
# Import pre-run profile grid (skip SS3 profile runs):
#   Rscript appendices/R/run-diagnostics.R \
#     --appendix-dir appendices \
#     --study sma-natl \
#     --import-profiles C:/path/to/ss-ll/00_Model_run_files \
#     --tasks first_pass_import
#
# Plot only from existing outputs:
#   Rscript appendices/R/run-diagnostics.R \
#     --appendix-dir appendices \
#     --study sma-natl \
#     --plot-only
#
# Environment: SS3_EXE, DIAG_SS_BASE, DIAG_PROFILE_M, DIAG_PROFILE_BETA

diag_appendix_dir_from_args <- function(args) {
  out <- "appendices"
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--appendix-dir" && i < length(args)) {
      out <- args[[i + 1L]]
    }
    i <- i + 1L
  }
  out
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- diag_appendix_dir_from_args(args)

  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript appendices/R/run-diagnostics.R --appendix-dir DIR [options]\n\n")
    cat("Options:\n")
    cat("  --study ID              Case study (default: sma-2026)\n")
    cat("  --ss-base PATH          SS3 scenario parent directory\n")
    cat("  --scenarios A,B,C       Comma-separated scenario folder names\n")
    cat("  --profile-scenario ID   Reference for M x beta grid (default: study reference)\n")
    cat("  --import-profiles PATH  Import profile grid instead of running SS3\n")
    cat("  --tasks PRESET|t1,t2    Task preset or comma-separated task ids\n")
    cat("  --plot-only             figures preset only (no SS3 runs)\n")
    cat("  --workers N             Parallel workers (default: detectCores - 2)\n")
    cat("  --no-parallel           Disable parallel SS3 profile runs\n")
    cat("  --parallel-curvess      Parallel curveSS (Linux; can hang on Windows)\n")
    cat("  --maxY N                curveSS maxY for profile production plot\n")
    cat("\n")
    tasks_path <- file.path(diags::appendixCoreDir(appendix_dir), "R", "diag-tasks.R")
    if (file.exists(tasks_path)) {
      source(tasks_path, local = TRUE)
      diag_task_help()
    }
    return(invisible(NULL))
  }

  tasks_path <- file.path(diags::appendixCoreDir(appendix_dir), "R", "diag-tasks.R")
  if (!file.exists(tasks_path)) {
    stop("Cannot find ", tasks_path, call. = FALSE)
  }
  source(tasks_path, local = FALSE)

  parsed <- parse_diag_args(args)
  if (is.null(parsed$appendix_dir)) {
    parsed$appendix_dir <- appendix_dir
  }
  ctx <- init_diag_context(parsed)
  run_diag_tasks(ctx)
}

if (sys.nframe() == 0L) {
  main()
}
