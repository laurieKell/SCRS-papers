#!/usr/bin/env Rscript
#
# 00 — i) list scenarios  ii) SS_outputs  iii) readFLSss3
#
#   Rscript appendices/R/00_scenarios.R --ss-base /path/to/results --workers 4
#
# In R:
#   sstruth::ssScenarios(base = "/path/to/results")

if (!requireNamespace("sstruth", quietly = TRUE)) {
  stop("Package 'sstruth' is required.", call. = FALSE)
}

sstruth::ssScenariosCli()
