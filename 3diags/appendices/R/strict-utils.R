# Fail-fast helpers for the 3diags diagnostic workflow.
# Set strict_errors = FALSE (or DIAG_STRICT_ERRORS=false) only for exploratory runs.

diag_strict <- function(cfg = NULL) {
  if (!is.null(cfg) && !is.null(cfg$strict_errors)) {
    return(isTRUE(cfg$strict_errors))
  }
  val <- tolower(Sys.getenv("DIAG_STRICT_ERRORS", unset = "true"))
  val %in% c("1", "true", "yes")
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package '", pkg, "' is not installed.", call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  invisible(TRUE)
}

require_ok <- function(condition, message, cfg = NULL) {
  ok <- if (length(condition) == 1L && is.numeric(condition)) {
    is.finite(condition) && condition > 0
  } else {
    isTRUE(condition)
  }
  if (ok) {
    return(invisible(TRUE))
  }
  if (diag_strict(cfg)) {
    stop(message, call. = FALSE)
  }
  FALSE
}

require_data <- function(x, label, cfg = NULL, min_rows = 1L) {
  ok <- !is.null(x) && (is.data.frame(x) || is.matrix(x)) && NROW(x) >= min_rows
  require_ok(
    ok,
    paste0(label, " is required but missing or empty"),
    cfg = cfg
  )
  ok
}

read_rds_strict <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Required RDS not found: ", path, " (", label, ")", call. = FALSE)
  }
  readRDS(path)
}

read_json_strict <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Required JSON not found: ", path, " (", label, ")", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read ", path, call. = FALSE)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

stop_if <- function(condition, message) {
  if (isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  invisible(NULL)
}

# In strict mode: stop with `message`. In soft mode: return FALSE (caller may return NULL).
soft_ok <- function(ok, message, cfg = NULL) {
  require_ok(ok, message, cfg = cfg)
}
