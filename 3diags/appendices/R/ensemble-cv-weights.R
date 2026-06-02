# Ensemble weights from cross-validation / hindcast prediction skill.
#
# Public API:
#   CVIndexmase() / CVIndexrho()           — single cv_data table → scores + weights
#   CVIndexmase_combine() / CVIndexrho_combine() — list of cv tables → combined weights
#
# cv_data columns: model, year, observed, predicted (optional: peel for ordering).

cv_validate_data <- function(cv_data) {
  req <- c("model", "year", "observed", "predicted")
  if (is.null(cv_data) || !NROW(cv_data)) {
    return(NULL)
  }
  if (!all(req %in% names(cv_data))) {
    stop("cv_data must contain: ", paste(req, collapse = ", "), call. = FALSE)
  }
  cv_data[is.finite(cv_data$observed) & is.finite(cv_data$predicted), , drop = FALSE]
}

cv_naive_error <- function(dat) {
  if ("peel" %in% names(dat)) {
    dat <- dat[order(dat$peel), , drop = FALSE]
  } else {
    dat <- dat[order(dat$year, seq_len(nrow(dat))), , drop = FALSE]
  }
  obs <- dat$observed
  pred <- dat$predicted
  if (length(obs) < 2L) {
    return(NA_real_)
  }
  err <- abs(obs - c(pred[1], pred[-length(pred)]))
  if (!any(is.finite(err)) || max(err, na.rm = TRUE) <= 0) {
    err <- abs(diff(obs))
  }
  if (!any(is.finite(err)) || max(err, na.rm = TRUE) <= 0) {
    err <- abs(obs - stats::median(obs, na.rm = TRUE))
  }
  if (!any(is.finite(err)) || max(err, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }
  mean(err, na.rm = TRUE)
}

cv_calc_mase <- function(cv_data) {
  dat <- cv_validate_data(cv_data)
  if (is.null(dat) || !NROW(dat)) {
    return(NULL)
  }
  mae <- stats::aggregate(
    abs(dat$observed - dat$predicted),
    by = list(model = dat$model),
    FUN = mean,
    na.rm = TRUE
  )
  names(mae)[2] <- "mae"
  naive <- lapply(split(dat, dat$model), cv_naive_error)
  mae$mean_naive_error <- vapply(mae$model, function(m) naive[[m]], numeric(1))
  mae$mase <- mae$mae / mae$mean_naive_error
  bad <- !is.finite(mae$mase) | !is.finite(mae$mean_naive_error) | mae$mean_naive_error <= 0
  mae$mase[bad] <- mae$mae[bad]
  mae[, c("model", "mae", "mean_naive_error", "mase")]
}

cv_calc_rho <- function(cv_data) {
  dat <- cv_validate_data(cv_data)
  if (is.null(dat) || !NROW(dat)) {
    return(NULL)
  }
  dat <- dat[dat$observed != 0, , drop = FALSE]
  if (!NROW(dat)) {
    return(NULL)
  }
  stats::aggregate(
    (dat$predicted - dat$observed) / dat$observed,
    by = list(model = dat$model),
    FUN = mean,
    na.rm = TRUE
  ) -> out
  names(out)[2] <- "rho"
  out$abs_rho <- abs(out$rho)
  out
}
cv_optimize_linear_weights <- function(values, model_names = NULL, abs_values = FALSE) {
  v <- as.numeric(values)
  n <- length(v)
  if (n < 1L) {
    return(numeric(0))
  }
  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_len(n))
  }
  if (n == 1L) {
    return(stats::setNames(1, model_names[[1]]))
  }
  v <- if (isTRUE(abs_values)) abs(v) else v
  if (any(!is.finite(v))) {
    v[!is.finite(v)] <- max(v[is.finite(v)], na.rm = TRUE) * 10
  }
  if (!any(is.finite(v))) {
    w <- rep(1 / n, n)
    return(stats::setNames(w, model_names))
  }

  solve_weights <- function() {
    if (requireNamespace("nloptr", quietly = TRUE)) {
      objective <- function(w) sum(w * v)
      constraint <- function(w) sum(w) - 1
      result <- nloptr::nloptr(
        x0 = rep(1 / n, n),
        eval_f = objective,
        lb = rep(0, n),
        ub = rep(1, n),
        eval_g_eq = constraint,
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 1e-8,
          maxeval = 10000
        )
      )
      w <- result$solution
      w[w < 0] <- 0
      if (sum(w) <= 0) {
        w <- rep(0, n)
        w[which.min(v)] <- 1
      } else {
        w <- w / sum(w)
      }
      return(list(weights = stats::setNames(w, model_names), convergence = result$status))
    }
    w <- rep(0, n)
    w[which.min(v)] <- 1
    list(weights = stats::setNames(w, model_names), convergence = NA_integer_)
  }

  out <- solve_weights()
  out$objective_value <- sum(out$weights * v)
  out
}

cv_optimize_mase_weights <- function(mase_tbl) {
  if (is.null(mase_tbl) || !NROW(mase_tbl)) {
    return(NULL)
  }
  ok <- is.finite(mase_tbl$mase)
  if (!any(ok)) {
    return(NULL)
  }
  tbl <- mase_tbl[ok, , drop = FALSE]
  cv_optimize_linear_weights(tbl$mase, model_names = as.character(tbl$model))
}

cv_optimize_rho_weights <- function(rho_tbl) {
  if (is.null(rho_tbl) || !NROW(rho_tbl)) {
    return(NULL)
  }
  val_col <- if ("abs_rho" %in% names(rho_tbl)) "abs_rho" else "rho"
  ok <- is.finite(rho_tbl[[val_col]])
  if (!any(ok)) {
    return(NULL)
  }
  tbl <- rho_tbl[ok, , drop = FALSE]
  cv_optimize_linear_weights(tbl[[val_col]], model_names = as.character(tbl$model), abs_values = FALSE)
}

#' CVindexmase: MASE scores and optional ensemble weights for one CV table
#'
#' MASE = mean(|observed - predicted|) / mean(|naive error|). Naive forecast uses
#' lagged predictions within each model (or peel order when \code{peel} is present).
#'
#' @param cv_data data.frame with \code{model}, \code{year}, \code{observed}, \code{predicted}.
#' @param optimize If \code{TRUE}, minimize \eqn{\sum w_i \mathrm{MASE}_i} subject to
#'   \eqn{w_i \ge 0}, \eqn{\sum w_i = 1} (requires \pkg{nloptr} for general case).
#' @return List with \code{index}, \code{scores}, and optionally \code{weights},
#'   \code{objective_value}, \code{convergence}.
CVIndexmase <- function(cv_data, optimize = TRUE) {
  scores <- cv_calc_mase(cv_data)
  if (is.null(scores)) {
    return(NULL)
  }
  out <- list(index = "CVindexmase", scores = scores)
  if (isTRUE(optimize)) {
    opt <- cv_optimize_mase_weights(scores)
    if (!is.null(opt)) {
      out$weights <- opt$weights
      out$objective_value <- opt$objective_value
      out$convergence <- opt$convergence
    }
  }
  out
}

#' CVindexrho: retrospective-style rho and optional ensemble weights for one CV table
#'
#' \eqn{\rho = \mathrm{mean}((\mathrm{predicted} - \mathrm{observed}) / \mathrm{observed})}
#' per model. Weights minimize \eqn{\sum w_i |\rho_i|}.
#'
#' @inheritParams CVIndexmase
CVIndexrho <- function(cv_data, optimize = TRUE) {
  scores <- cv_calc_rho(cv_data)
  if (is.null(scores)) {
    return(NULL)
  }
  out <- list(index = "CVindexrho", scores = scores)
  if (isTRUE(optimize)) {
    opt <- cv_optimize_rho_weights(scores)
    if (!is.null(opt)) {
      out$weights <- opt$weights
      out$objective_value <- opt$objective_value
      out$convergence <- opt$convergence
    }
  }
  out
}

#' Combine CVindexmase weights across multiple CV tables (e.g. CPUE, catch, depletion)
#'
#' @param cv_list Named list of cv_data frames; \code{NULL} entries are skipped.
#' @param combine_method \code{"average"}, \code{"product"}, or \code{"min_var"}.
CVIndexmase_combine <- function(cv_list, combine_method = "average", optimize = TRUE) {
  nms <- names(cv_list)
  if (is.null(nms) || !length(nms)) {
    nms <- paste0("metric_", seq_along(cv_list))
  }
  keep <- !vapply(cv_list, is.null, logical(1))
  cv_list <- cv_list[keep]
  nms <- nms[keep]
  if (!length(cv_list)) {
    return(NULL)
  }
  by_metric <- stats::setNames(
    lapply(cv_list, function(x) CVIndexmase(x, optimize = optimize)),
    nms
  )
  by_metric <- by_metric[!vapply(by_metric, is.null, logical(1))]
  if (!length(by_metric)) {
    return(NULL)
  }
  weight_objs <- lapply(by_metric, function(x) {
    if (is.null(x$weights)) {
      return(NULL)
    }
    list(weights = x$weights)
  })
  weight_objs <- weight_objs[!vapply(weight_objs, is.null, logical(1))]
  combined <- if (length(weight_objs)) {
    cv_combine_weights(weight_objs, combine_method)
  } else {
    NULL
  }
  scores <- lapply(names(by_metric), function(nm) {
    s <- by_metric[[nm]]$scores
    s$metric <- nm
    s
  })
  list(
    index = "CVindexmase",
    by_metric = by_metric,
    scores = do.call(rbind, scores),
    weights = combined,
    combine_method = combine_method
  )
}

#' Combine CVindexrho weights across multiple CV tables
#' @inheritParams CVIndexmase_combine
CVIndexrho_combine <- function(cv_list, combine_method = "average", optimize = TRUE) {
  nms <- names(cv_list)
  if (is.null(nms) || !length(nms)) {
    nms <- paste0("metric_", seq_along(cv_list))
  }
  keep <- !vapply(cv_list, is.null, logical(1))
  cv_list <- cv_list[keep]
  nms <- nms[keep]
  if (!length(cv_list)) {
    return(NULL)
  }
  by_metric <- stats::setNames(
    lapply(cv_list, function(x) CVIndexrho(x, optimize = optimize)),
    nms
  )
  by_metric <- by_metric[!vapply(by_metric, is.null, logical(1))]
  if (!length(by_metric)) {
    return(NULL)
  }
  weight_objs <- lapply(by_metric, function(x) {
    if (is.null(x$weights)) {
      return(NULL)
    }
    list(weights = x$weights)
  })
  weight_objs <- weight_objs[!vapply(weight_objs, is.null, logical(1))]
  combined <- if (length(weight_objs)) {
    cv_combine_weights(weight_objs, combine_method)
  } else {
    NULL
  }
  scores <- lapply(names(by_metric), function(nm) {
    s <- by_metric[[nm]]$scores
    s$metric <- nm
    s
  })
  list(
    index = "CVindexrho",
    by_metric = by_metric,
    scores = do.call(rbind, scores),
    weights = combined,
    combine_method = combine_method
  )
}

cv_combine_weights <- function(weight_list, method = "average") {
  weight_list <- weight_list[!vapply(weight_list, is.null, logical(1))]
  if (!length(weight_list)) {
    return(NULL)
  }
  mat <- do.call(rbind, lapply(weight_list, function(x) {
    as.numeric(x$weights)
  }))
  model_names <- names(weight_list[[1]]$weights)
  if (is.null(model_names) || length(model_names) != ncol(mat)) {
    model_names <- paste0("Model_", seq_len(ncol(mat)))
  }
  colnames(mat) <- model_names
  if (nrow(mat) == 1L) {
    w <- mat[1, ]
  } else if (method == "average") {
    w <- colMeans(mat)
  } else if (method == "product") {
    w <- apply(mat, 2, prod)
  } else if (method == "min_var") {
    vw <- apply(mat, 2, stats::var)
    vw[!is.finite(vw) | vw <= 0] <- 1
    w <- colMeans(mat / vw)
  } else {
    stop("Unknown combine method: ", method, call. = FALSE)
  }
  w <- w / sum(w)
  stats::setNames(w, model_names)
}

cv_data_from_retro_one_step <- function(one_step, metric = c("catch", "f", "depletion")) {
  metric <- match.arg(metric)
  if (is.null(one_step) || !NROW(one_step)) {
    return(NULL)
  }
  peel_ord <- one_step$terminalYr
  if (is.null(peel_ord)) {
    peel_ord <- seq_len(nrow(one_step))
  }
  if (metric == "catch") {
    data.frame(
      model = one_step$scenario,
      year = one_step$forecastYr,
      peel = peel_ord,
      observed = one_step$catchEstBaseTime,
      predicted = one_step$catchFore,
      stringsAsFactors = FALSE
    )
  } else if (metric == "f") {
    data.frame(
      model = one_step$scenario,
      year = one_step$forecastYr,
      peel = peel_ord,
      observed = one_step$estFBaseTime,
      predicted = one_step$forecastF,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

cv_data_from_retro_ratio_horizon <- function(ratio_horizon, horizon = 1L) {
  if (is.null(ratio_horizon) || !NROW(ratio_horizon)) {
    return(NULL)
  }
  dat <- ratio_horizon[ratio_horizon$horizon == horizon, , drop = FALSE]
  if (!NROW(dat)) {
    return(NULL)
  }
  data.frame(
    model = dat$scenario,
    year = dat$forecastYr,
    observed = dat$refBratio,
    predicted = dat$foreBratio,
    stringsAsFactors = FALSE
  )
}

cv_data_from_loo_cpue <- function(params, scenario = NULL) {
  cv <- load_collected_crossvalidation(params, scenario = scenario)
  if (is.null(cv) || is.null(cv$results) || !NROW(cv$results)) {
    return(NULL)
  }
  if (!all(c("observed", "predicted", "year", "model") %in% names(cv$results))) {
    return(NULL)
  }
  cv$results[, c("model", "year", "observed", "predicted")]
}

build_scenario_cv_inputs <- function(params, scenarios = NULL) {
  scenarios <- scenarios %||% params$scenarios
  retro <- load_collected_retrospectives(params)
  if (is.null(retro) && retro_dirs_available(params)) {
    retro <- backend_retro_backtest_data(params)
  }
  if (is.null(retro)) {
    if (diag_strict(params)) {
      stop("Retrospective data required for ensemble CV weights.", call. = FALSE)
    }
    return(NULL)
  }

  list(
    cpue = cv_data_from_loo_cpue(params),
    catch = cv_data_from_retro_one_step(retro$oneStep, metric = "catch"),
    f = cv_data_from_retro_one_step(retro$oneStep, metric = "f"),
    depletion = cv_data_from_retro_ratio_horizon(retro$ratioHorizon, horizon = 1L)
  )
}

compute_ensemble_cv_weights <- function(params, combine_method = "average", scenarios = NULL) {
  inputs <- build_scenario_cv_inputs(params, scenarios = scenarios)
  if (is.null(inputs)) {
    return(NULL)
  }

  metrics <- list(
    cpue = inputs$cpue,
    catch = inputs$catch,
    f = inputs$f,
    depletion = inputs$depletion
  )
  metrics <- metrics[!vapply(metrics, is.null, logical(1))]

  cv_index_mase <- CVIndexmase_combine(metrics, combine_method = combine_method)
  cv_index_rho <- CVIndexrho_combine(metrics, combine_method = combine_method)

  mase_by_metric <- lapply(metrics, cv_calc_mase)
  rho_by_metric <- lapply(metrics, cv_calc_rho)
  mase_weights <- lapply(mase_by_metric, cv_optimize_mase_weights)
  rho_weights <- lapply(rho_by_metric, cv_optimize_rho_weights)
  combined_mase <- cv_index_mase$weights
  combined_rho <- cv_index_rho$weights

  summary_rows <- lapply(names(metrics), function(nm) {
    mase_tbl <- mase_by_metric[[nm]]
    rho_tbl <- rho_by_metric[[nm]]
    if (is.null(mase_tbl) && is.null(rho_tbl)) {
      return(NULL)
    }
    models <- unique(c(
      if (!is.null(mase_tbl)) mase_tbl$model else character(0),
      if (!is.null(rho_tbl)) rho_tbl$model else character(0)
    ))
    out <- data.frame(model = models, metric = nm, stringsAsFactors = FALSE)
    if (!is.null(mase_tbl)) {
      out <- merge(out, mase_tbl, by = "model", all.x = TRUE)
    } else {
      out$mase <- NA_real_
    }
    if (!is.null(rho_tbl)) {
      out <- merge(out, rho_tbl, by = "model", all.x = TRUE)
    } else {
      out$rho <- NA_real_
    }
    out
  })
  summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]
  summary_table <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

  weight_table <- NULL
  if (!is.null(summary_table) && NROW(summary_table)) {
    models <- sort(unique(as.character(summary_table$model)))
    weight_table <- data.frame(model = models, stringsAsFactors = FALSE)
    for (nm in names(mase_weights)) {
      wobj <- mase_weights[[nm]]
      col <- paste0("weight_mase_", nm)
      weight_table[[col]] <- if (!is.null(wobj)) wobj$weights[models] else NA_real_
    }
    for (nm in names(rho_weights)) {
      wobj <- rho_weights[[nm]]
      col <- paste0("weight_rho_", nm)
      weight_table[[col]] <- if (!is.null(wobj)) wobj$weights[models] else NA_real_
    }
    if (!is.null(combined_mase)) {
      weight_table$weight_cvindexmase <- combined_mase[models]
    }
    if (!is.null(combined_rho)) {
      weight_table$weight_cvindexrho <- combined_rho[models]
    }
  }

  list(
    inputs = inputs,
    cv_index_mase = cv_index_mase,
    cv_index_rho = cv_index_rho,
    mase = mase_by_metric,
    rho = rho_by_metric,
    weights_mase = mase_weights,
    weights_rho = rho_weights,
    weights_combined_mase = combined_mase,
    weights_combined_rho = combined_rho,
    summary = summary_table,
    weights = weight_table,
    combine_method = combine_method
  )
}

collect_ensemble_cv_weights <- function(params, combine_method = "average", scenarios = NULL) {
  cfg <- mvln_cfg(params)
  out <- compute_ensemble_cv_weights(
    params,
    combine_method = combine_method %||% cfg$weight_combine %||% "average",
    scenarios = scenarios
  )
  if (is.null(out)) {
    if (diag_strict(params)) {
      stop("Ensemble CV weights could not be computed.", call. = FALSE)
    }
    return(NULL)
  }
  root <- longrun_job_root(params, "ensemble-weights")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(out$summary) && NROW(out$summary)) {
    utils::write.csv(out$summary, file.path(root, "cv-metrics.csv"), row.names = FALSE)
  }
  if (!is.null(out$weights) && NROW(out$weights)) {
    utils::write.csv(out$weights, file.path(root, "ensemble_weights.csv"), row.names = FALSE)
  }
  cmp <- ensemble_weight_comparison(params, scenarios = scenarios, cv = out)
  if (!is.null(cmp) && NROW(cmp)) {
    utils::write.csv(cmp, file.path(root, "ensemble_weight_comparison.csv"), row.names = FALSE)
    out$weight_comparison <- cmp
  }
  saveRDS(out, file.path(root, "collected.rds"))
  out$root <- root
  out
}

load_collected_ensemble_cv_weights <- function(params) {
  root <- longrun_job_root(params, "ensemble-weights")
  path <- file.path(root, "collected.rds")
  if (file.exists(path)) {
    if (diag_strict(params)) {
      return(read_rds_strict(path, "ensemble-weights/collected.rds"))
    }
    return(readRDS(path))
  }
  csv <- file.path(root, "ensemble_weights.csv")
  if (file.exists(csv)) {
    return(list(root = root, weights = utils::read.csv(csv, stringsAsFactors = FALSE)))
  }
  if (diag_strict(params)) {
    stop("Ensemble weights not collected: ", path, call. = FALSE)
  }
  NULL
}

#' Total negative log-likelihood (SS3 \code{likelihoods_used TOTAL}) per scenario
scenario_total_nll <- function(params, scenarios = NULL) {
  scenarios <- as.character(scenarios %||% params$scenarios)
  paths <- scenario_paths(params)
  names(paths) <- basename(paths)
  paths <- paths[names(paths) %in% scenarios]
  if (!length(paths)) {
    return(NULL)
  }
  nll <- vapply(paths, function(p) {
    rep <- backend_read_output(params, p)
    if (is.null(rep$likelihoods_used)) {
      if (diag_strict(params)) {
        stop("Likelihoods missing for ", p, call. = FALSE)
      }
      return(NA_real_)
    }
    as.numeric(rep$likelihoods_used["TOTAL", "values"])
  }, numeric(1))
  stats::setNames(nll, names(paths))
}

#' ICES-style information weights from total likelihood (Akaike / benchmark practice)
#'
#' For models with the same number of estimated parameters, weights are
#' \eqn{w_i \propto \exp(-\Delta\mathrm{AIC}/2)} with
#' \eqn{\Delta\mathrm{AIC}_i = 2\,\mathrm{NLL}_i - \min_j 2\,\mathrm{NLL}_j}.
#' Uses \code{sstruth::ensembleWeightsSoftmaxNeg2()}.
ensemble_weights_ices <- function(scenarios, nll = NULL, params = NULL) {
  scenarios <- as.character(scenarios)
  if (is.null(nll)) {
    if (is.null(params)) {
      stop("Provide nll or params.", call. = FALSE)
    }
    nll <- scenario_total_nll(params, scenarios)
  }
  nll <- nll[scenarios]
  if (any(!is.finite(nll))) {
    stop("Non-finite total NLL for scenarios: ", paste(scenarios[!is.finite(nll)], collapse = ", "), call. = FALSE)
  }
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    delta <- 2 * nll - min(2 * nll)
    w <- exp(-0.5 * delta)
    return(stats::setNames(w / sum(w), scenarios))
  }
  stats::setNames(
    sstruth::ensembleWeightsSoftmaxNeg2(2 * nll),
    scenarios
  )
}

#' Wide table comparing ensemble weight schemes
#'
#' Columns: \code{scenario}, \code{w_equal}, \code{w_ices}, \code{w_cvindexmase},
#' \code{w_cvindexrho}, plus \code{total_nll} and \code{delta_aic} for transparency.
ensemble_weight_comparison <- function(params, scenarios = NULL, combine_method = "average", cv = NULL) {
  scenarios <- as.character(scenarios %||% params$scenarios)
  w_equal <- ensemble_scenario_weights(scenarios)
  nll <- scenario_total_nll(params, scenarios)
  w_ices <- ensemble_weights_ices(scenarios, nll = nll)

  if (is.null(cv)) {
    cv <- load_collected_ensemble_cv_weights(params)
    if (is.null(cv) || is.null(cv$cv_index_mase)) {
      cv <- compute_ensemble_cv_weights(params, combine_method = combine_method, scenarios = scenarios)
    }
  }
  w_mase <- if (!is.null(cv) && !is.null(cv$weights_combined_mase)) {
    cv$weights_combined_mase[scenarios]
  } else {
    rep(NA_real_, length(scenarios))
  }
  w_rho <- if (!is.null(cv) && !is.null(cv$weights_combined_rho)) {
    cv$weights_combined_rho[scenarios]
  } else {
    rep(NA_real_, length(scenarios))
  }
  names(w_mase) <- scenarios
  names(w_rho) <- scenarios

  aic <- 2 * nll
  delta_aic <- aic - min(aic, na.rm = TRUE)

  out <- data.frame(
    scenario = scenarios,
    total_nll = as.numeric(nll[scenarios]),
    delta_aic = as.numeric(delta_aic[scenarios]),
    w_equal = as.numeric(w_equal[scenarios]),
    w_ices = as.numeric(w_ices[scenarios]),
    w_cvindexmase = as.numeric(w_mase[scenarios]),
    w_cvindexrho = as.numeric(w_rho[scenarios]),
    stringsAsFactors = FALSE
  )
  attr(out, "cv") <- cv
  out
}

ensemble_weight_comparison_table <- function(comparison) {
  if (is.null(comparison) || !NROW(comparison)) {
    return(NULL)
  }
  out <- comparison[, c("scenario", "w_equal", "w_ices", "w_cvindexmase", "w_cvindexrho"), drop = FALSE]
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) signif(x, 4))
  out
}

resolve_ensemble_weights <- function(params, scenarios) {
  cfg <- mvln_cfg(params)
  scenarios <- as.character(scenarios)
  if (isTRUE(cfg$equal_weights)) {
    return(ensemble_scenario_weights(scenarios))
  }
  cv <- load_collected_ensemble_cv_weights(params)
  if (is.null(cv)) {
    cv <- collect_ensemble_cv_weights(params)
  }
  method <- cfg$weight_method %||% "cvindexmase"
  w <- NULL
  if (!is.null(cv)) {
    if (method == "cvindexrho" && !is.null(cv$weights_combined_rho)) {
      w <- cv$weights_combined_rho[scenarios]
    } else if (!is.null(cv$weights_combined_mase)) {
      w <- cv$weights_combined_mase[scenarios]
    } else if (!is.null(cv$weights) && "weight_cvindexmase" %in% names(cv$weights)) {
      ix <- match(scenarios, cv$weights$model)
      w <- stats::setNames(cv$weights$weight_cvindexmase[ix], scenarios)
    }
  }
  if (!is.null(w)) {
    names(w) <- scenarios
    if (any(!is.finite(w))) {
      w[!is.finite(w)] <- 0
    }
    if (sum(w) > 0) {
      return(ensemble_scenario_weights(scenarios, weights = w))
    }
  }
  require_ok(
    FALSE,
    "CV ensemble weights unavailable; run appendices/R/run-catch-backtests.R --collect and ensemble collect.",
    cfg = params
  )
  ensemble_scenario_weights(scenarios)
}
