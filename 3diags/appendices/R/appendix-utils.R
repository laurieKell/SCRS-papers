# Shared appendix helpers for 3diags diagnostic templates and collectors.

.strict_utils_path <- function() {
  candidates <- c(
    file.path(getwd(), "appendices", "R", "strict-utils.R"),
    file.path(getwd(), "R", "strict-utils.R")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) hits[[1L]] else NULL
}

if (!exists("diag_strict", mode = "function", inherits = FALSE)) {
  path <- .strict_utils_path()
  if (!is.null(path)) {
    source(path, local = FALSE)
  }
}

write_content <- function(params, fragment, lines) {
  path <- file.path(params$appendix_dir, "content", fragment)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x
}

kable_fragment <- function(x, caption = NULL, label = NULL, booktabs = TRUE) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required to format appendix tables.", call. = FALSE)
  }
  if (!is.null(label)) {
    label <- sub("^tab:", "", label)
  }
  out <- knitr::kable(
    x,
    format = "latex",
    booktabs = booktabs,
    caption = caption,
    label = label,
    escape = TRUE
  )
  paste(out, collapse = "\n")
}

figure_fragment <- function(params, image_path, caption, label, width = 0.92) {
  rel <- normalizePath(image_path, winslash = "/", mustWork = TRUE)
  proj <- normalizePath(params$project_dir, winslash = "/", mustWork = TRUE)
  rel_to_project <- sub(
    paste0("^", gsub("([\\W])", "\\\\\\1", proj)),
    "",
    rel
  )
  rel_to_project <- sub("^/", "", rel_to_project)

  cap <- gsub("%", "\\%", caption, fixed = TRUE)
  c(
    "\\begin{figure}[H]",
    "  \\centering",
    sprintf("  \\includegraphics[width=%s\\linewidth]{%s}", width, rel_to_project),
    sprintf("  \\caption{%s}", cap),
    sprintf("  \\label{%s}", label),
    "\\end{figure}"
  )
}

merge_params <- function(overrides = list(), appendix_dir = NULL, study_id = NULL) {
  cfg_path <- if (!is.null(appendix_dir)) {
    file.path(normalizePath(appendix_dir, winslash = "/", mustWork = TRUE), "config-case-study.R")
  } else {
    normalizePath("appendices/config-case-study.R", winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(cfg_path)) {
    stop("Cannot find config-case-study.R at ", cfg_path, call. = FALSE)
  }
  source(cfg_path, local = TRUE)
  study_id <- study_id %||% overrides$study_id %||% Sys.getenv("DIAG_STUDY", unset = "sma-2026")
  defaults <- case_study_defaults(
    appendix_dir = appendix_dir,
    study_id = study_id,
    overrides = overrides
  )
  defaults
}

load_appendix_helpers <- function(params) {
  core <- diags::appendixCoreDir(params$appendix_dir)
  strict_path <- file.path(core, "R", "strict-utils.R")
  if (file.exists(strict_path)) {
    source(strict_path, local = parent.frame())
  }
  study_path <- file.path(params$appendix_dir, "R", "study-config.R")
  if (file.exists(study_path)) {
    source(study_path, local = parent.frame())
  }
  backend_path <- file.path(core, "R", "backend.R")
  if (file.exists(backend_path)) {
    source(backend_path, local = parent.frame())
  }
  source(file.path(params$appendix_dir, "config-case-study.R"), local = parent.frame())
  ss_mvln <- file.path(core, "R", "ss-mvln.R")
  if (file.exists(ss_mvln)) {
    source(ss_mvln, local = parent.frame())
  }
  cfg_mvln <- file.path(params$appendix_dir, "config-mvln.R")
  if (file.exists(cfg_mvln)) {
    source(cfg_mvln, local = parent.frame())
  }
  ens_w <- file.path(core, "R", "ensemble-cv-weights.R")
  if (file.exists(ens_w)) {
    source(ens_w, local = parent.frame())
  }
  invisible(NULL)
}

resolve_appendix_params <- function(params) {
  defaults <- merge_params(
    overrides = list(run_code = isTRUE(params$run_code)),
    appendix_dir = if (!is.null(params$appendix_dir) && nzchar(params$appendix_dir)) {
      params$appendix_dir
    } else {
      NULL
    }
  )
  cfg <- as.list(params)
  for (nm in names(defaults)) {
    val <- cfg[[nm]]
    if (is.null(val) || (is.character(val) && length(val) == 1L && !nzchar(val))) {
      cfg[[nm]] <- defaults[[nm]]
    }
  }
  cfg
}

library_register <- function(packages) {
  stopifnot(is.data.frame(packages))
  required_cols <- c("package", "role", "tier")
  missing_cols <- setdiff(required_cols, names(packages))
  if (length(missing_cols)) {
    stop("library_register() requires columns: ", paste(required_cols, collapse = ", "))
  }
  packages$installed <- vapply(
    packages$package,
    function(pkg) {
      if (pkg %in% c("base", "grDevices", "utils", "stats")) {
        return("built-in")
      }
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        "not installed"
      }
    },
    character(1)
  )
  packages
}

write_libraries <- function(params, strand, packages) {
  libs <- library_register(packages)
  write_content(
    params,
    file.path(strand, "libraries.tex"),
    c(
      kable_fragment(
        libs[, c("package", "role", "tier", "installed")],
        caption = paste("Required R packages for", strand, "diagnostics."),
        label = paste0("tab:", strand, "-libraries")
      ),
      ""
    )
  )
  invisible(libs)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  if (length(x) == 1L && is.na(x)) {
    return(y)
  }
  x
}

safe_require <- function(pkg) {
  if (diag_strict()) {
    require_pkg(pkg)
    return(TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

read_ss3_output <- function(path, covar = FALSE, cfg = NULL) {
  if (!is.null(cfg)) {
    return(backend_read_output(cfg, path, covar = covar))
  }
  if (!safe_require("r4ss")) {
    stop("Package 'r4ss' is required when run_code = TRUE.", call. = FALSE)
  }
  r4ss::SS_output(
    dir = path,
    verbose = FALSE,
    printstats = FALSE,
    forecast = FALSE,
    covar = isTRUE(covar),
    hidewarn = TRUE,
    NoCompOK = TRUE
  )
}

summarise_ss3_run <- function(path, rep = NULL, cfg = NULL) {
  if (!is.null(cfg)) {
    return(backend_summarise_run(cfg, path, rep = rep))
  }
  rep <- rep %||% read_ss3_output(path)
  if (!safe_require("r4ss")) {
    stop("Package 'r4ss' is required.", call. = FALSE)
  }
  info <- r4ss::SSsummarize(list(rep), verbose = FALSE)
  total_ll <- info$likelihoods[
    info$likelihoods$Label == "TOTAL",
    "model1",
    drop = TRUE
  ]
  data.frame(
    scenario = basename(path),
    n_params = info$npars,
    total_ll = round(total_ll, 3),
    max_grad = signif(info$maxgrad, 3),
    log_det_hessian = if (!is.null(rep$log_det_hessian)) {
      signif(rep$log_det_hessian, 4)
    } else {
      NA_real_
    },
    n_iterations = rep$N_iterations %||% NA_integer_,
    n_warnings = rep$Nwarnings %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

likelihood_components_table <- function(rep) {
  ll <- rep$likelihoods_used
  if (is.null(ll) || !NROW(ll)) {
    return(NULL)
  }
  out <- data.frame(
    component = rownames(ll),
    log_likelihood = round(ll[, "values", drop = TRUE], 3),
    lambda = ll[, "lambdas", drop = TRUE],
    stringsAsFactors = FALSE
  )
  out <- out[out$component != "TOTAL", , drop = FALSE]
  out[order(-abs(out$log_likelihood)), , drop = FALSE]
}

parameters_near_bounds_table <- function(rep, tol = 0.02) {
  pars <- rep$parameters
  if (is.null(pars) || !NROW(pars)) {
    return(NULL)
  }
  est <- subset(pars, Phase > 0 & Status %in% c("act", "OK"))
  if (!NROW(est)) {
    return(NULL)
  }
  near <- est[
    (is.finite(est$Min) & is.finite(est$Max) & est$Max > est$Min &
       ((est$Value - est$Min) / (est$Max - est$Min) <= tol |
          (est$Max - est$Value) / (est$Max - est$Min) <= tol)),
    ,
    drop = FALSE
  ]
  if (!NROW(near)) {
    return(data.frame(
      label = character(0),
      value = numeric(0),
      min = numeric(0),
      max = numeric(0),
      status = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    label = near$Label,
    value = signif(near$Value, 4),
    min = near$Min,
    max = near$Max,
    status = near$Status,
    stringsAsFactors = FALSE
  )
}

advice_parameters_table <- function(rep) {
  pars <- rep$parameters
  if (is.null(pars) || !NROW(pars)) {
    return(NULL)
  }
  keep <- grepl(
    "^(SR_|Steep|R0|virg|NatM|Natural_M|M_|Beta|sel|Sel|Q_|q_|F_|break)",
    pars$Label,
    ignore.case = TRUE
  )
  keep <- keep & !grepl("RecrDev|_Dev_|GP_|Fem|Male|Wtlen|Mat|Eggs|VonBert|CV_|L_at", pars$Label)
  out <- pars[keep, , drop = FALSE]
  if (!NROW(out)) {
    return(NULL)
  }
  cols <- intersect(
    names(out),
    c("Label", "Value", "Min", "Max", "Phase", "Status", "Parm_StDev")
  )
  out <- out[, cols, drop = FALSE]
  names(out)[names(out) == "Parm_StDev"] <- "se"
  out[order(out$Label), , drop = FALSE]
}

recruitment_deviations_table <- function(rep) {
  pars <- rep$parameters
  if (is.null(pars) || !NROW(pars)) {
    return(NULL)
  }
  dev <- pars[grepl("RecrDev", pars$Label), , drop = FALSE]
  if (!NROW(dev)) {
    return(NULL)
  }
  year <- suppressWarnings(as.integer(sub(".*_(\\d{4})$", "\\1", dev$Label)))
  data.frame(
    year = year,
    deviation = dev$Value,
    status = dev$Status,
    stringsAsFactors = FALSE
  )
}

composition_ess_table <- function(rep) {
  comp <- rep$Length_Comp_Fit_Summary
  if (is.null(comp) || !NROW(comp)) {
    return(NULL)
  }
  cols <- intersect(
    names(comp),
    c(
      "Fleet_name", "Fleet", "N", "mean_Nsamp_in", "mean_Nsamp_adj",
      "mean_effN", "HarMean", "Curr_Var_Adj", "Recommend_var_adj"
    )
  )
  out <- comp[, cols, drop = FALSE]
  if ("mean_effN" %in% names(out)) {
    out$mean_effN <- signif(out$mean_effN, 4)
  }
  if ("HarMean" %in% names(out)) {
    out$HarMean <- signif(out$HarMean, 4)
  }
  out
}

index_tuning_table <- function(rep) {
  idx <- rep$index_variance_tuning_check
  if (is.null(idx) || !NROW(idx)) {
    return(NULL)
  }
  cols <- intersect(
    names(idx),
    c("fleetname", "Fleet", "N", "Npos", "RMSE", "Input+VarAdj+extra", "VarAdj", "New_VarAdj")
  )
  out <- idx[, cols, drop = FALSE]
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) signif(x, 4))
  out
}

recruitment_process_table <- function(rep) {
  pieces <- list()
  if (!is.null(rep$sigma_R_info) && NROW(rep$sigma_R_info)) {
    pieces$sigma_R <- rep$sigma_R_info
  }
  if (!is.null(rep$rmse_table) && NROW(rep$rmse_table)) {
    pieces$rmse <- rep$rmse_table
  }
  if (!length(pieces)) {
    return(NULL)
  }
  pieces
}

plot_likelihood_components <- function(rep) {
  tab <- likelihood_components_table(rep)
  if (!require_ok(!is.null(tab) && NROW(tab), "plot_likelihood_components: likelihood table empty")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  tab$component <- factor(tab$component, levels = rev(tab$component))
  ggplot2::ggplot(tab, ggplot2::aes(x = component, y = log_likelihood)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Log-likelihood component") +
    ggplot2::theme_minimal()
}

plot_recruitment_deviations <- function(rep) {
  tab <- recruitment_deviations_table(rep)
  if (!require_ok(!is.null(tab) && NROW(tab), "plot_recruitment_deviations: recruitment deviations missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  tab <- tab[order(tab$year), , drop = FALSE]
  ggplot2::ggplot(tab, ggplot2::aes(x = year, y = deviation)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_point(colour = "steelblue", size = 1.2) +
    ggplot2::labs(x = "Year", y = "Recruitment deviation") +
    ggplot2::theme_minimal()
}

plot_cpue_fit <- function(rep, fleets = NULL, max_fleets = 4L) {
  cpue <- rep$cpue
  if (!require_ok(!is.null(cpue) && NROW(cpue), "plot_cpue_fit: CPUE data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (is.null(fleets)) {
    fleets <- unique(cpue$Fleet)[seq_len(min(max_fleets, length(unique(cpue$Fleet))))]
  }
  dat <- subset(cpue, Fleet %in% fleets & Use %in% c(1, "1", TRUE))
  if (!require_ok(NROW(dat), "plot_cpue_fit: no CPUE rows for selected fleets")) return(NULL)
  dat$series <- dat$Fleet_name
  ggplot2::ggplot(dat, ggplot2::aes(x = Yr)) +
    ggplot2::geom_line(ggplot2::aes(y = Obs, colour = series)) +
    ggplot2::geom_line(ggplot2::aes(y = Exp, colour = series), linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Year",
      y = "Index (log scale)",
      colour = "Fleet",
      caption = "Solid: observed; dashed: expected."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

plot_cpue_residuals <- function(rep, fleets = NULL, max_fleets = 4L) {
  cpue <- rep$cpue
  if (!require_ok(!is.null(cpue) && NROW(cpue), "plot_cpue_residuals: CPUE data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (is.null(fleets)) {
    fleets <- unique(cpue$Fleet)[seq_len(min(max_fleets, length(unique(cpue$Fleet))))]
  }
  dat <- subset(cpue, Fleet %in% fleets & Use %in% c(1, "1", TRUE))
  if (!require_ok(NROW(dat), "plot_cpue_residuals: no CPUE rows for selected fleets")) return(NULL)
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~Fleet_name, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Year", y = "CPUE residual") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

qic_individual_stats <- function(year, residual) {
  ord <- order(year, na.last = NA)
  rr <- as.numeric(residual[ord])
  rr <- rr[is.finite(rr)]
  sgn <- ifelse(rr >= 0, 1L, -1L)
  n <- length(sgn)
  if (n < 2L) {
    return(data.frame(
      n.crossings = NA_real_,
      n.crossings.min = NA_real_,
      longest.run = NA_real_,
      longest.run.max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    n.crossings = sum(diff(sgn) != 0),
    n.crossings.min = floor((n - 1) * 0.25),
    longest.run = max(base::rle(sgn)$lengths),
    longest.run.max = max(2, ceiling(log(n, base = 2) + 2)),
    stringsAsFactors = FALSE
  )
}

cpue_runs_qic_table <- function(diag_df, scenario = "reference") {
  if (is.null(diag_df) || !NROW(diag_df) || !safe_require("plyr")) {
    return(NULL)
  }
  req <- c("year", "residual", "name")
  if (!all(req %in% names(diag_df))) {
    return(NULL)
  }
  dat <- diag_df[is.finite(diag_df$residual), req, drop = FALSE]
  dat$scenario <- scenario
  stats <- plyr::ddply(dat, plyr::.(scenario, name), function(x) {
    qic_individual_stats(x$year, x$residual)
  })
  summary <- plyr::ddply(
    stats,
    plyr::.(scenario, name),
    function(x) {
      data.frame(
        crossings = all(x$n.crossings > x$n.crossings.min, na.rm = TRUE),
        runs = all(x$longest.run < x$longest.run.max, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
  out <- merge(stats, summary, by = c("scenario", "name"), all.x = TRUE)
  out$Pass <- out$crossings & out$runs
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) signif(x, 4))
  out
}

plot_cpue_runs_qic <- function(diag_df, scenario = "reference") {
  if (!require_ok(!is.null(diag_df) && NROW(diag_df), "plot_cpue_runs_qic: diagnostic data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  req <- c("year", "residual", "name")
  if (!require_ok(all(req %in% names(diag_df)), "plot_cpue_runs_qic: required columns missing")) return(NULL)
  rns <- diag_df[is.finite(diag_df$residual), , drop = FALSE]
  rns$scenario <- scenario
  if (!safe_require("plyr")) return(NULL)
  runs <- cpue_runs_qic_table(rns, scenario = scenario)
  if (!require_ok(!is.null(runs) && NROW(runs), "plot_cpue_runs_qic: runs table empty")) return(NULL)
  neg <- subset(rns, residual < 0)
  pos <- subset(rns, residual > 0)
  yr_breaks <- pretty(rns$year, n = 8)

  ggplot2::ggplot(rns, ggplot2::aes(x = year, y = residual)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_line(colour = "grey60", linewidth = 0.5) +
    ggplot2::geom_point(
      data = neg,
      colour = "red",
      size = 1.6
    ) +
    ggplot2::geom_linerange(
      data = neg,
      ggplot2::aes(ymin = 0, ymax = residual),
      colour = "red",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = pos,
      colour = "black",
      size = 1.6
    ) +
    ggplot2::geom_linerange(
      data = pos,
      ggplot2::aes(ymin = 0, ymax = residual),
      colour = "black",
      linewidth = 0.5
    ) +
    ggplot2::geom_rect(
      data = runs,
      ggplot2::aes(
        fill = Pass,
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.25,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "red", `TRUE` = "green"),
      name = "Run test pass"
    ) +
    ggplot2::facet_grid(name ~ scenario, scales = "free", space = "free_x") +
    ggplot2::scale_x_continuous(breaks = yr_breaks) +
    ggplot2::labs(
      x = "Year",
      y = "CPUE residual",
      caption = "Individual-chart runs test (crossings and longest-run criteria; cf. qic(..., chart='i'))."
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.text.x = ggplot2::element_text(angle = 0),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

cpue_runstest_table <- function(rep, mixing = "less") {
  if (!require_ok(requireNamespace("ss3diags", quietly = TRUE), "cpue_runstest_table: package 'ss3diags' required")) return(NULL)
  if (!require_ok(!is.null(rep$cpue), "cpue_runstest_table: CPUE data missing")) return(NULL)
  tbl <- ss3diags::SSrunstest(rep, mixing = mixing, quants = "cpue", verbose = FALSE)
  if (!require_ok(!is.null(tbl) && NROW(tbl), "cpue_runstest_table: SSrunstest returned no rows")) return(NULL)
  num_cols <- vapply(tbl, is.numeric, logical(1))
  tbl[num_cols] <- lapply(tbl[num_cols], function(x) signif(x, 4))
  tbl
}

plot_cpue_runstest <- function(rep, mixing = "less", ncol = 3L) {
  if (!safe_require("ggplot2")) return(NULL)
  if (!require_ok(requireNamespace("ss3diags", quietly = TRUE), "plot_cpue_runstest: package 'ss3diags' required")) return(NULL)
  cpue <- rep$cpue
  if (!require_ok(!is.null(cpue) && NROW(cpue), "plot_cpue_runstest: CPUE data missing")) return(NULL)
  fleet_col <- if ("Fleet_name" %in% names(cpue)) {
    "Fleet_name"
  } else if ("Name" %in% names(cpue)) {
    "Name"
  } else {
    return(NULL)
  }
  cpue$residuals <- ifelse(
    is.na(cpue$Obs) | is.na(cpue$Like),
    NA_real_,
    log(cpue$Obs) - log(cpue$Exp)
  )
  cpue <- cpue[is.finite(cpue$residuals), , drop = FALSE]
  if (!NROW(cpue)) {
    return(NULL)
  }

  fleets <- unique(cpue[[fleet_col]])
  point_rows <- list()
  band_rows <- list()

  for (fl in fleets) {
    resid <- cpue[cpue[[fleet_col]] == fl, , drop = FALSE]
    rr <- as.numeric(resid$residuals)
    if (sum(is.finite(rr)) < 3L) {
      next
    }
    rt <- ss3diags::ssruns_sig3(rr, type = "resid", mixing = mixing)
    lims <- rt$sig3lim
    passed <- !is.na(rt$p.runs) && rt$p.runs >= 0.05
    yr <- resid$Yr
    ti <- if ("Time" %in% names(resid)) resid$Time else yr
    fleet_label <- paste0(
      fl,
      "\n(p = ",
      ifelse(is.na(rt$p.runs), "NA", signif(rt$p.runs, 3)),
      "; ",
      ifelse(passed, "pass", "fail"),
      ")"
    )
    point_rows[[length(point_rows) + 1L]] <- data.frame(
      fleet_label = fleet_label,
      year = yr,
      time = ti,
      residual = rr,
      lo = lims[1],
      hi = lims[2],
      passed = passed,
      stringsAsFactors = FALSE
    )
    band_rows[[length(band_rows) + 1L]] <- data.frame(
      fleet_label = fleet_label,
      xmin = min(yr, na.rm = TRUE) - 1,
      xmax = max(yr, na.rm = TRUE) + 1,
      ymin = lims[1],
      ymax = lims[2],
      passed = passed,
      stringsAsFactors = FALSE
    )
  }

  if (!length(point_rows)) {
    return(NULL)
  }
  pts <- do.call(rbind, point_rows)
  bands <- do.call(rbind, band_rows)
  bands$fill <- ifelse(bands$passed, "Passed", "Failed")
  pts$outside <- pts$residual < pts$lo | pts$residual > pts$hi

  ggplot2::ggplot(pts, ggplot2::aes(x = time, y = residual)) +
    ggplot2::geom_rect(
      data = bands,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill
      ),
      inherit.aes = FALSE,
      alpha = 0.35,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = c(Passed = "#00AA0040", Failed = "#CC000040"),
      name = "Runs test"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_segment(
      ggplot2::aes(x = time, xend = time, y = 0, yend = residual),
      colour = "grey30",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = outside),
      shape = 21,
      fill = "white",
      size = 1.8,
      stroke = 0.35
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#CC0000", `FALSE` = "grey20"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~fleet_label, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      x = "Year",
      y = "Log index residual",
      title = "CPUE residual runs tests",
      caption = "Green/red band = 3-sigma limits; p from ss3diags::ssruns_sig3 / randtests::runs.test."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 8))
}

plot_length_comp_residuals <- function(rep, fleet = 1L) {
  comp <- rep$len_comp_fit_table
  if (!require_ok(!is.null(comp) && NROW(comp), "plot_length_comp_residuals: length comp data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  dat <- subset(comp, Fleet == fleet)
  if (!require_ok(NROW(dat), paste0("plot_length_comp_residuals: no data for fleet ", fleet))) return(NULL)
  fleet_name <- dat$Fleet_Name[1]
  ggplot2::ggplot(dat, ggplot2::aes(x = Yr, y = All_delta)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_point(colour = "steelblue", size = 1.1) +
    ggplot2::labs(
      x = "Year",
      y = "Observed - expected mean length",
      title = paste("Length composition residual:", fleet_name)
    ) +
    ggplot2::theme_minimal()
}

append_figure_tex <- function(params, lines, image_path, caption, label, width = 0.92) {
  if (is.null(image_path) || !file.exists(image_path)) {
    require_ok(
      FALSE,
      paste0("Figure file missing: ", caption),
      cfg = params
    )
    return(lines)
  }
  c(
    lines,
    figure_fragment(params, image_path, caption, label, width = width),
    ""
  )
}

save_base_plot <- function(params, plot_fun, strand, filename,
                           width = 8, height = 6, dpi = 200) {
  fig_dir <- file.path(params$appendix_dir, "figs", strand)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(fig_dir, filename)
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (diag_strict(params)) {
    plot_fun()
  } else {
    try(plot_fun(), silent = TRUE)
  }
  if (file.exists(path) && file.info(path)$size > 1000L) {
    path
  } else {
    require_ok(FALSE, paste0("Base plot not produced: ", file.path(strand, filename)), cfg = params)
    NULL
  }
}

scenario_labels <- function(params, scenarios = NULL) {
  scenarios <- scenarios %||% params$scenarios
  diags::scenarioLabels(scenarios, scenario_register = params$scenario_register)
}

scenario_display_label <- function(id, params = NULL) {
  reg <- if (!is.null(params)) params$scenario_register else NULL
  diags::scenarioDisplayLabel(id, scenario_register = reg)
}

relabel_pe_ids <- function(pe, params) {
  diags::relabelPeIds(pe, scenario_register = params$scenario_register)
}

load_ss3_pe <- function(params, paths = NULL) {
  backend_load_pe(params, paths = paths, for_plots = TRUE)
}

load_ss3_fls <- function(params, paths = NULL) {
  backend_load_fls(params, paths = paths)
}

load_sstruth_internal <- function(path, cfg = NULL) {
  if (is.null(cfg)) {
    stop("load_sstruth_internal() requires cfg for backend dispatch.", call. = FALSE)
  }
  if (!require_ok(requireNamespace("sstruth", quietly = TRUE), "Package 'sstruth' is required.", cfg)) {
    return(NULL)
  }
  if (!require_ok(dir.exists(path), paste0("SS3 scenario path not found: ", path), cfg)) {
    return(NULL)
  }
  safe_require("ggplot2")
  backend_load_internal(cfg, path)
}

ssU_compat <- function(path, cfg = NULL) {
  if (!is.null(cfg)) {
    out <- backend_load_internal(cfg, path)
    return(out$ssu)
  }
  ss3_require_r4ss <- function() safe_require("r4ss")
  if (!ss3_require_r4ss()) {
    stop("Package 'r4ss' is required for CPUE ssU summaries.", call. = FALSE)
  }
  safe_require("plyr")
  if (!require_ok(requireNamespace("diags", quietly = TRUE), "Package 'diags' is required for ssU summaries.", cfg = NULL)) {
    return(NULL)
  }
  diags_fn <- utils::getFromNamespace("diagsFn", "diags")
  cpue <- read_ss3_output(path)$cpue
  res <- subset(cpue, !is.na(Dev), select = c(Fleet_name, Yr, Seas, Obs, Exp, Dev, SE))
  names(res) <- c("name", "year", "season", "obs", "hat", "residual", "se")
  subset(plyr::ddply(res, plyr::.(name, season), diags_fn), !is.na(residual))
}

ssu_summary_table <- function(ssu) {
  if (is.null(ssu) || !is.data.frame(ssu) || !NROW(ssu)) {
    return(NULL)
  }
  cols <- intersect(
    names(ssu),
    c("name", "year", "season", "obs", "hat", "residual", "se", "std_residual", "runs.p", "acf")
  )
  out <- ssu[, cols, drop = FALSE]
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) signif(x, 4))
  out
}

ssu_fleet_summary <- function(ssu) {
  if (is.null(ssu) || !is.data.frame(ssu) || !NROW(ssu)) {
    return(NULL)
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    return(NULL)
  }
  plyr::ddply(
    ssu,
    plyr::.(name),
    function(x) {
      data.frame(
        n = NROW(x),
        mean_abs_residual = signif(mean(abs(x$residual), na.rm = TRUE), 4),
        max_abs_std = if ("std_residual" %in% names(x)) {
          signif(max(abs(x$std_residual), na.rm = TRUE), 4)
        } else {
          NA_real_
        },
        runs_p = if ("runs.p" %in% names(x)) {
          signif(min(x$runs.p, na.rm = TRUE), 4)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )
}

diag_fleet_summary <- function(diag) {
  if (is.null(diag) || !is.data.frame(diag) || !NROW(diag)) {
    return(NULL)
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    return(NULL)
  }
  plyr::ddply(
    diag,
    plyr::.(name),
    function(x) {
      data.frame(
        n = NROW(x),
        mean_abs_residual = signif(mean(abs(x$residual), na.rm = TRUE), 4),
        max_abs_std = if ("std_residual" %in% names(x)) {
          signif(max(abs(x$std_residual), na.rm = TRUE), 4)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )
}

plot_index_fit <- function(diag_df) {
  if (!require_ok(!is.null(diag_df) && NROW(diag_df), "plot_index_fit: diagnostic data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  ggplot2::ggplot(diag_df, ggplot2::aes(x = year)) +
    ggplot2::geom_line(ggplot2::aes(y = obs, colour = name), linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = hat, colour = name), linewidth = 0.6, linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Year",
      y = "Index (log scale)",
      colour = "Fleet",
      caption = "Solid: observed; dashed: expected (sstruth::diagSS)."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

plot_length_comp_from_sslen <- function(ss_len, fleet = 1L) {
  fit <- ss_len$fit
  if (!require_ok(!is.null(fit) && NROW(fit), "plot_length_comp_from_sslen: ssLen fit missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  dat <- fit[fit$Fleet == fleet, , drop = FALSE]
  if (!require_ok(NROW(dat), paste0("plot_length_comp_from_sslen: no data for fleet ", fleet))) return(NULL)
  fleet_name <- dat$Fleet_Name[1]
  delta_col <- if ("All_delta" %in% names(dat)) {
    "All_delta"
  } else if ("F_delta" %in% names(dat)) {
    "F_delta"
  } else {
    NULL
  }
  if (!require_ok(!is.null(delta_col), "plot_length_comp_from_sslen: no residual column in ssLen fit")) return(NULL)
  ggplot2::ggplot(dat, ggplot2::aes(x = Yr, y = .data[[delta_col]])) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_point(colour = "steelblue", size = 1.1) +
    ggplot2::labs(
      x = "Year",
      y = "Observed - expected mean length",
      title = paste("Length composition residual:", fleet_name),
      caption = "From sstruth::ssLen()."
    ) +
    ggplot2::theme_minimal()
}

load_ss3_case_study <- function(params) {
  backend_load_case_study(params)
}

load_assessment_robustness <- function(cfg) {
  retro <- load_collected_retrospectives(cfg)
  if (is.null(retro)) {
    backend_load_robustness(cfg, retro = NULL)
  } else {
    backend_load_robustness(cfg, retro = retro)
  }
}

save_appendix_figure <- function(params, plot, strand, filename,
                                 width = 8, height = 6, dpi = 200) {
  if (is.null(plot)) {
    require_ok(
      FALSE,
      paste0("Figure not produced: ", file.path(strand, filename)),
      cfg = params
    )
    return(invisible(NULL))
  }
  fig_dir <- file.path(params$appendix_dir, "figs", strand)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(fig_dir, filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  invisible(path)
}

summarise_pe_scenarios <- function(pe) {
  if (is.null(pe) || is.null(pe$tseries) || is.null(pe$refpts)) {
    return(NULL)
  }
  ts <- merge(pe$tseries, pe$refpts, by = "id", all.x = TRUE)
  last <- stats::aggregate(
    year ~ id,
    data = ts,
    FUN = function(x) max(x, na.rm = TRUE)
  )
  ts <- merge(ts, last, by = "id", suffixes = c("", ".last"))
  ts <- subset(ts, year == year.last)
  data.frame(
    scenario = ts$id,
    terminal_year = ts$year,
    ssb_bmsy = signif(ts$ssb / ts$bmsy, 3),
    catch_msy = signif(ts$yield / ts$msy, 3),
    harvest_fmsy = signif((ts$yield / ts$ssb) / (ts$msy / ts$bmsy), 3),
    stringsAsFactors = FALSE
  )
}

plot_pe_ts <- function(pe) {
  if (!require_ok(!is.null(pe), "plot_pe_ts: process-error object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("reshape2")) return(NULL)
  ts <- transform(
    merge(pe$tseries, pe$refpts, by = "id"),
    ssb = ssb / bmsy,
    catch = yield / msy,
    harvest = yield / ssb / (msy / bmsy)
  )
  ts <- reshape2::melt(
    ts[, c("id", "year", "ssb", "catch", "harvest")],
    id.vars = c("id", "year")
  )
  ts$variable <- factor(
    ts$variable,
    levels = c("ssb", "catch", "harvest"),
    labels = c("SSB / B_MSY", "Catch / MSY", "Harvest (rel. MSY)")
  )
  ggplot2::ggplot(ts) +
    ggplot2::geom_line(ggplot2::aes(year, value, colour = id), linewidth = 0.35) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      colour = "darkgreen",
      linewidth = 0.35
    ) +
    ggplot2::facet_grid(variable ~ ., scales = "free_y") +
    ggplot2::scale_color_discrete(name = "SS3 scenario") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlab("Year") +
    ggplot2::ylab(NULL)
}

plot_selmat_fls <- function(fls, labels = NULL) {
  if (!require_ok(!is.null(fls), "plot_selmat_fls: FLStocks object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("plyr")) return(NULL)
  if (!requireNamespace("FLCandy", quietly = TRUE)) {
    return(NULL)
  }
  safe_require("dplyr")
  dat <- plyr::ldply(fls, function(x) FLCandy:::matSel(x))
  if (!is.null(labels)) {
    dat$.id <- labels[dat$.id]
  } else {
    dat$.id <- gsub("Scenario-", "Scenario ", dat$.id)
  }
  ggplot2::ggplot(dat, ggplot2::aes(x = age, y = sel, colour = .id)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = mat), linetype = 2, linewidth = 1) +
    ggplot2::labs(
      x = "Age",
      y = "Proportion",
      colour = "Scenario",
      caption = "Solid: selectivity; dashed: maturity."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

curve_ss_production_col <- function(ts) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::curveSsProductionCol(ts))
  }
  for (col in c("sp", "P_ssb", "pf")) {
    if (col %in% names(ts)) {
      return(col)
    }
  }
  NULL
}

curve_ss_terminal_rows <- function(ts) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::curveSsTerminalRows(ts))
  }
  if (is.null(ts) || !NROW(ts) || !"id" %in% names(ts)) {
    return(NULL)
  }
  do.call(rbind, lapply(split(ts, ts$id, drop = TRUE), function(x) {
    x[which.max(x$year), , drop = FALSE]
  }))
}

pe_residual_tseries <- function(pe, pe_mode = NULL) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    mode <- pe_mode %||% Sys.getenv("DIAG_PE_MODE", unset = "relative")
    ts <- sstruth::peResidualTseries(pe, pe_mode = mode)
    ts$pe2 <- ts$pe_resid
    return(ts)
  }
  if (is.null(pe) || is.null(pe$tseries) || !NROW(pe$tseries)) {
    return(NULL)
  }
  ts <- pe$tseries
  if (!"pe2" %in% names(ts) && "pe" %in% names(ts)) {
    ts$pe2 <- ts$pe
  }
  if (!"sp" %in% names(ts)) {
    prod_col <- curve_ss_production_col(ts)
    if (!is.null(prod_col)) {
      ts$sp <- ts[[prod_col]]
    }
  }
  if (!"pe2" %in% names(ts)) {
    return(NULL)
  }
  ts$pe_resid <- ts$pe2
  ts
}

.pe_plot_limits <- function(x, ylim = NULL) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::pePlotLimits(x, ylim = ylim))
  }
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(-1, 1))
  }
  range(x, na.rm = TRUE)
}

plot_pe_residual_time <- function(pe, pe_mode = NULL, ylim = NULL) {
  if (!require_ok(!is.null(pe), "plot_pe_residual_time: process-error object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("plyr")) return(NULL)
  ts <- pe_residual_tseries(pe, pe_mode = pe_mode)
  if (is.null(ts)) {
    return(NULL)
  }
  resid <- ts$pe_resid %||% ts$pe2
  ylim <- .pe_plot_limits(resid, ylim = ylim)
  ts_f <- ts[is.finite(resid) & abs(resid) < 1e4, , drop = FALSE]
  star <- NULL
  if (requireNamespace("FLRebuild", quietly = TRUE) && NROW(ts_f)) {
    star <- plyr::ddply(
      ts_f,
      plyr::.(id),
      with,
      FLRebuild:::rodFn(data = pe_resid[!is.na(pe_resid)], year = year[!is.na(pe_resid)])
    )
  }
  ylab <- if (requireNamespace("sstruth", quietly = TRUE)) {
    sstruth::processErrorLabel(pe_mode %||% Sys.getenv("DIAG_PE_MODE", unset = "relative"))
  } else {
    "Process residual"
  }
  p <- ggplot2::ggplot(ts, ggplot2::aes(x = year, y = pe_resid)) +
    ggplot2::facet_grid(. ~ id, scales = "free") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::labs(x = "Year", y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(ylim = ylim)
  if (!is.null(star) && NROW(star)) {
    p <- p +
      ggplot2::geom_polygon(
        ggplot2::aes(year, data, group = regime),
        fill = "blue",
        alpha = 0.2,
        data = star
      )
  }
  p
}

plot_pe_residual_acf <- function(pe, pe_mode = NULL) {
  if (!require_ok(!is.null(pe), "plot_pe_residual_acf: process-error object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("plyr")) return(NULL)
  ts <- pe_residual_tseries(pe, pe_mode = pe_mode)
  if (is.null(ts)) {
    return(NULL)
  }
  resid <- ts$pe_resid %||% ts$pe2
  ts_f <- ts[is.finite(resid) & abs(resid) < 1e4, , drop = FALSE]
  acf_df <- plyr::ddply(
    ts_f,
    plyr::.(id),
    with,
    cbind(
      acf = stats::acf(pe_resid[!is.na(pe_resid)], lag.max = 25, plot = FALSE)$acf,
      lag = 0:25
    )
  )
  crit <- stats::qnorm(0.975) / sqrt(nrow(pe$tseries))
  ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::facet_grid(. ~ id, scales = "free") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_hline(
      yintercept = c(-crit, crit),
      linetype = "dashed",
      colour = "grey60"
    ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
}

plot_pe_phase <- function(pe, ylab = "Production") {
  if (!require_ok(!is.null(pe), "plot_pe_phase: process-error object missing")) return(NULL)
  if (!require_ok(!is.null(pe$tseries) && !is.null(pe$curve) && !is.null(pe$refpts),
                  "plot_pe_phase: curveSS tseries/curve/refpts required")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  ts <- pe_residual_tseries(pe)
  if (is.null(ts) || !"sp" %in% names(ts)) {
    return(NULL)
  }
  terminal <- curve_ss_terminal_rows(ts)
  ggplot2::ggplot(ts, ggplot2::aes(x = ssb, y = sp)) +
    ggplot2::facet_grid(. ~ id) +
    ggplot2::geom_path(ggplot2::aes(x = ssb, y = sp), colour = "grey", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(fill = year), alpha = 0.35, size = 1.6, shape = 21) +
    ggplot2::geom_line(
      data = pe$curve,
      ggplot2::aes(x = ssb, y = yield, group = id),
      colour = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = pe$refpts,
      ggplot2::aes(x = bmsy, y = msy),
      colour = "darkgreen",
      size = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(
      data = terminal,
      ggplot2::aes(x = ssb, y = sp),
      colour = "blue",
      size = 2.5
    ) +
    ggplot2::geom_text(
      data = terminal,
      ggplot2::aes(x = ssb * 0.82, y = sp, label = year),
      colour = "blue",
      size = 3
    ) +
    ggplot2::labs(x = "SSB", y = ylab, caption = "Grey: empirical production path; black: curveSS yield curve; green: MSY.") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

plot_pe_yield_phase <- function(pe) {
  if (!require_ok(!is.null(pe), "plot_pe_yield_phase: process-error object missing")) return(NULL)
  if (!require_ok(!is.null(pe$tseries) && !is.null(pe$curve) && !is.null(pe$refpts),
                  "plot_pe_yield_phase: curveSS tseries/curve/refpts required")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  ts <- pe$tseries
  terminal <- curve_ss_terminal_rows(ts)
  ggplot2::ggplot(ts, ggplot2::aes(x = ssb, y = yield)) +
    ggplot2::facet_grid(. ~ id) +
    ggplot2::geom_path(colour = "grey", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(fill = year), alpha = 0.35, size = 1.6, shape = 21) +
    ggplot2::geom_line(
      data = pe$curve,
      ggplot2::aes(x = ssb, y = yield, group = id),
      colour = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = pe$refpts,
      ggplot2::aes(x = bmsy, y = msy),
      colour = "darkgreen",
      size = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(data = terminal, colour = "blue", size = 2.5) +
    ggplot2::labs(x = "SSB", y = "Yield", caption = "From curveSS().") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

production_facet_cols <- function(pe) {
  if (!is.null(pe$refpts) && all(c("M", "beta") %in% names(pe$refpts))) {
    return(c("M", "beta"))
  }
  "id"
}

production_facet_limits <- function(crv, rfs, facet_cols, traj = NULL) {
  if (!requireNamespace("plyr", quietly = TRUE)) {
    return(rfs)
  }
  lim <- plyr::ddply(crv, facet_cols, function(df) {
    data.frame(
      b0 = max(df$ssb, na.rm = TRUE),
      curve_ymax = max(df$yield, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  rfs2 <- merge(rfs, lim, by = facet_cols, all.x = TRUE, sort = FALSE)
  if (!is.null(traj) && NROW(traj)) {
    traj_lim <- plyr::ddply(traj, facet_cols, function(df) {
      data.frame(
        traj_xmax = max(df$ssb, na.rm = TRUE),
        traj_ymax = max(df$yield, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    rfs2 <- merge(rfs2, traj_lim, by = facet_cols, all.x = TRUE, sort = FALSE)
  } else {
    rfs2$traj_xmax <- NA_real_
    rfs2$traj_ymax <- NA_real_
  }
  rfs2$x_cap <- pmax(rfs2$b0, rfs2$traj_xmax, rfs2$bmsy, na.rm = TRUE)
  rfs2$y_cap <- pmax(rfs2$curve_ymax, rfs2$msy, rfs2$traj_ymax, na.rm = TRUE)
  rfs2
}

production_panel_limits_table <- function(rfs, facet_cols, xlim = NULL, ylim = NULL) {
  cols <- unique(c(facet_cols, "x_cap", "y_cap"))
  cols <- intersect(cols, names(rfs))
  caps <- unique(rfs[, cols, drop = FALSE])
  caps$x_min <- if (!is.null(xlim) && length(xlim) == 2L) xlim[1L] else 0
  caps$x_max <- if (!is.null(xlim) && length(xlim) == 2L) {
    xlim[2L]
  } else {
    caps$x_cap
  }
  caps$y_min <- if (!is.null(ylim) && length(ylim) == 2L) ylim[1L] else 0
  caps$y_max <- if (!is.null(ylim) && length(ylim) == 2L) {
    ylim[2L]
  } else {
    caps$y_cap
  }
  caps
}

production_facet_formula_lhs <- function(facet_cols, row) {
  parts <- vapply(facet_cols, function(col) {
    val <- row[[col]]
    if (is.character(val) || is.factor(val)) {
      paste0(col, ' == "', gsub('"', '\\\\"', as.character(val)), '"')
    } else {
      paste0(col, " == ", val)
    }
  }, character(1))
  paste(parts, collapse = " & ")
}

production_facetted_scale_formulas <- function(caps, facet_cols, axis = c("x", "y")) {
  axis <- match.arg(axis)
  lo_col <- if (axis == "x") "x_min" else "y_min"
  hi_col <- if (axis == "x") "x_max" else "y_max"
  scale_fn <- if (axis == "x") "scale_x_continuous" else "scale_y_continuous"
  lapply(seq_len(nrow(caps)), function(i) {
    row <- caps[i, , drop = FALSE]
    lhs <- production_facet_formula_lhs(facet_cols, row)
    as.formula(paste0(
      lhs,
      " ~ ",
      scale_fn,
      "(limits = c(",
      row[[lo_col]], ", ", row[[hi_col]],
      "), expand = c(0, 0), labels = scales::comma)"
    ))
  })
}

production_facet_scales_arg <- function(facet_cols, auto_x, auto_y) {
  if (auto_x && auto_y) {
    "free"
  } else if (auto_x) {
    "free_x"
  } else if (auto_y) {
    "free_y"
  } else {
    "fixed"
  }
}

production_add_axis_limits <- function(p, rfs, facet_cols, xlim = NULL, ylim = NULL) {
  global_x <- !is.null(xlim) && length(xlim) == 2L
  global_y <- !is.null(ylim) && length(ylim) == 2L
  auto_x <- !global_x
  auto_y <- !global_y
  use_fps <- (auto_x || auto_y) && requireNamespace("ggh4x", quietly = TRUE)

  if (use_fps) {
    caps <- production_panel_limits_table(rfs, facet_cols, xlim = xlim, ylim = ylim)
    if (global_y) {
      p <- p + ggplot2::scale_y_continuous(
        limits = ylim,
        expand = ggplot2::expansion(c(0, 0)),
        labels = scales::comma
      )
    }
    if (global_x) {
      p <- p + ggplot2::scale_x_continuous(
        limits = xlim,
        expand = ggplot2::expansion(c(0, 0)),
        labels = scales::comma
      )
    }
    fps_x <- if (auto_x) production_facetted_scale_formulas(caps, facet_cols, "x") else NULL
    fps_y <- if (auto_y) production_facetted_scale_formulas(caps, facet_cols, "y") else NULL
    return(p + ggh4x::facetted_pos_scales(x = fps_x, y = fps_y))
  }

  if (global_x || global_y) {
    coord_args <- list(expand = FALSE)
    if (global_x) {
      coord_args$xlim <- xlim
    }
    if (global_y) {
      coord_args$ylim <- ylim
    }
    return(p + do.call(ggplot2::coord_cartesian, coord_args))
  }

  p + ggplot2::coord_cartesian(expand = TRUE)
}

production_clip_triangle <- function(trngl, rfs, facet_cols) {
  if (is.null(trngl) || !NROW(trngl) || !all(c("x", "y") %in% names(trngl))) {
    return(trngl)
  }
  cap_cols <- unique(c(facet_cols, "x_cap", "y_cap"))
  caps <- unique(rfs[, cap_cols, drop = FALSE])
  out <- merge(trngl, caps, by = facet_cols, all.x = TRUE, sort = FALSE)
  out$x <- pmin(out$x, out$x_cap)
  out$y <- pmin(out$y, out$y_cap)
  out$x_cap <- NULL
  out$y_cap <- NULL
  out
}

production_curve_yield_at <- function(df, ssb_target) {
  ok <- is.finite(df$ssb) & is.finite(df$yield)
  df <- df[ok, , drop = FALSE]
  if (NROW(df) < 2L) {
    return(NA_real_)
  }
  df <- df[order(df$ssb), , drop = FALSE]
  if (anyDuplicated(df$ssb)) {
    df <- stats::aggregate(yield ~ ssb, df, max)
  }
  stats::approx(df$ssb, df$yield, xout = ssb_target, rule = 2)$y
}

production_b40_path <- function(crv, by) {
  if (!requireNamespace("plyr", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(crv) || !NROW(crv)) {
    return(NULL)
  }
  out <- plyr::ddply(crv, by, function(df) {
    b0 <- max(df$ssb, na.rm = TRUE)
    if (!is.finite(b0) || b0 <= 0) {
      return(NULL)
    }
    target <- b0 * 0.4
    y40 <- production_curve_yield_at(df, target)
    if (!is.finite(y40)) {
      dif <- abs(df$ssb - target)
      idx <- order(dif)[1L]
      target <- df$ssb[idx]
      y40 <- df$yield[idx]
    }
    data.frame(
      ssb = c(target, target, b0),
      yield = c(0, y40, y40),
      .path_order = 1:3,
      stringsAsFactors = FALSE
    )
  })
  if (is.null(out) || !NROW(out)) {
    return(NULL)
  }
  grp_cols <- intersect(by, names(out))
  if (length(grp_cols)) {
    out$.b40_group <- interaction(out[grp_cols], drop = TRUE)
  } else {
    out$.b40_group <- seq_len(nrow(out))
  }
  out[order(out$.b40_group, out$.path_order), , drop = FALSE]
}

subset_production_pe <- function(pe, M_filter = NULL, beta_filter = NULL, id_filter = NULL) {
  if (is.null(pe)) {
    return(NULL)
  }
  subset_df <- function(df) {
    if (is.null(df) || !is.data.frame(df) || !NROW(df)) {
      return(df)
    }
    out <- df
    if (!is.null(M_filter) && "M" %in% names(out)) {
      out <- out[out$M %in% M_filter, , drop = FALSE]
    }
    if (!is.null(beta_filter) && "beta" %in% names(out)) {
      out <- out[out$beta %in% beta_filter, , drop = FALSE]
    }
    if (!is.null(id_filter) && "id" %in% names(out)) {
      out <- out[out$id %in% id_filter, , drop = FALSE]
    }
    out
  }
  list(
    tseries = subset_df(pe$tseries),
    curve = subset_df(pe$curve),
    refpts = subset_df(pe$refpts),
    triangle = subset_df(pe$triangle),
    derived = if (!is.null(pe$derived)) subset_df(pe$derived) else NULL
  )
}

load_profile_pe <- function(cfg, prof = NULL, maxY = 2.5, parallel = TRUE, workers = NULL) {
  prof <- prof %||% collect_likelihood_profiles(cfg)
  if (is.null(prof) || is.null(prof$manifest) || !NROW(prof$manifest)) {
    return(NULL)
  }
  mf <- prof$manifest
  mf <- mf[dir.exists(mf$dir), , drop = FALSE]
  if (!NROW(mf)) {
    return(NULL)
  }
  if (!requireNamespace("FLRebuild", quietly = TRUE)) {
    stop("Package 'FLRebuild' is required for profile production functions.", call. = FALSE)
  }
  if (requireNamespace("magrittr", quietly = TRUE)) {
    suppressPackageStartupMessages(library(magrittr, quietly = TRUE))
  }
  if (!exists("parallel_lapply", mode = "function", inherits = TRUE)) {
    longrun_path <- file.path(cfg$appendix_dir, "R", "longrun-utils.R")
    if (file.exists(longrun_path)) {
      source(longrun_path, local = FALSE)
    }
  }

  run_one <- function(i) {
    out <- tryCatch(
      FLRebuild::curveSS(mf$dir[i], maxY = maxY),
      error = function(e) NULL
    )
    if (is.null(out)) {
      return(NULL)
    }
    lapply(out, function(df) {
      if (is.data.frame(df)) {
        cbind(M = mf$M[i], beta = mf$beta[i], df, stringsAsFactors = FALSE)
      } else {
        df
      }
    })
  }
  runs <- parallel_lapply(
    seq_len(nrow(mf)),
    run_one,
    parallel = parallel,
    workers = workers
  )
  runs <- runs[!vapply(runs, is.null, logical(1))]
  if (!length(runs)) {
    return(NULL)
  }
  parts <- names(runs[[1L]])
  combined <- lapply(parts, function(part) {
    do.call(rbind, lapply(runs, `[[`, part))
  })
  names(combined) <- parts
  combined
}

production_trajectory_df <- function(pe) {
  ts <- pe$tseries
  if (is.null(ts) || !is.data.frame(ts) || !NROW(ts)) {
    return(NULL)
  }
  req <- c("year", "ssb", "yield")
  if (!all(req %in% names(ts))) {
    return(NULL)
  }
  keep <- unique(c(intersect(c("M", "beta", "id"), names(ts)), req))
  ts <- ts[is.finite(ts$ssb) & is.finite(ts$yield) & is.finite(ts$year), keep, drop = FALSE]
  if (!"id" %in% names(ts)) {
    ts$id <- "."
  }
  if (!NROW(ts)) {
    return(NULL)
  }
  ord <- intersect(c("M", "beta", "id", "year"), names(ts))
  ts[do.call(order, ts[ord]), , drop = FALSE]
}

production_terminal_rows <- function(traj, facet_cols = "id") {
  if (is.null(traj) || !NROW(traj) || !"year" %in% names(traj)) {
    return(NULL)
  }
  grp <- intersect(c(facet_cols, "id"), names(traj))
  grp <- setdiff(grp, ".")
  if (!length(grp)) {
    grp <- intersect(c("M", "beta"), names(traj))
  }
  if (!length(grp)) {
    return(traj[which.max(traj$year), , drop = FALSE])
  }
  if (length(grp) == 1L) {
    split_key <- traj[[grp[1L]]]
  } else {
    split_key <- interaction(traj[grp], drop = TRUE)
  }
  do.call(rbind, lapply(split(traj, split_key, drop = TRUE), function(x) {
    x[which.max(x$year), , drop = FALSE]
  }))
}

production_trajectory_group <- function(traj, facet_cols) {
  grp <- intersect(c(facet_cols, "id"), names(traj))
  grp <- setdiff(grp, ".")
  if (length(grp) <= 1L) {
    if (length(grp) == 1L) {
      return(traj[[grp[1L]]])
    }
    return(seq_len(nrow(traj)))
  }
  interaction(traj[grp], drop = TRUE)
}

plot_production_functions <- function(
  pe,
  facet_cols = NULL,
  M_filter = NULL,
  beta_filter = NULL,
  id_filter = NULL,
  xlim = NULL,
  ylim = NULL,
  title = "Production Function",
  show_trajectories = NULL,
  trajectory_yield_scale = 1
) {
  if (!require_ok(!is.null(pe), "plot_production_functions: curveSS object missing")) {
    return(NULL)
  }
  if (!require_ok(
    !is.null(pe$refpts) && !is.null(pe$curve) && NROW(pe$refpts) && NROW(pe$curve),
    "plot_production_functions: curveSS refpts/curve required"
  )) {
    return(NULL)
  }
  if (!safe_require("ggplot2")) {
    return(NULL)
  }
  if (!safe_require("scales")) {
    return(NULL)
  }

  pe <- subset_production_pe(pe, M_filter = M_filter, beta_filter = beta_filter, id_filter = id_filter)
  facet_cols <- facet_cols %||% production_facet_cols(pe)
  rfs <- pe$refpts
  crv <- pe$curve
  trngl <- pe$triangle
  if (!NROW(rfs) || !NROW(crv)) {
    return(NULL)
  }

  if (is.null(show_trajectories)) {
    show_trajectories <- length(facet_cols) == 1L && identical(facet_cols, "id") ||
      length(facet_cols) == 2L && identical(facet_cols, c("M", "beta"))
  }
  traj <- if (isTRUE(show_trajectories)) production_trajectory_df(pe) else NULL
  if (!is.null(traj) && NROW(traj) && is.finite(trajectory_yield_scale) &&
      !identical(trajectory_yield_scale, 1)) {
    traj$yield <- traj$yield * trajectory_yield_scale
  }
  terminal <- if (!is.null(traj)) production_terminal_rows(traj, facet_cols) else NULL

  rfs <- production_facet_limits(crv, rfs, facet_cols, traj = traj)
  global_x <- !is.null(xlim) && length(xlim) == 2L
  global_y <- !is.null(ylim) && length(ylim) == 2L
  auto_x <- !global_x
  auto_y <- !global_y
  use_fps <- (auto_x || auto_y) && requireNamespace("ggh4x", quietly = TRUE)
  facet_scales <- production_facet_scales_arg(facet_cols, auto_x, auto_y)

  b40 <- production_b40_path(crv, facet_cols)
  if (is.null(b40)) {
    return(NULL)
  }

  if (length(facet_cols) == 2L && identical(facet_cols, c("M", "beta"))) {
    facet <- ggplot2::facet_grid(
      beta ~ M,
      labeller = ggplot2::labeller(
        M = function(x) paste("M *", x),
        beta = function(x) paste("beta =", x)
      ),
      scales = facet_scales,
      space = "free"
    )
  } else if (length(facet_cols) == 1L && facet_cols == "id") {
    facet <- ggplot2::facet_wrap(
      ~id,
      scales = facet_scales,
      ncol = min(3L, length(unique(rfs$id)))
    )
  } else {
    return(NULL)
  }

  if (!use_fps && !is.null(trngl) && NROW(trngl)) {
    trngl <- production_clip_triangle(trngl, rfs, facet_cols)
  }

  p <- ggplot2::ggplot(rfs) + facet
  if (!is.null(trngl) && NROW(trngl) && all(c("x", "y") %in% names(trngl))) {
    p <- p + ggplot2::geom_polygon(
      data = trngl,
      ggplot2::aes(x = x, y = y),
      fill = "orange",
      alpha = 0.5
    )
  }
  p <- p +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 0, xmax = bmsy, ymin = msy, ymax = Inf),
      fill = "#D55E00",
      alpha = 0.3
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = bmsy, xmax = Inf, ymin = 0, ymax = msy),
      fill = "#009E73",
      alpha = 0.3
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 0, xmax = bmsy, ymin = 0, ymax = msy),
      fill = "#F0E442",
      alpha = 0.3
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = bmsy, xmax = Inf, ymin = msy, ymax = Inf),
      fill = "#F0E442",
      alpha = 0.3
    )

  if (!use_fps) {
    p <- p +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::scale_y_continuous(labels = scales::comma)
  }

  if (!is.null(traj) && NROW(traj)) {
    traj$.traj_group <- production_trajectory_group(traj, facet_cols)
    p <- p +
      ggplot2::geom_path(
        data = traj,
        ggplot2::aes(x = ssb, y = yield, group = .traj_group),
        colour = "grey45",
        linewidth = 0.45,
        alpha = 0.85
      ) +
      ggplot2::geom_point(
        data = traj,
        ggplot2::aes(x = ssb, y = yield, colour = year),
        size = 1.8,
        alpha = 0.55
      ) +
      ggplot2::scale_color_viridis_c(name = "Year", option = "C")
  }

  p <- p +
    ggplot2::geom_line(
      data = crv,
      ggplot2::aes(x = ssb, y = yield),
      colour = "blue",
      linewidth = 0.8
    ) +
    ggplot2::geom_abline(
      ggplot2::aes(slope = msy / bmsy, intercept = 0),
      colour = "grey10",
      linetype = 3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = bmsy, y = msy),
      shape = 21,
      fill = "white",
      size = 3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = bmsy, y = msy, label = "MSY"),
      hjust = -0.2,
      vjust = -0.2,
      size = 3.5
    ) +
    ggplot2::geom_path(
      data = b40,
      ggplot2::aes(x = ssb, y = yield, group = .b40_group),
      colour = "green",
      linewidth = 0.8
    )

  if (!is.null(terminal) && NROW(terminal)) {
    p <- p +
      ggplot2::geom_point(
        data = terminal,
        ggplot2::aes(x = ssb, y = yield),
        colour = "blue",
        size = 2.8,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = terminal,
        ggplot2::aes(x = ssb * 0.82, y = yield, label = year),
        colour = "blue",
        size = 3,
        inherit.aes = FALSE
      )
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = "Spawning stock biomass (t)",
      y = "Production / yield (t)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9, margin = ggplot2::margin(t = 8)),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
      axis.line = ggplot2::element_line(colour = "black")
    )

  production_add_axis_limits(p, rfs, facet_cols, xlim = xlim, ylim = ylim)
}

#' Build Kobe phase-plot trajectories (B/BMSY vs F/FMSY) from a curveSS object
kobe_trajectory_df <- function(pe) {
  if (is.null(pe) || is.null(pe$tseries) || is.null(pe$refpts)) {
    return(NULL)
  }
  ts <- pe$tseries
  rf <- pe$refpts
  by_cols <- intersect(c("id", "M", "beta"), intersect(names(ts), names(rf)))
  if (!length(by_cols)) {
    by_cols <- intersect("id", c(names(ts), names(rf)))
  }
  if (!length(by_cols)) {
    return(NULL)
  }
  ts <- merge(ts, rf, by = by_cols, suffixes = c("", ".ref"))
  if ("bmsy.ref" %in% names(ts) && !"bmsy" %in% names(ts)) {
    ts$bmsy <- ts$bmsy.ref
  }
  if ("msy.ref" %in% names(ts) && !"msy" %in% names(ts)) {
    ts$msy <- ts$msy.ref
  }
  if (!all(c("ssb", "bmsy") %in% names(ts))) {
    return(NULL)
  }
  ts$stock <- ts$ssb / ts$bmsy
  if ("ffmsy" %in% names(ts)) {
    ts$harvest <- ts$ffmsy
  } else if ("bratio" %in% names(ts) && all(c("F", "fmsy") %in% names(ts))) {
    ts$harvest <- ts$F / ts$fmsy
  } else if (all(c("yield", "msy") %in% names(ts))) {
    ts$harvest <- ts$yield / ts$ssb / (ts$msy / ts$bmsy)
  } else {
    return(NULL)
  }
  keep <- unique(c(intersect(c("M", "beta", "id"), names(ts)), "year", "stock", "harvest"))
  ts <- ts[
    is.finite(ts$stock) & is.finite(ts$harvest) & is.finite(ts$year),
    keep,
    drop = FALSE
  ]
  if (!"id" %in% names(ts)) {
    ts$id <- "."
  }
  if (!NROW(ts)) {
    return(NULL)
  }
  ord <- intersect(c("M", "beta", "id", "year"), names(ts))
  ts[do.call(order, ts[ord]), , drop = FALSE]
}

#' Equilibrium locus in Kobe space (F/FMSY vs B/BMSY) from curveSS equil_yield
kobe_equilibrium_df <- function(pe, facet_cols = NULL) {
  crv <- pe$curve
  rfs <- pe$refpts
  if (is.null(crv) || !is.data.frame(crv) || !NROW(crv)) {
    return(NULL)
  }
  if (is.null(rfs) || !is.data.frame(rfs) || !NROW(rfs)) {
    return(NULL)
  }
  if (!all(c("ssb", "yield") %in% names(crv))) {
    return(NULL)
  }
  facet_cols <- facet_cols %||% production_facet_cols(pe)
  by_cols <- intersect(facet_cols, intersect(names(crv), names(rfs)))
  ref_cols <- unique(c(by_cols, "bmsy", "msy"))
  ref_cols <- intersect(ref_cols, names(rfs))

  if (length(by_cols)) {
    df <- merge(crv, unique(rfs[, ref_cols, drop = FALSE]), by = by_cols, all.x = TRUE)
  } else if ("id" %in% names(crv) && "id" %in% names(rfs)) {
    df <- merge(crv, unique(rfs[, intersect(c("id", "bmsy", "msy"), names(rfs))]), by = "id", all.x = TRUE)
  } else if (NROW(rfs) == 1L) {
    df <- cbind(crv, bmsy = rfs$bmsy[1L], msy = rfs$msy[1L])
  } else {
    return(NULL)
  }

  fill_msy <- function(x) {
    msy <- x$msy[1L]
    bmsy <- x$bmsy[1L]
    if (!is.finite(msy) || msy <= 0) {
      ok <- is.finite(x$yield) & is.finite(x$ssb)
      if (!any(ok)) {
        return(x)
      }
      ix <- which.max(x$yield[ok])
      idx <- which(ok)[ix]
      msy <- x$yield[idx]
      if (!is.finite(bmsy) || bmsy <= 0) {
        bmsy <- x$ssb[idx]
      }
    }
    if (!is.finite(bmsy) || bmsy <= 0) {
      ok <- is.finite(x$ssb) & x$ssb > 0
      if (any(ok)) {
        bmsy <- max(x$ssb[ok], na.rm = TRUE)
      }
    }
    x$msy <- msy
    x$bmsy <- bmsy
    x
  }

  split_cols <- intersect(facet_cols, names(df))
  if (length(split_cols) && requireNamespace("plyr", quietly = TRUE)) {
    df <- plyr::ddply(df, split_cols, fill_msy)
  } else {
    df <- fill_msy(df)
  }

  df$stock <- df$ssb / df$bmsy
  if ("ffmsy" %in% names(df) && any(is.finite(df$ffmsy))) {
    df$harvest <- df$ffmsy
  } else if (all(c("F", "fmsy") %in% names(df)) && any(is.finite(df$fmsy) & df$fmsy > 0)) {
    df$harvest <- df$F / df$fmsy
  } else {
    df$harvest <- df$yield / df$ssb / (df$msy / df$bmsy)
  }
  df <- df[
    is.finite(df$stock) & is.finite(df$harvest) & is.finite(df$ssb) & df$ssb > 0,
    ,
    drop = FALSE
  ]
  if (!NROW(df)) {
    return(NULL)
  }
  ord <- intersect(c(facet_cols, "stock"), names(df))
  df[do.call(order, df[ord]), , drop = FALSE]
}

#' Kobe phase plots (kobe::kobePhase) aligned with production-function facets
plot_kobe_phase_functions <- function(
  pe,
  facet_cols = NULL,
  M_filter = NULL,
  beta_filter = NULL,
  id_filter = NULL,
  xlim = NULL,
  ylim = NULL,
  title = "Kobe phase plot",
  show_trajectories = NULL,
  show_equilibrium = TRUE,
  equilibrium_colour = "steelblue",
  equilibrium_linewidth = 0.9,
  equilibrium_linetype = "solid",
  show_msy = TRUE,
  facet_scales = "fixed",
  trajectory_linewidth = 0.8,
  trajectory_colour = "grey25"
) {
  if (!require_ok(!is.null(pe), "plot_kobe_phase_functions: curveSS object missing")) {
    return(NULL)
  }
  if (!safe_require("ggplot2")) {
    return(NULL)
  }
  if (!safe_require("kobe")) {
    return(NULL)
  }

  pe <- subset_production_pe(pe, M_filter = M_filter, beta_filter = beta_filter, id_filter = id_filter)
  facet_cols <- facet_cols %||% production_facet_cols(pe)
  traj <- kobe_trajectory_df(pe)
  if (is.null(traj) || !NROW(traj)) {
    return(NULL)
  }

  if (is.null(show_trajectories)) {
    show_trajectories <- length(facet_cols) == 1L && identical(facet_cols, "id") ||
      length(facet_cols) == 2L && identical(facet_cols, c("M", "beta"))
  }
  terminal <- if (isTRUE(show_trajectories)) {
    production_terminal_rows(traj, facet_cols)
  } else {
    NULL
  }

  if (length(facet_cols) == 2L && identical(facet_cols, c("M", "beta"))) {
    facet <- ggplot2::facet_grid(
      beta ~ M,
      labeller = ggplot2::labeller(
        M = function(x) paste("M *", x),
        beta = function(x) paste("beta =", x)
      )
    )
  } else if (length(facet_cols) == 1L && facet_cols == "id") {
    facet <- ggplot2::facet_wrap(
      ~id,
      scales = facet_scales %||% "fixed",
      ncol = min(3L, length(unique(traj$id)))
    )
  } else {
    return(NULL)
  }

  phase_args <- list(object = traj[, c("stock", "harvest", "year"), drop = FALSE])
  if (!is.null(xlim) && length(xlim) == 2L) {
    phase_args$xlim <- xlim
  }
  if (!is.null(ylim) && length(ylim) == 2L) {
    phase_args$ylim <- ylim
  }
  p <- do.call(kobe::kobePhase, phase_args)
  p <- p + facet

  if (isTRUE(show_trajectories)) {
    traj$.group <- production_trajectory_group(traj, facet_cols)
    p <- p +
      ggplot2::geom_path(
        ggplot2::aes(x = stock, y = harvest, group = .group),
        data = traj,
        colour = trajectory_colour %||% "grey25",
        linewidth = trajectory_linewidth %||% 0.8,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = traj,
        ggplot2::aes(x = stock, y = harvest, group = .group),
        colour = trajectory_colour %||% "grey25",
        size = 0.6,
        alpha = 0.35,
        inherit.aes = FALSE
      )
    if (!is.null(terminal) && NROW(terminal)) {
      p <- p +
        ggplot2::geom_point(
          data = terminal,
          ggplot2::aes(x = stock, y = harvest),
          colour = "blue",
          size = 2.5,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = terminal,
          ggplot2::aes(x = stock, y = harvest, label = year),
          colour = "blue",
          size = 3,
          hjust = -0.1,
          vjust = -0.2,
          inherit.aes = FALSE
        )
    }
  }

  if (isTRUE(show_equilibrium)) {
    eq <- kobe_equilibrium_df(pe, facet_cols)
    if (!is.null(eq) && NROW(eq)) {
      eq$.eq_group <- production_trajectory_group(eq, facet_cols)
      eq$locus <- "Equilibrium"
      p <- p +
        ggplot2::geom_line(
          data = eq,
          ggplot2::aes(
            x = stock,
            y = harvest,
            group = .eq_group,
            linetype = locus
          ),
          colour = equilibrium_colour,
          linewidth = equilibrium_linewidth,
          inherit.aes = FALSE
        ) +
        ggplot2::scale_linetype_manual(
          name = NULL,
          values = c(Equilibrium = equilibrium_linetype),
          breaks = "Equilibrium"
        )
      if (isTRUE(show_msy)) {
        rfs <- pe$refpts
        msy_cols <- intersect(facet_cols, names(rfs))
        if (length(msy_cols)) {
          msy_pt <- unique(rfs[, msy_cols, drop = FALSE])
        } else {
          msy_pt <- data.frame(row = 1L)
        }
        msy_pt$stock <- 1
        msy_pt$harvest <- 1
        p <- p +
          ggplot2::geom_point(
            data = msy_pt,
            ggplot2::aes(x = stock, y = harvest),
            shape = 21,
            fill = "white",
            colour = equilibrium_colour,
            size = 3,
            inherit.aes = FALSE
          )
      }
    }
  }

  if (!identical(facet_scales, "free_y") && !identical(facet_scales, "free")) {
    p <- p + ggplot2::coord_equal()
  }

  p +
    ggplot2::labs(
      title = title,
      x = expression(B / B[MSY]),
      y = expression(F / F[MSY])
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

plot_pe_residual_fit <- function(pe, pe_mode = NULL, ylim = NULL) {
  if (!require_ok(!is.null(pe), "plot_pe_residual_fit: process-error object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  ts <- pe_residual_tseries(pe, pe_mode = pe_mode)
  if (is.null(ts)) {
    return(NULL)
  }
  resid <- ts$pe_resid %||% ts$pe2
  ylim <- .pe_plot_limits(resid, ylim = ylim)
  ylab <- if (requireNamespace("sstruth", quietly = TRUE)) {
    sstruth::processErrorLabel(pe_mode %||% Sys.getenv("DIAG_PE_MODE", unset = "relative"))
  } else {
    "Process residual"
  }
  ggplot2::ggplot(ts, ggplot2::aes(x = ssb, y = pe_resid)) +
    ggplot2::facet_grid(. ~ id, scales = "free") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.6, colour = "steelblue") +
    ggplot2::geom_smooth(se = FALSE, colour = "darkorange", linewidth = 0.8) +
    ggplot2::geom_vline(
      data = pe$refpts,
      ggplot2::aes(xintercept = bmsy),
      colour = "red",
      linetype = 2,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(x = "SSB", y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(ylim = ylim)
}

plot_pe_residual_hist <- function(pe, pe_mode = NULL, ylim = NULL) {
  if (!require_ok(!is.null(pe), "plot_pe_residual_hist: process-error object missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  ts <- pe_residual_tseries(pe, pe_mode = pe_mode)
  if (is.null(ts)) {
    return(NULL)
  }
  resid <- ts$pe_resid %||% ts$pe2
  ylim <- .pe_plot_limits(resid, ylim = ylim)
  ts_f <- ts[is.finite(resid) & abs(resid) < 1e4, , drop = FALSE]
  if (!NROW(ts_f)) {
    return(NULL)
  }
  ylab <- if (requireNamespace("sstruth", quietly = TRUE)) {
    sstruth::processErrorLabel(pe_mode %||% Sys.getenv("DIAG_PE_MODE", unset = "relative"))
  } else {
    "Process residual"
  }
  ggplot2::ggplot(ts_f, ggplot2::aes(x = pe_resid)) +
    ggplot2::facet_grid(. ~ id) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 40,
      fill = "steelblue",
      colour = "white",
      alpha = 0.7
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggplot2::stat_function(
      fun = stats::dnorm,
      args = list(
        mean = mean(ts_f$pe_resid, na.rm = TRUE),
        sd = stats::sd(ts_f$pe_resid, na.rm = TRUE)
      ),
      colour = "red",
      linewidth = 0.8
    ) +
    ggplot2::labs(x = ylab, y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = ylim)
}

plot_pe_four_panel <- function(pe, pe_mode = NULL, ylim = NULL) {
  mode <- pe_mode %||% Sys.getenv("DIAG_PE_MODE", unset = "relative")
  if (requireNamespace("sstruth", quietly = TRUE)) {
    out <- sstruth::plotSsPeFourPanel(pe, pe_mode = mode, ylim = ylim)
    return(out$panel)
  }
  p_time <- plot_pe_residual_time(pe, pe_mode = pe_mode, ylim = ylim)
  p_fit <- plot_pe_residual_fit(pe, pe_mode = pe_mode, ylim = ylim)
  p_acf <- plot_pe_residual_acf(pe, pe_mode = pe_mode)
  p_hist <- plot_pe_residual_hist(pe, pe_mode = pe_mode, ylim = ylim)
  if (is.null(p_time) || is.null(p_fit) || is.null(p_acf) || is.null(p_hist)) {
    return(NULL)
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(
      p_time / p_fit / p_acf / p_hist +
        patchwork::plot_layout(heights = c(1.2, 1.2, 0.9, 0.9))
    )
  }
  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(p_time, p_fit, p_acf, p_hist, ncol = 1, align = "v"))
  }
  p_time
}

curve_ss_paper_figure_specs <- function() {
  list(
    list(
      file = "sma-fig-ts-1.png",
      plot = plot_pe_ts,
      width = 8,
      height = 8
    ),
    list(
      file = "sma-fig-pe-phase-1.png",
      plot = plot_pe_phase,
      width = 9,
      height = 4
    ),
    list(
      file = "sma-fig-res-time-1.png",
      plot = plot_pe_residual_time,
      width = 9,
      height = 4
    ),
    list(
      file = "sma-fig-res-acf-1.png",
      plot = plot_pe_residual_acf,
      width = 9,
      height = 3
    ),
    list(
      file = "4panel.png",
      plot = plot_pe_four_panel,
      width = 9,
      height = 12
    )
  )
}

save_curve_ss_paper_figures <- function(cfg, pe, figures = NULL) {
  if (is.null(pe)) {
    return(invisible(character(0)))
  }
  fig_dir <- file.path(cfg$project_dir, "figs")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- curve_ss_paper_figure_specs()
  if (!is.null(figures) && length(figures)) {
    specs <- specs[vapply(specs, function(x) x$file %in% figures, logical(1))]
  } else if (!is.null(cfg$paper_figures) && length(cfg$paper_figures)) {
    specs <- specs[vapply(specs, function(x) x$file %in% cfg$paper_figures, logical(1))]
  }
  saved <- character(0)
  for (spec in specs) {
    plt <- spec$plot(pe)
    if (is.null(plt)) {
      next
    }
    path <- file.path(fig_dir, spec$file)
    ggplot2::ggsave(
      path,
      plot = plt,
      width = spec$width,
      height = spec$height,
      dpi = 200
    )
    saved <- c(saved, path)
  }
  invisible(saved)
}

write_process_error_figures <- function(cfg, pe, output_lines = character(0)) {
  if (is.null(pe)) {
    return(output_lines)
  }
  sc_prod <- production_plot_opts(cfg, "scenarios")
  sc_kobe <- kobe_plot_opts(cfg, "scenarios")
  pe_figs <- list(
    list(
      plot = plot_production_functions(
        pe,
        title = "Production function by scenario",
        show_trajectories = sc_prod$show_trajectories %||% TRUE,
        xlim = sc_prod$xlim,
        ylim = sc_prod$ylim,
        trajectory_yield_scale = sc_prod$trajectory_yield_scale
      ),
      file = "production-functions-scenarios.png",
      caption = paste0(
        "Production functions by SS3 scenario (curveSS; historical SSB--yield trajectory; ",
        "Kobe quadrants, MSY, 40\\% $B_\\mathrm{max}$ path)."
      ),
      label = "fig:internal-production-functions-scenarios",
      width = sc_prod$width %||% 9,
      height = sc_prod$height %||% 6,
      fig_width = 0.95
    ),
    list(
      plot = plot_kobe_phase_functions(
        pe,
        title = "Kobe phase plot by scenario",
        show_trajectories = sc_kobe$show_trajectories,
        xlim = sc_kobe$xlim,
        ylim = sc_kobe$ylim,
        show_equilibrium = sc_kobe$show_equilibrium %||% TRUE,
        facet_scales = sc_kobe$facet_scales %||% "fixed",
        trajectory_linewidth = sc_kobe$trajectory_linewidth,
        trajectory_colour = sc_kobe$trajectory_colour
      ),
      file = "kobe-phase-scenarios.png",
      caption = "Kobe phase plot by SS3 scenario (curveSS; equilibrium F/FMSY vs B/BMSY locus).",
      label = "fig:internal-kobe-phase-scenarios",
      width = sc_kobe$width %||% 9,
      height = sc_kobe$height %||% 6,
      fig_width = 0.95
    ),
    list(
      plot = plot_pe_four_panel(pe),
      file = "pe-four-panel.png",
      caption = "Four-panel process-error diagnostics (curveSS).",
      label = "fig:internal-pe-four-panel",
      width = 9,
      height = 12,
      fig_width = 0.92
    )
  )
  for (item in pe_figs) {
    fig_path <- save_appendix_figure(
      cfg,
      item$plot,
      strand = "internal",
      filename = item$file,
      width = item$width,
      height = item$height
    )
    if (!is.null(fig_path) && file.exists(fig_path)) {
      output_lines <- append_figure_tex(
        cfg,
        output_lines,
        fig_path,
        item$caption,
        item$label,
        width = item$fig_width
      )
    }
  }
  output_lines
}

backtest_params_from_cfg <- function(params) {
  longrun_params_from_cfg(params)
}

retro_dirs_available <- function(params) {
  backend_retro_dirs_available(params)
}

sp_ssb_from_timeseries <- function(rep, cfg = NULL) {
  if (!is.null(cfg)) {
    ext <- backend_load_external(cfg, backend_reference_path(cfg))
    return(ext$sp)
  }
  ts <- rep[["timeseries"]]
  if (is.null(ts) || !NROW(ts)) {
    stop("SS3 timeseries required for SP/SSB plot.", call. = FALSE)
  }
  ts <- ts[ts$Era %in% c("TIME", "INIT"), , drop = FALSE]
  if ("Seas" %in% names(ts)) {
    ts <- ts[ts$Seas == min(ts$Seas, na.rm = TRUE), , drop = FALSE]
  }
  ts <- ts[order(ts$Yr), , drop = FALSE]
  dead_cols <- grep("^dead\\(B\\):", names(ts), value = TRUE)
  catch <- if (length(dead_cols)) {
    rowSums(ts[, dead_cols, drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, NROW(ts))
  }
  yrs <- ts$Yr
  ssb <- ts$SpawnBio
  ny <- length(yrs)
  sp_ssb <- rep(NA_real_, ny)
  if (ny > 1L) {
    sp_ssb[seq_len(ny - 1L)] <- ssb[2:ny] - ssb[seq_len(ny - 1L)] + catch[seq_len(ny - 1L)]
  }
  data.frame(
    year = yrs,
    ssb = ssb,
    catchSsb = catch,
    spSsb = sp_ssb,
    stringsAsFactors = FALSE
  )
}

ss_pe_diagnose_compat <- function(path, shapePt = NULL, cfg = NULL) {
  if (!is.null(cfg)) {
    ext <- backend_load_external(cfg, path)
    return(ext$ss_pe)
  }
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required.", call. = FALSE)
  }
  cs <- sstruth::curveSS(path)
  if (is.null(shapePt)) {
    if (!is.null(cs$refpts) && "bmsy" %in% names(cs$refpts)) {
      k_ssb <- max(cs$tseries$ssb, na.rm = TRUE)
      shapePt <- as.numeric(cs$refpts$bmsy[[1L]]) / k_ssb
    }
    if (!is.numeric(shapePt) || length(shapePt) != 1L || !is.finite(shapePt) || shapePt <= 0 || shapePt >= 1) {
      shapePt <- 0.4
    }
  }
  pe <- sstruth::ssPeCompare(
    tseries = cs$tseries,
    eqlYield = cs$curve,
    shapePt = shapePt,
    makePlots = FALSE
  )
  pe$diagnostics
}

load_sstruth_external <- function(path, cfg = NULL) {
  if (is.null(cfg)) {
    stop("load_sstruth_external() requires cfg for backend dispatch.", call. = FALSE)
  }
  if (!require_ok(requireNamespace("sstruth", quietly = TRUE), "Package 'sstruth' is required.", cfg)) {
    return(NULL)
  }
  if (!require_ok(dir.exists(path), paste0("SS3 scenario path not found: ", path), cfg)) {
    return(NULL)
  }
  safe_require("ggplot2")
  backend_load_external(cfg, path)
}

load_sstruth_robustness <- function(params) {
  load_assessment_robustness(params)
}

ss_pe_diagnostics_table <- function(ss_pe) {
  if (is.null(ss_pe) || !NROW(ss_pe)) {
    return(NULL)
  }
  num_cols <- vapply(ss_pe, is.numeric, logical(1))
  out <- ss_pe
  out[num_cols] <- lapply(out[num_cols], function(x) signif(x, 4))
  out
}

plot_cpue_ts_external <- function(diag_df) {
  if (!require_ok(!is.null(diag_df) && NROW(diag_df), "plot_cpue_ts_external: diagnostic data missing")) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  dat <- diag_df[!is.na(diag_df$obs), , drop = FALSE]
  if (!require_ok(NROW(dat), "plot_cpue_ts_external: no observed CPUE values")) return(NULL)
  dat$obsStd <- ave(dat$obs, dat$name, FUN = function(x) as.numeric(scale(x)))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = year, y = obsStd, group = name)) +
    ggplot2::geom_line(ggplot2::aes(y = hat), colour = "grey60", linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = obsStd), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(y = obsStd), size = 1.2) +
    ggplot2::facet_grid(name ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x = "Year",
      y = "Standardised CPUE",
      caption = "Fleet-wise z-scores; grey = model expectation (sstruth::diagSS)."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  if (requireNamespace("mpb", quietly = TRUE)) {
    scaled <- tryCatch(
      mpb:::scale(x = dat$year, y = dat$obsStd, name = dat$name),
      error = function(e) NULL
    )
    if (!is.null(scaled) && NROW(scaled)) {
      p <- p + ggplot2::geom_line(
        data = scaled,
        ggplot2::aes(x = x, y = data),
        colour = "black",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
    }
  } else {
    p <- p + ggplot2::geom_smooth(
      ggplot2::aes(y = obsStd),
      method = "loess",
      se = FALSE,
      colour = "black",
      linewidth = 0.7,
      span = 0.6
    )
  }
  p
}

plot_sp_ssb_external <- function(sp_df, cfg = NULL) {
  if (!require_ok(!is.null(sp_df) && NROW(sp_df), "plot_sp_ssb_external: SP/SSB data missing")) return(NULL)
  if (!is.null(cfg)) {
    return(backend_plot_sp_ssb(cfg, sp_df))
  }
  if (!require_ok(requireNamespace("sstruth", quietly = TRUE), "plot_sp_ssb_external: package 'sstruth' required")) return(NULL)
  sstruth::plotSpSsb(sp_df, addSmooth = TRUE, showYears = FALSE)
}

plot_retro_traj_msy <- function(retro_data) {
  if (!require_ok(
    !is.null(retro_data) && !is.null(retro_data$trajMsyLong) && NROW(retro_data$trajMsyLong),
    "plot_retro_traj_msy: retrospective trajectory data missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("dplyr")) return(NULL)
  traj <- retro_data$trajMsyLong
  traj$metricLab <- factor(
    traj$metric,
    levels = c("bratio", "ffmsy", "catchMsyRatio"),
    labels = c("SSB / BMSY", "F / FMSY", "Catch / MSY")
  )
  ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dotted",
      colour = "grey50",
      linewidth = 0.35
    ) +
    ggplot2::geom_ribbon(
      data = traj,
      ggplot2::aes(
        x = year,
        ymin = ciLo,
        ymax = ciHi,
        fill = peel,
        group = interaction(peel, phase)
      ),
      alpha = 0.22,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = traj,
      ggplot2::aes(
        year,
        value,
        colour = peel,
        linetype = phase,
        group = interaction(peel, phase)
      ),
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    ggplot2::facet_grid(metricLab ~ scenario, scales = "free_y") +
    ggplot2::scale_linetype_manual(
      name = "Period",
      values = c(historical = "solid", forecast = "dotdash"),
      labels = c("Historical (TIME)", "Forecast (FORE)")
    ) +
    ggplot2::labs(
      title = "Retrospective trajectories vs MSY benchmarks",
      x = "Year",
      y = NULL,
      colour = "Peel",
      fill = "Peel",
      caption = "From sstruth::retroBacktestData()."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

plot_retro_ratio_horizon <- function(retro_data, fore_col, ref_col, ylab, title, cfg = NULL) {
  if (!require_ok(
    !is.null(retro_data) && !is.null(retro_data$ratioHorizon) && NROW(retro_data$ratioHorizon),
    "plot_retro_ratio_horizon: ratioHorizon data missing"
  )) return(NULL)
  if (!is.null(cfg)) {
    return(backend_plot_retro_ratios(cfg, retro_data, fore_col, ref_col, ylab, title))
  }
  if (!require_ok(requireNamespace("sstruth", quietly = TRUE), "plot_retro_ratio_horizon: package 'sstruth' required")) return(NULL)
  sstruth::plotRatios(
    retro_data$ratioHorizon,
    foreCol = fore_col,
    refCol = ref_col,
    ylab = ylab,
    title = title
  )
}

profile_root <- function(params, scenario = NULL) {
  scenario <- scenario %||% params$reference_scenario
  source(file.path(params$appendix_dir, "config-profile-grid.R"), local = TRUE)
  profile_data_root(params$appendix_dir, scenario, data_dir = params$data_dir)
}

profile_status_path <- function(params, scenario = NULL) {
  file.path(profile_root(params, scenario), "status.json")
}

write_profile_status <- function(root, status, extra = list()) {
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

read_profile_status <- function(params, scenario = NULL) {
  root <- profile_root(params, scenario)
  path_json <- file.path(root, "status.json")
  path_rds <- file.path(root, "status.rds")
  if (file.exists(path_json)) {
    if (diag_strict(params)) {
      return(read_json_strict(path_json, "profile status.json"))
    }
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      return(jsonlite::fromJSON(path_json, simplifyVector = TRUE))
    }
  }
  if (file.exists(path_rds)) {
    if (diag_strict(params)) {
      return(read_rds_strict(path_rds, "profile status.rds"))
    }
    return(readRDS(path_rds))
  }
  NULL
}

profile_component_exclude <- function() {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::profileLikelihoodExclude())
  }
  c("InitEQ_Regime", "Crash_Pen", "Equil_catch", "Parm_devs", "TOTAL")
}

extract_likelihood_components <- function(ss_rep) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::extractLikelihoodComponents(ss_rep))
  }
  ll <- ss_rep$likelihoods_used
  if (is.null(ll) || !NROW(ll)) {
    return(NULL)
  }
  comps <- rownames(ll)
  vals <- as.numeric(ll[, "values", drop = TRUE])
  data.frame(
    component = comps,
    log_likelihood = vals,
    stringsAsFactors = FALSE
  )
}

read_profile_run_dir <- function(run_dir, cfg = NULL) {
  if (!is.null(cfg)) {
    return(backend_read_profile_output(run_dir))
  }
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::readSs3Run(run_dir))
  }
  rdata <- file.path(run_dir, "ss.RData")
  if (file.exists(rdata)) {
    env <- new.env(parent = emptyenv())
    load(rdata, envir = env)
    if (exists("ss", envir = env, inherits = FALSE)) {
      return(env$ss)
    }
  }
  if (file.exists(file.path(run_dir, "Report.sso"))) {
    return(read_ss3_output(run_dir))
  }
  NULL
}

copy_profile_inputs <- function(from_dir, to_dir) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::copySsInputs(from_dir, to_dir))
  }
  copy_ss_inputs(from_dir, to_dir)
}

setup_profile_scenario <- function(src_dir, dest_dir, M, beta, ctlfile) {
  if (!requireNamespace("sstruth", quietly = TRUE)) {
    stop("Package 'sstruth' is required for profile setup.", call. = FALSE)
  }
  sstruth::setupSs3ProfileScenario(src_dir, dest_dir, M, beta, ctlfile)
}

run_profile_scenario <- function(dest_dir, ss_exe, skip_finished = TRUE, show_in_console = FALSE) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    res <- sstruth::runSs3Dir(
      dest_dir,
      ss_exe = ss_exe,
      skip_finished = skip_finished,
      show_in_console = show_in_console,
      save_output = TRUE
    )
    if (isTRUE(diag_strict()) && !isTRUE(res$ok) && !isTRUE(res$skipped)) {
      stop("SS3 profile run failed in ", dest_dir, call. = FALSE)
    }
    return(res)
  }
  if (!requireNamespace("r4ss", quietly = TRUE)) {
    stop("Package 'r4ss' is required.", call. = FALSE)
  }
  if (isTRUE(skip_finished) && file.exists(file.path(dest_dir, "Report.sso"))) {
    return(list(ok = TRUE, skipped = TRUE, dir = dest_dir))
  }
  ok <- tryCatch({
    r4ss::run(
      dest_dir,
      exe = ss_exe,
      skipfinished = isTRUE(skip_finished),
      show_in_console = isTRUE(show_in_console)
    )
    TRUE
  }, error = function(e) {
    if (diag_strict()) {
      stop("SS3 run failed in ", dest_dir, ": ", conditionMessage(e), call. = FALSE)
    }
    message("SS3 run failed in ", dest_dir, ": ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) {
    ss_rep <- read_profile_run_dir(dest_dir)
    if (!is.null(ss_rep)) {
      ss <- ss_rep
      save(ss, file = file.path(dest_dir, "ss.RData"))
    } else if (diag_strict()) {
      stop("SS3 output read failed after profile run in ", dest_dir, call. = FALSE)
    }
  } else if (diag_strict()) {
    stop("SS3 profile run failed in ", dest_dir, call. = FALSE)
  }
  list(ok = isTRUE(ok), skipped = FALSE, dir = dest_dir)
}

import_profile_grid <- function(import_root, dest_root) {
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::importProfileGrid(import_root, dest_root))
  }
  stop("Package 'sstruth' is required to import profile grids.", call. = FALSE)
}

collect_likelihood_profiles <- function(params, scenario = NULL) {
  root <- profile_root(params, scenario)
  if (requireNamespace("sstruth", quietly = TRUE)) {
    return(sstruth::collectProfileLikelihoods(root))
  }
  manifest_path <- file.path(root, "manifest.csv")
  comp_path <- file.path(root, "likelihood-by-component.csv")
  if (file.exists(comp_path)) {
    comp <- utils::read.csv(comp_path, stringsAsFactors = FALSE)
    manifest <- if (file.exists(manifest_path)) {
      utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    } else {
      NULL
    }
    return(list(root = root, manifest = manifest, components = comp))
  }

  scen_path <- file.path(root, "grid.csv")
  if (!file.exists(scen_path)) {
    return(NULL)
  }
  scen <- utils::read.csv(scen_path, stringsAsFactors = FALSE)
  rows <- list()
  manifest <- list()
  for (i in seq_len(nrow(scen))) {
    run_dir <- scen$dir[i]
    if (!dir.exists(run_dir)) {
      next
    }
    ss_rep <- read_profile_run_dir(run_dir, cfg = params)
    if (is.null(ss_rep)) {
      next
    }
    ll <- extract_likelihood_components(ss_rep)
    if (is.null(ll)) {
      next
    }
    ll$M <- scen$M[i]
    ll$beta <- scen$beta[i]
    ll$run_dir <- run_dir
    rows[[length(rows) + 1L]] <- ll
    manifest[[length(manifest) + 1L]] <- data.frame(
      M = scen$M[i],
      beta = scen$beta[i],
      dir = run_dir,
      total_ll = ll$log_likelihood[ll$component == "TOTAL"],
      ok = TRUE,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(NULL)
  }
  comp <- do.call(rbind, rows)
  manifest_df <- do.call(rbind, manifest)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(manifest_df, manifest_path, row.names = FALSE)
  utils::write.csv(comp, comp_path, row.names = FALSE)
  list(root = root, manifest = manifest_df, components = comp)
}

plot_likelihood_profile_total <- function(prof) {
  if (!require_ok(
    !is.null(prof) && !is.null(prof$components) && NROW(prof$components),
    "plot_likelihood_profile_total: profile data missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  dat <- prof$components[prof$components$component == "TOTAL", , drop = FALSE]
  if (!require_ok(NROW(dat), "plot_likelihood_profile_total: TOTAL component missing")) return(NULL)
  dat$beta_lab <- paste0("beta = ", dat$beta)
  ggplot2::ggplot(dat, ggplot2::aes(x = M, y = log_likelihood, colour = beta_lab, group = beta_lab)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "M multiplier",
      y = "Total log-likelihood",
      colour = expression(beta),
      caption = "Likelihood profile over fixed M and SR_surv_Beta (SR_surv_Beta grid)."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

plot_likelihood_profile_components <- function(prof) {
  if (!require_ok(
    !is.null(prof) && !is.null(prof$components) && NROW(prof$components),
    "plot_likelihood_profile_components: profile data missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  exclude <- profile_component_exclude()
  dat <- prof$components[!prof$components$component %in% exclude, , drop = FALSE]
  if (!require_ok(NROW(dat), "plot_likelihood_profile_components: no component rows after filtering")) return(NULL)
  dat$beta_lab <- paste0("beta = ", dat$beta)
  dat$component_lab <- gsub("_", " ", dat$component)
  ggplot2::ggplot(dat, ggplot2::aes(x = M, y = log_likelihood, colour = beta_lab, group = beta_lab)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~component_lab, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      x = "M multiplier",
      y = "Log-likelihood component",
      colour = expression(beta),
      caption = "Component likelihood profiles (ss-ll/Rmd/sma_ll.Rmd LL2 pattern)."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 8)
    )
}

profile_runs_complete <- function(params, scenario = NULL) {
  prof <- collect_likelihood_profiles(params, scenario)
  !is.null(prof) && !is.null(prof$components) && NROW(prof$components) > 0L
}

launch_likelihood_profiles <- function(params, background = TRUE, scenario = NULL) {
  scenario <- scenario %||% params$reference_scenario
  script <- normalizePath(
    file.path(params$appendix_dir, "R", "run-likelihood-profiles.R"),
    winslash = "/",
    mustWork = TRUE
  )
  root <- profile_root(params, scenario)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(root, "run.log")
  args <- c(
    "--vanilla",
    script,
    "--appendix-dir", params$appendix_dir,
    "--scenario", scenario,
    "--all"
  )
  write_profile_status(root, "queued", list(scenario = scenario, log = log_path))
  if (isTRUE(background)) {
    system2(
      normalizePath(Sys.which("Rscript"), winslash = "/", mustWork = TRUE),
      args,
      stdout = log_path,
      stderr = log_path,
      wait = FALSE
    )
    message("Likelihood profile job queued in background; log: ", log_path)
  } else {
    system2(
      normalizePath(Sys.which("Rscript"), winslash = "/", mustWork = TRUE),
      args,
      stdout = log_path,
      stderr = log_path
    )
    message("Likelihood profile job finished; log: ", log_path)
  }
  invisible(log_path)
}

longrun_job_root <- function(params, job, scenario = NULL) {
  source(file.path(params$appendix_dir, "config-longruns.R"), local = TRUE)
  longrun_data_root(params$appendix_dir, job, scenario, data_dir = params$data_dir)
}

load_collected_retrospectives <- function(params) {
  root <- longrun_job_root(params, "retrospectives")
  path <- file.path(root, "collected.rds")
  if (!file.exists(path)) {
    if (diag_strict(params)) {
      stop("Retrospectives not collected: ", path, call. = FALSE)
    }
    return(NULL)
  }
  if (diag_strict(params)) {
    return(read_rds_strict(path, "retrospectives/collected.rds"))
  }
  tryCatch(readRDS(path), error = function(e) {
    stop("Failed to read ", path, ": ", conditionMessage(e), call. = FALSE)
  })
}

collect_retrospectives <- function(params, scenarios = NULL, import_from = NULL) {
  root <- longrun_job_root(params, "retrospectives")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  scenarios <- scenarios %||% params$retro_scenarios
  bt <- backend_backtest_params(params)

  retro <- backend_collect_retrospectives(
    params, scenarios = scenarios, import_from = import_from
  )

  saveRDS(retro, file.path(root, "collected.rds"))
  if (diag_strict(params)) {
    require_data(retro$oneStep, "Retrospective one-step summary", params)
  }
  if (!is.null(retro$oneStep) && NROW(retro$oneStep)) {
    utils::write.csv(retro$oneStep, file.path(root, "one-step.csv"), row.names = FALSE)
  }
  if (!is.null(retro$ratioHorizon) && NROW(retro$ratioHorizon)) {
    utils::write.csv(retro$ratioHorizon, file.path(root, "ratio-horizon.csv"), row.names = FALSE)
  }
  manifest <- data.frame(
    scenario = scenarios,
    retro_dir = file.path(params$ss_base, scenarios, bt$retroSubdir),
    ok = dir.exists(file.path(params$ss_base, scenarios, bt$retroSubdir)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(root, "manifest.csv"), row.names = FALSE)
  list(root = root, data = retro, manifest = manifest)
}

load_collected_catch_backtests <- function(params) {
  root <- longrun_job_root(params, "catch-backtests")
  path <- file.path(root, "collected.rds")
  csv_path <- file.path(root, "forecast-summary.csv")
  if (file.exists(path)) {
    if (diag_strict(params)) {
      return(read_rds_strict(path, "catch-backtests/collected.rds"))
    }
    return(tryCatch(readRDS(path), error = function(e) {
      stop("Failed to read ", path, ": ", conditionMessage(e), call. = FALSE)
    }))
  }
  if (file.exists(csv_path)) {
    return(list(root = root, summary = utils::read.csv(csv_path, stringsAsFactors = FALSE)))
  }
  if (diag_strict(params)) {
    stop("Catch back-tests not collected: ", csv_path, call. = FALSE)
  }
  NULL
}

collect_catch_backtests <- function(params, scenarios = NULL, lr = NULL) {
  root <- longrun_job_root(params, "catch-backtests")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (is.null(lr)) {
    source(file.path(params$appendix_dir, "config-longruns.R"), local = TRUE)
    source(file.path(params$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
    lr <- longrun_defaults_from_env(longrun_defaults())
  }

  collected <- backend_collect_catch_backtests(params, scenarios = scenarios, lr = lr)
  if (is.null(collected)) {
    return(NULL)
  }
  summary_df <- collected$summary
  manifest_df <- collected$manifest
  utils::write.csv(summary_df, file.path(root, "forecast-summary.csv"), row.names = FALSE)
  utils::write.csv(manifest_df, file.path(root, "manifest.csv"), row.names = FALSE)
  out <- list(root = root, summary = summary_df, manifest = manifest_df)
  saveRDS(out, file.path(root, "collected.rds"))
  out
}

load_collected_crossvalidation <- function(params, scenario = NULL) {
  scenario <- scenario %||% params$reference_scenario
  root <- longrun_job_root(params, "crossvalidation", scenario)
  path <- file.path(root, "results.csv")
  if (file.exists(path)) {
    return(list(
      root = root,
      results = utils::read.csv(path, stringsAsFactors = FALSE),
      reference = NULL
    ))
  }
  rds <- file.path(root, "collected.rds")
  if (file.exists(rds)) {
    if (diag_strict(params)) {
      return(read_rds_strict(rds, "crossvalidation/collected.rds"))
    }
    return(readRDS(rds))
  }
  NULL
}

collect_crossvalidation <- function(params, scenario = NULL) {
  scenario <- scenario %||% params$reference_scenario
  source(file.path(params$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  root <- longrun_job_root(params, "crossvalidation", scenario)
  grid_path <- file.path(root, "grid.csv")
  if (!file.exists(grid_path)) {
    return(NULL)
  }
  grid <- utils::read.csv(grid_path, stringsAsFactors = FALSE)
  ref_dir <- reference_path(params)
  ref_rep <- backend_read_output(params, ref_dir)
  ref_ll <- if (!is.null(ref_rep) && !is.null(ref_rep$likelihoods_used)) {
    as.numeric(ref_rep$likelihoods_used["TOTAL", "values"])
  } else {
    NA_real_
  }
  ref_stat <- terminal_bratio_f(ref_rep)

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    run_dir <- grid$dir[i]
    ss_rep <- read_profile_run_dir(run_dir, cfg = params)
    if (is.null(ss_rep)) {
      return(NULL)
    }
    ll <- if (!is.null(ss_rep$likelihoods_used)) {
      as.numeric(ss_rep$likelihoods_used["TOTAL", "values"])
    } else {
      NA_real_
    }
    st <- terminal_bratio_f(ss_rep)
    data.frame(
      fleet_id = grid$fleet_id[i],
      fleet_name = grid$fleet_name[i],
      run_dir = run_dir,
      total_ll = ll,
      delta_ll = ll - ref_ll,
      bratio = st$bratio,
      delta_bratio = st$bratio - ref_stat$bratio,
      ffmsy = st$ffmsy,
      delta_ffmsy = st$ffmsy - ref_stat$ffmsy,
      ok = TRUE,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(NULL)
  }
  results <- do.call(rbind, rows)
  ref_row <- data.frame(
    fleet_id = NA_integer_,
    fleet_name = "reference",
    run_dir = ref_dir,
    total_ll = ref_ll,
    delta_ll = 0,
    bratio = ref_stat$bratio,
    delta_bratio = 0,
    ffmsy = ref_stat$ffmsy,
    delta_ffmsy = 0,
    ok = TRUE,
    stringsAsFactors = FALSE
  )
  utils::write.csv(results, file.path(root, "results.csv"), row.names = FALSE)
  out <- list(root = root, results = results, reference = ref_row)
  saveRDS(out, file.path(root, "collected.rds"))
  out
}

longrun_runs_complete <- function(params, job, scenario = NULL) {
  root <- longrun_job_root(params, job, scenario)
  path_json <- file.path(root, "status.json")
  if (file.exists(path_json) && requireNamespace("jsonlite", quietly = TRUE)) {
    st <- if (diag_strict(params)) {
      read_json_strict(path_json, "longrun job status.json")
    } else {
      jsonlite::fromJSON(path_json, simplifyVector = TRUE)
    }
    if (!is.null(st) && identical(st$status, "complete")) {
      return(TRUE)
    }
  }
  if (identical(job, "crossvalidation")) {
    return(file.exists(file.path(root, "results.csv")))
  }
  file.exists(file.path(root, "collected.rds"))
}

launch_retrospectives <- function(params, background = TRUE, run = TRUE, collect = TRUE) {
  source(file.path(params$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  root <- longrun_job_root(params, "retrospectives")
  args <- c(if (isTRUE(run) && isTRUE(collect)) "--all" else {
    c(if (isTRUE(run)) "--run", if (isTRUE(collect)) "--collect")
  })
  script <- file.path(params$appendix_dir, "R", "run-retrospectives.R")
  if (isTRUE(background)) {
    launch_longrun_job(script, params$appendix_dir, args, root, "Retrospectives")
  } else {
    system2(Sys.which("Rscript"), c("--vanilla", script, "--appendix-dir", params$appendix_dir, args))
  }
}

launch_catch_backtests <- function(params, background = TRUE, run = TRUE, collect = TRUE) {
  source(file.path(params$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  root <- longrun_job_root(params, "catch-backtests")
  args <- c(if (isTRUE(run) && isTRUE(collect)) "--all" else {
    c(if (isTRUE(run)) "--run", if (isTRUE(collect)) "--collect")
  })
  script <- file.path(params$appendix_dir, "R", "run-catch-backtests.R")
  if (isTRUE(background)) {
    launch_longrun_job(script, params$appendix_dir, args, root, "Catch back-tests")
  } else {
    system2(Sys.which("Rscript"), c("--vanilla", script, "--appendix-dir", params$appendix_dir, args))
  }
}

launch_crossvalidation <- function(params, background = TRUE, scenario = NULL) {
  source(file.path(params$appendix_dir, "R", "longrun-utils.R"), local = TRUE)
  scenario <- scenario %||% params$reference_scenario
  root <- longrun_job_root(params, "crossvalidation", scenario)
  script <- file.path(params$appendix_dir, "R", "run-crossvalidation.R")
  args <- c("--scenario", scenario, "--all")
  if (isTRUE(background)) {
    launch_longrun_job(script, params$appendix_dir, args, root, "LOO cross-validation")
  } else {
    system2(Sys.which("Rscript"), c("--vanilla", script, "--appendix-dir", params$appendix_dir, args))
  }
}

mvln_cfg <- function(params) {
  defaults <- if (exists("mvln_defaults", mode = "function")) {
    mvln_defaults()
  } else {
    list(mc = 500L, seed = 1224L, terminal_years = NULL, equal_weights = TRUE)
  }
  defaults
}

mvln_scenario_draws <- function(path, mc = 500L, seed = NULL) {
  if (!require_ok(
    requireNamespace("mvtnorm", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE),
    "mvln_scenario_draws: packages 'mvtnorm' and 'reshape2' required"
  )) return(NULL)
  rep <- read_ss3_output(path, covar = TRUE)
  if (!require_ok(
    !is.null(rep[["CoVar"]]) && NROW(rep[["CoVar"]]),
    paste0("No covariance matrix for MVLN in ", path)
  )) return(NULL)
  draws <- ssmvln(rep[["CoVar"]], rep[["derived_quants"]], mc = mc, seed = seed)
  if (!require_ok(!is.null(draws), paste0("MVLN sampling failed for ", path))) return(NULL)
  kobe <- mvln_long_to_kobe(draws, fmsy_from_rep(rep))
  list(
    scenario = basename(path),
    path = path,
    rep = rep,
    draws = draws,
    kobe = kobe,
    fmsy = fmsy_from_rep(rep),
    endyr = rep$endyr %||% max(draws$year, na.rm = TRUE)
  )
}

ensemble_scenario_weights <- function(scenarios, weights = NULL) {
  scenarios <- as.character(scenarios)
  if (is.null(weights)) {
    w <- rep(1 / length(scenarios), length(scenarios))
  } else {
    w <- as.numeric(weights)
    if (length(w) != length(scenarios)) {
      stop("ensemble weights length must match scenarios.", call. = FALSE)
    }
    if (any(!is.finite(w)) || sum(w) <= 0) {
      w <- rep(1 / length(scenarios), length(scenarios))
    } else {
      w <- w / sum(w)
    }
  }
  stats::setNames(w, scenarios)
}

load_mvln_ensemble <- function(params, mc = NULL, seed = NULL) {
  if (!require_ok(requireNamespace("mvtnorm", quietly = TRUE), "Package 'mvtnorm' is required for MVLN ensemble.", params)) {
    return(NULL)
  }
  cfg <- mvln_cfg(params)
  if (is.null(mc)) {
    mc <- cfg$mc
  }
  if (is.null(seed)) {
    seed <- cfg$seed
  }
  paths <- scenario_paths(params)
  paths <- paths[dir.exists(paths)]
  if (!require_ok(length(paths) >= 2L, "MVLN ensemble needs at least two scenario directories.", params)) {
    return(NULL)
  }

  rows <- lapply(paths, function(p) {
    mvln_scenario_draws(p, mc = mc, seed = seed)
  })
  if (diag_strict(params)) {
    failed <- vapply(rows, is.null, logical(1))
    if (any(failed)) {
      stop(
        "MVLN draws failed for: ",
        paste(basename(paths)[failed], collapse = ", "),
        call. = FALSE
      )
    }
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!require_ok(length(rows) >= 2L, "MVLN draws available for fewer than two scenarios.", params)) {
    return(NULL)
  }

  kobe <- lapply(rows, function(x) {
    if (is.null(x$kobe) || !NROW(x$kobe)) {
      return(NULL)
    }
    cbind(scenario = x$scenario, x$kobe, stringsAsFactors = FALSE)
  })
  kobe <- kobe[!vapply(kobe, is.null, logical(1))]
  if (!length(kobe)) {
    return(NULL)
  }
  kobe_df <- do.call(rbind, kobe)

  ref <- params$reference_scenario
  ref_row <- rows[vapply(rows, function(x) identical(x$scenario, ref), logical(1))]
  terminal_year <- if (!is.null(cfg$terminal_years) && length(cfg$terminal_years)) {
    max(cfg$terminal_years)
  } else if (length(ref_row)) {
    ref_row[[1]]$endyr
  } else {
    terminal_yrs <- stats::aggregate(year ~ scenario, data = kobe_df, FUN = max)
    min(terminal_yrs$year, na.rm = TRUE)
  }

  terminal <- kobe_df[kobe_df$year == terminal_year, , drop = FALSE]
  if (!NROW(terminal)) {
    terminal_year <- max(kobe_df$year, na.rm = TRUE)
    terminal <- kobe_df[kobe_df$year == terminal_year, , drop = FALSE]
  }

  summary_terminal <- do.call(rbind, lapply(split(terminal, terminal$scenario), function(dat) {
    data.frame(
      scenario = dat$scenario[[1]],
      terminal_year = terminal_year,
      stock_median = stats::median(dat$stock, na.rm = TRUE),
      stock_p05 = stats::quantile(dat$stock, 0.05, na.rm = TRUE),
      stock_p95 = stats::quantile(dat$stock, 0.95, na.rm = TRUE),
      harvest_median = stats::median(dat$harvest, na.rm = TRUE),
      harvest_p05 = stats::quantile(dat$harvest, 0.05, na.rm = TRUE),
      harvest_p95 = stats::quantile(dat$harvest, 0.95, na.rm = TRUE),
      p_overfished = mean(dat$stock < 1, na.rm = TRUE),
      p_overfishing = mean(dat$harvest > 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  weights <- resolve_ensemble_weights(params, summary_terminal$scenario)
  summary_terminal$weight <- weights[summary_terminal$scenario]

  list(
    kobe = kobe_df,
    terminal = terminal,
    terminal_year = terminal_year,
    summary = summary_terminal,
    weights = weights,
    reference = ref,
    mc = mc,
    scenarios = rows
  )
}

mvln_summary_table <- function(mvln) {
  if (is.null(mvln) || is.null(mvln$summary) || !NROW(mvln$summary)) {
    return(NULL)
  }
  out <- mvln$summary
  num_cols <- c(
    "stock_median", "stock_p05", "stock_p95",
    "harvest_median", "harvest_p05", "harvest_p95",
    "p_overfished", "p_overfishing", "weight"
  )
  for (col in num_cols) {
    if (col %in% names(out)) {
      out[[col]] <- signif(out[[col]], 3)
    }
  }
  out
}

plot_mvln_kobe_ensemble <- function(mvln, max_points = 400L) {
  if (!require_ok(
    !is.null(mvln) && !is.null(mvln$terminal) && NROW(mvln$terminal),
    "plot_mvln_kobe_ensemble: MVLN terminal draws missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  dat <- mvln$terminal
  if (NROW(dat) > max_points) {
    set.seed(mvln$mc %||% 1L)
    dat <- dat[sample.int(NROW(dat), max_points), , drop = FALSE]
  }
  dat$scenario <- factor(
    dat$scenario,
    levels = c(mvln$reference, setdiff(unique(dat$scenario), mvln$reference))
  )
  ref_dat <- dat[dat$scenario == mvln$reference, , drop = FALSE]
  ggplot2::ggplot(dat, ggplot2::aes(x = stock, y = harvest, colour = scenario)) +
    ggplot2::geom_point(alpha = 0.25, size = 1.2) +
    ggplot2::geom_point(
      data = ref_dat,
      ggplot2::aes(x = stock, y = harvest),
      colour = "black",
      alpha = 0.35,
      size = 1.4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    ggplot2::labs(
      title = "MVLN terminal status: best assessment vs scenarios",
      subtitle = paste0("Terminal year ", mvln$terminal_year, "; harvest = F / FMSY"),
      x = expression(B / B[MSY]),
      y = expression(F / F[MSY]),
      colour = "Scenario",
      caption = "Hessian-based MVLN draws from SS3 covariance (ssmvln)."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

plot_mvln_terminal_compare <- function(mvln) {
  if (!require_ok(
    !is.null(mvln) && !is.null(mvln$summary) && NROW(mvln$summary),
    "plot_mvln_terminal_compare: MVLN summary missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("reshape2")) return(NULL)
  dat <- mvln$summary
  dat$scenario <- factor(
    dat$scenario,
    levels = c(mvln$reference, setdiff(dat$scenario, mvln$reference))
  )
  long <- reshape2::melt(
    dat[, c("scenario", "stock_median", "stock_p05", "stock_p95",
            "harvest_median", "harvest_p05", "harvest_p95")],
    id.vars = "scenario"
  )
  long$metric <- ifelse(grepl("^stock", long$variable), "SSB / BMSY", "F / FMSY")
  long$stat <- ifelse(grepl("_median$", long$variable), "median",
    ifelse(grepl("_p05$", long$variable), "p05", "p95"))
  med <- long[long$stat == "median", , drop = FALSE]
  rib <- reshape2::dcast(long, scenario + metric ~ stat, value.var = "value")
  ggplot2::ggplot() +
    ggplot2::geom_errorbar(
      data = rib,
      ggplot2::aes(
        x = scenario,
        ymin = p05,
        ymax = p95,
        colour = metric
      ),
      width = 0.15,
      linewidth = 0.8,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = scenario, y = value, colour = metric),
      size = 2.5,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::labs(
      title = "Terminal MVLN status by scenario",
      subtitle = paste0("Year ", mvln$terminal_year, "; ribbons = 5--95% MVLN interval"),
      x = NULL,
      y = NULL,
      colour = NULL,
      caption = "Reference scenario compared with structural sensitivities."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

plot_mvln_status_ts <- function(mvln) {
  if (!require_ok(
    !is.null(mvln) && !is.null(mvln$kobe) && NROW(mvln$kobe),
    "plot_mvln_status_ts: MVLN kobe time series missing"
  )) return(NULL)
  if (!safe_require("ggplot2")) return(NULL)
  if (!safe_require("dplyr")) return(NULL)
  yrs <- sort(unique(mvln$kobe$year))
  yrs <- tail(yrs, min(15L, length(yrs)))
  dat <- mvln$kobe[mvln$kobe$year %in% yrs, , drop = FALSE]
  q <- dplyr::group_by(dat, scenario, year)
  q <- dplyr::summarise(
    q,
    stock_med = stats::median(stock, na.rm = TRUE),
    stock_lo = stats::quantile(stock, 0.05, na.rm = TRUE),
    stock_hi = stats::quantile(stock, 0.95, na.rm = TRUE),
    harvest_med = stats::median(harvest, na.rm = TRUE),
    harvest_lo = stats::quantile(harvest, 0.05, na.rm = TRUE),
    harvest_hi = stats::quantile(harvest, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
  q$scenario <- factor(
    q$scenario,
    levels = c(mvln$reference, setdiff(unique(q$scenario), mvln$reference))
  )
  stock_plot <- ggplot2::ggplot(q, ggplot2::aes(x = year, colour = scenario, fill = scenario)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = stock_lo, ymax = stock_hi),
      alpha = 0.15,
      colour = NA
    ) +
    ggplot2::geom_line(ggplot2::aes(y = stock_med), linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(
      title = "SSB / BMSY (MVLN)",
      x = "Year",
      y = expression(B / B[MSY]),
      caption = "5--95% MVLN intervals by scenario."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  harvest_plot <- ggplot2::ggplot(q, ggplot2::aes(x = year, colour = scenario, fill = scenario)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = harvest_lo, ymax = harvest_hi),
      alpha = 0.15,
      colour = NA
    ) +
    ggplot2::geom_line(ggplot2::aes(y = harvest_med), linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(
      title = "F / FMSY (MVLN)",
      x = "Year",
      y = expression(F / F[MSY])
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(stock_plot, harvest_plot, ncol = 1))
  }
  stock_plot
}
