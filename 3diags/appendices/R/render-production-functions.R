#!/usr/bin/env Rscript
# Generate curveSS production-function figures for a study.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  appendix_dir <- "appendices"
  study_id <- Sys.getenv("DIAG_STUDY", unset = "sma-2026")
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--appendix-dir" && i < length(args)) {
      appendix_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] == "--study" && i < length(args)) {
      study_id <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }

  appendix_dir <- normalizePath(appendix_dir, winslash = "/", mustWork = TRUE)
  source(file.path(appendix_dir, "R", "appendix-utils.R"))
  source(file.path(appendix_dir, "config-case-study.R"))

  cfg <- case_study_defaults(appendix_dir, study_id)
  load_appendix_helpers(cfg)

  fig_dir <- file.path(appendix_dir, "figs", "internal")
  paper_dir <- file.path(cfg$project_dir, "figs")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

  saved <- character(0)

  ss3 <- load_assessment_case_study(cfg)
  if (!is.null(ss3$pe)) {
    sc_opts <- production_plot_opts(cfg, "scenarios")
    p_sc <- plot_production_functions(
      ss3$pe,
      title = "Production function by scenario",
      show_trajectories = sc_opts$show_trajectories %||% TRUE,
      xlim = sc_opts$xlim,
      ylim = sc_opts$ylim,
      trajectory_yield_scale = sc_opts$trajectory_yield_scale
    )
    if (!is.null(p_sc)) {
      path <- file.path(fig_dir, "production-functions-scenarios.png")
      ggplot2::ggsave(
        path,
        p_sc,
        width = sc_opts$width %||% 9,
        height = sc_opts$height %||% 6,
        dpi = 200
      )
      saved <- c(saved, path)
      message("Wrote ", path)
    }
    kb_sc <- kobe_plot_opts(cfg, "scenarios")
    p_kobe_sc <- plot_kobe_phase_functions(
      ss3$pe,
      title = "Kobe phase plot by scenario",
      show_trajectories = kb_sc$show_trajectories,
      xlim = kb_sc$xlim,
      ylim = kb_sc$ylim,
      show_equilibrium = kb_sc$show_equilibrium %||% TRUE,
      facet_scales = kb_sc$facet_scales %||% "fixed",
      trajectory_linewidth = kb_sc$trajectory_linewidth,
      trajectory_colour = kb_sc$trajectory_colour
    )
    if (!is.null(p_kobe_sc)) {
      path <- file.path(fig_dir, "kobe-phase-scenarios.png")
      ggplot2::ggsave(
        path,
        p_kobe_sc,
        width = kb_sc$width %||% 9,
        height = kb_sc$height %||% 6,
        dpi = 200
      )
      saved <- c(saved, path)
      message("Wrote ", path)
    }
  }

  prof <- collect_likelihood_profiles(cfg)
  if (!is.null(prof)) {
    prof_pe <- load_profile_pe(cfg, prof, maxY = 2.5)
    if (!is.null(prof_pe)) {
      grid <- profile_grid_from_cfg(cfg)
      m_vals <- grid$M %||% seq(0.5, 1.5, 0.25)
      pr_opts <- production_plot_opts(cfg, "profile")
      p_pr <- plot_production_functions(
        prof_pe,
        M_filter = m_vals,
        title = "Production function (M x SR_surv_Beta profile grid)",
        show_trajectories = pr_opts$show_trajectories,
        xlim = pr_opts$xlim,
        ylim = pr_opts$ylim,
        trajectory_yield_scale = pr_opts$trajectory_yield_scale
      )
      if (!is.null(p_pr)) {
        path <- file.path(fig_dir, "production-functions-profile.png")
        ggplot2::ggsave(path, p_pr, width = 9, height = 6, dpi = 200)
        saved <- c(saved, path)
        message("Wrote ", path)
      }
      kb_pr <- kobe_plot_opts(cfg, "profile")
      p_kobe_pr <- plot_kobe_phase_functions(
        prof_pe,
        M_filter = m_vals,
        title = "Kobe phase plot (M x SR_surv_Beta profile grid)",
        show_trajectories = kb_pr$show_trajectories,
        xlim = kb_pr$xlim,
        ylim = kb_pr$ylim,
        show_equilibrium = kb_pr$show_equilibrium %||% TRUE
      )
      if (!is.null(p_kobe_pr)) {
        path <- file.path(fig_dir, "kobe-phase-profile.png")
        ggplot2::ggsave(path, p_kobe_pr, width = 9, height = 6, dpi = 200)
        saved <- c(saved, path)
        message("Wrote ", path)
      }
    }
  }

  if (!length(saved)) {
    stop("No production-function figures were produced.", call. = FALSE)
  }
  invisible(saved)
}

if (sys.nframe() == 0L) {
  main()
}
