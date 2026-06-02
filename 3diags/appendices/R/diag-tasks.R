# Modular diagnostic task registry: SS3 setup, parallel runs, and figure exports.
#
# Used by appendices/R/run-diagnostics.R. Tasks are composable; presets bundle
# common first-pass workflows (profiles + PE figures).

diag_task_presets <- function() {
  list(
    first_pass = c(
      "profiles.setup",
      "profiles.run",
      "profiles.collect",
      "pe.load",
      "plot.production_scenarios",
      "plot.production_profile",
      "plot.kobe_scenarios",
      "plot.kobe_profile",
      "plot.ll_components",
      "plot.pe_four_panel"
    ),
    first_pass_import = c(
      "profiles.import",
      "pe.load",
      "plot.production_scenarios",
      "plot.production_profile",
      "plot.kobe_scenarios",
      "plot.kobe_profile",
      "plot.ll_components",
      "plot.pe_four_panel"
    ),
    figures = c(
      "pe.load",
      "plot.production_scenarios",
      "plot.production_profile",
      "plot.kobe_scenarios",
      "plot.kobe_profile",
      "plot.ll_components",
      "plot.pe_four_panel"
    ),
    figures_scenarios = c(
      "pe.load",
      "plot.production_scenarios",
      "plot.kobe_scenarios",
      "plot.pe_four_panel"
    ),
    process_error = c(
      "pe.load",
      "plot.production_scenarios",
      "plot.kobe_scenarios",
      "plot.pe_four_panel"
    ),
    profiles = c("profiles.setup", "profiles.run", "profiles.collect")
  )
}

parse_diag_args <- function(args) {
  out <- list(
    appendix_dir = NULL,
    study_id = Sys.getenv("DIAG_STUDY", unset = "sma-2026"),
    ss_base = NULL,
    scenarios = NULL,
    profile_scenario = NULL,
    import_profiles = NULL,
    tasks = "first_pass",
    plot_only = FALSE,
    parallel = TRUE,
    parallel_curvess = NULL,
    workers = NULL,
    maxY = 2.5
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
    } else if (key == "--scenarios" && i < length(args)) {
      out$scenarios <- strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
      out$scenarios <- trimws(out$scenarios)
      i <- i + 2L
    } else if (key == "--profile-scenario" && i < length(args)) {
      out$profile_scenario <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--import-profiles" && i < length(args)) {
      out$import_profiles <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--tasks" && i < length(args)) {
      out$tasks <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--workers" && i < length(args)) {
      out$workers <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--maxY" && i < length(args)) {
      out$maxY <- as.numeric(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--plot-only") {
      out$plot_only <- TRUE
      i <- i + 1L
    } else if (key == "--no-parallel") {
      out$parallel <- FALSE
      i <- i + 1L
    } else if (key == "--parallel-curvess") {
      out$parallel_curvess <- TRUE
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}

resolve_diag_task_list <- function(args) {
  if (isTRUE(args$plot_only)) {
    return(diag_task_presets()[["figures"]])
  }
  if (!is.null(args$import_profiles)) {
    return(diag_task_presets()[["first_pass_import"]])
  }
  raw <- args$tasks
  if (is.null(raw) || !nzchar(raw)) {
    return(diag_task_presets()[["first_pass"]])
  }
  key <- gsub("-", "_", tolower(raw))
  presets <- diag_task_presets()
  if (key %in% names(presets)) {
    return(presets[[key]])
  }
  strsplit(raw, ",", fixed = TRUE)[[1]]
}

init_diag_context <- function(args) {
  if (is.null(args$appendix_dir)) {
    stop("--appendix-dir is required.", call. = FALSE)
  }
  appendix_dir <- normalizePath(args$appendix_dir, winslash = "/", mustWork = TRUE)
  core <- diags::appendixCoreDir(appendix_dir)
  source(file.path(core, "R", "longrun-utils.R"), local = FALSE)
  source(file.path(core, "R", "appendix-utils.R"), local = FALSE)
  source(file.path(appendix_dir, "config-case-study.R"), local = FALSE)
  source(file.path(appendix_dir, "config-profile-grid.R"), local = FALSE)
  source(file.path(core, "R", "run-likelihood-profiles.R"), local = FALSE)

  cfg <- longrun_cfg(appendix_dir, args$study_id)
  load_appendix_helpers(cfg)
  source(file.path(core, "R", "backend.R"), local = FALSE)

  if (!is.null(args$ss_base) && nzchar(args$ss_base)) {
    base <- normalizePath(args$ss_base, winslash = "/", mustWork = FALSE)
    cfg$ss_base <- base
    cfg$assessment_base <- base
  }
  if (!is.null(args$scenarios) && length(args$scenarios)) {
    cfg$scenarios <- as.character(args$scenarios)
  }
  profile_scenario <- args$profile_scenario %||% cfg$profile_scenario %||% cfg$reference_scenario
  cfg$profile_scenario <- profile_scenario
  cfg$reference_scenario <- profile_scenario

  scenario_paths(cfg)
  if (!is.null(args$import_profiles) && nzchar(args$import_profiles)) {
    cfg$profile_import_from <- normalizePath(args$import_profiles, winslash = "/", mustWork = FALSE)
  }

  list(
    args = args,
    cfg = cfg,
    profile_scenario = profile_scenario,
    cache = new.env(parent = emptyenv()),
    results = list()
  )
}

ctx_cache <- function(ctx, key, value = NULL) {
  if (!missing(value)) {
    ctx$cache[[key]] <- value
    return(value)
  }
  if (exists(key, envir = ctx$cache, inherits = FALSE)) {
    return(get(key, envir = ctx$cache))
  }
  NULL
}

diag_parallel_curvess <- function(ctx) {
  if (!is.null(ctx$args$parallel_curvess)) {
    return(isTRUE(ctx$args$parallel_curvess))
  }
  isTRUE(ctx$args$parallel) && identical(Sys.info()[["sysname"]], "Linux")
}

profile_grid_scenario_df <- function(ctx) {
  cached <- ctx_cache(ctx, "profile_grid")
  if (!is.null(cached)) {
    return(cached)
  }
  grid_cfg <- utils::modifyList(profile_grid_defaults(), profile_grid_from_cfg(ctx$cfg))
  if (nzchar(Sys.getenv("DIAG_PROFILE_M", unset = ""))) {
    grid_cfg$M <- as.numeric(strsplit(Sys.getenv("DIAG_PROFILE_M"), ",", fixed = TRUE)[[1]])
  }
  if (nzchar(Sys.getenv("DIAG_PROFILE_BETA", unset = ""))) {
    grid_cfg$beta <- as.numeric(strsplit(Sys.getenv("DIAG_PROFILE_BETA"), ",", fixed = TRUE)[[1]])
  }
  if (nzchar(Sys.getenv("SS3_EXE", unset = ""))) {
    grid_cfg$ss_exe <- Sys.getenv("SS3_EXE")
  }
  scen <- profile_grid_scenarios(M = grid_cfg$M, beta = grid_cfg$beta)
  root <- profile_root(ctx$cfg, ctx$profile_scenario)
  scen$dir <- file.path(root, as.character(scen$M), as.character(scen$beta))
  ctx_cache(ctx, "profile_grid", list(cfg = grid_cfg, scen = scen, root = root))
}

task_profiles_setup <- function(ctx) {
  pg <- profile_grid_scenario_df(ctx)
  src_dir <- reference_path(ctx$cfg)
  if (!dir.exists(src_dir)) {
    stop("Reference scenario not found: ", src_dir, call. = FALSE)
  }
  root <- pg$root
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  files_dir <- file.path(root, "files")
  write_profile_status(root, "setup", list(n_scenarios = nrow(pg$scen)))
  copy_profile_inputs(src_dir, files_dir)
  for (i in seq_len(nrow(pg$scen))) {
    setup_profile_scenario(
      src_dir = files_dir,
      dest_dir = pg$scen$dir[i],
      M = pg$scen$M[i],
      beta = pg$scen$beta[i],
      ctlfile = pg$cfg$ctlfile
    )
  }
  utils::write.csv(
    pg$scen[, c("i", "M", "beta", "dir")],
    file.path(root, "grid.csv"),
    row.names = FALSE
  )
  message("[profiles.setup] ", nrow(pg$scen), " cells under ", root)
  invisible(pg)
}

task_profiles_run <- function(ctx) {
  pg <- profile_grid_scenario_df(ctx)
  scen <- pg$scen
  write_profile_status(pg$root, "running", list(n_scenarios = nrow(scen)))
  run_one <- function(i) {
    run_profile_scenario(
      scen$dir[i],
      ss_exe = pg$cfg$ss_exe,
      skip_finished = pg$cfg$skip_finished,
      show_in_console = pg$cfg$show_in_console
    )
  }
  results <- parallel_lapply(
    seq_len(nrow(scen)),
    run_one,
    parallel = ctx$args$parallel,
    workers = ctx$args$workers %||% pg$cfg$workers
  )
  n_ok <- sum(vapply(results, function(x) isTRUE(x$ok), logical(1)))
  write_profile_status(pg$root, "ran", list(n_ok = n_ok, n_total = nrow(scen)))
  message("[profiles.run] ", n_ok, "/", nrow(scen), " succeeded")
  invisible(results)
}

task_profiles_collect <- function(ctx) {
  prof <- collect_likelihood_profiles(ctx$cfg, ctx$profile_scenario)
  if (is.null(prof)) {
    stop(
      "No profile results under ",
      profile_root(ctx$cfg, ctx$profile_scenario),
      call. = FALSE
    )
  }
  ctx_cache(ctx, "prof", prof)
  write_profile_status(
    profile_root(ctx$cfg, ctx$profile_scenario),
    "complete",
    list(
      n_manifest = nrow(prof$manifest),
      n_component_rows = nrow(prof$components)
    )
  )
  message(
    "[profiles.collect] ",
    nrow(prof$manifest),
    " runs; ",
    nrow(prof$components),
    " component rows"
  )
  invisible(prof)
}

task_profiles_import <- function(ctx) {
  import_from <- ctx$args$import_profiles %||% ctx$cfg$profile_import_from
  if (is.null(import_from) || !nzchar(import_from)) {
    stop("Specify --import-profiles or study jobs$likelihood_profiles$import_from.", call. = FALSE)
  }
  root <- profile_root(ctx$cfg, ctx$profile_scenario)
  write_profile_status(root, "importing", list(from = import_from))
  prof <- import_profile_grid(import_from, root)
  ctx_cache(ctx, "prof", prof)
  write_profile_status(
    root,
    "complete",
    list(
      imported_from = import_from,
      n_manifest = nrow(prof$manifest),
      n_component_rows = nrow(prof$components)
    )
  )
  message("[profiles.import] from ", import_from)
  invisible(prof)
}

task_pe_load <- function(ctx) {
  paths <- backend_scenario_paths(ctx$cfg)
  paths <- paths[dir.exists(paths)]
  has_report <- vapply(paths, function(p) file.exists(file.path(p, "Report.sso")), logical(1))
  if (any(!has_report)) {
    skipped <- basename(paths[!has_report])
    message("[pe.load] skipping (no Report.sso): ", paste(skipped, collapse = ", "))
    paths <- paths[has_report]
  }
  if (!length(paths)) {
    stop("No converged SS3 scenarios found under ", ctx$cfg$ss_base, call. = FALSE)
  }
  if (!length(paths)) {
    stop("No scenario directories found under ", ctx$cfg$ss_base, call. = FALSE)
  }
  pe <- backend_load_pe(
    ctx$cfg,
    paths,
    for_plots = TRUE,
    parallel = diag_parallel_curvess(ctx),
    workers = ctx$args$workers
  )
  if (!is.null(pe)) {
    pe <- relabel_pe_ids(pe, ctx$cfg)
  }
  ss3 <- list(paths = paths, pe = pe, fls = NULL)
  ctx_cache(ctx, "ss3", ss3)
  ctx_cache(ctx, "pe", pe)
  message("[pe.load] curveSS for ", length(paths), " scenario(s)")
  invisible(ss3)
}

diag_fig_dir <- function(ctx) {
  file.path(ctx$cfg$appendix_dir, "figs", "internal")
}

save_diag_figure <- function(ctx, plot, filename, width = 9, height = 6) {
  if (is.null(plot)) {
    message("[plot] skipped (NULL): ", filename)
    return(invisible(NULL))
  }
  path <- save_appendix_figure(
    ctx$cfg,
    plot,
    strand = "internal",
    filename = filename,
    width = width,
    height = height
  )
  if (!is.null(path) && file.exists(path)) {
    message("[plot] wrote ", path)
    ctx$results[[filename]] <- path
  }
  invisible(path)
}

task_plot_production_scenarios <- function(ctx) {
  pe <- ctx_cache(ctx, "pe")
  if (is.null(pe)) {
    ss3 <- ctx_cache(ctx, "ss3")
    pe <- if (!is.null(ss3)) ss3$pe else NULL
  }
  if (is.null(pe)) {
    stop("PE data missing; run pe.load first.", call. = FALSE)
  }
  sc_opts <- production_plot_opts(ctx$cfg, "scenarios")
  p <- plot_production_functions(
    pe,
    title = "Production function by scenario",
    show_trajectories = sc_opts$show_trajectories %||% TRUE,
    xlim = sc_opts$xlim,
    ylim = sc_opts$ylim,
    trajectory_yield_scale = sc_opts$trajectory_yield_scale
  )
  save_diag_figure(
    ctx,
    p,
    "production-functions-scenarios.png",
    width = sc_opts$width %||% 9,
    height = sc_opts$height %||% 6
  )
}

task_plot_production_profile <- function(ctx) {
  prof <- ctx_cache(ctx, "prof") %||% collect_likelihood_profiles(ctx$cfg, ctx$profile_scenario)
  if (is.null(prof)) {
    stop("Profile data missing; run profiles.collect or profiles.import first.", call. = FALSE)
  }
  ctx_cache(ctx, "prof", prof)
  prof_pe <- load_profile_pe(
    ctx$cfg,
    prof,
    maxY = ctx$args$maxY,
    parallel = diag_parallel_curvess(ctx),
    workers = ctx$args$workers
  )
  if (is.null(prof_pe)) {
    stop("load_profile_pe returned NULL.", call. = FALSE)
  }
  grid <- profile_grid_from_cfg(ctx$cfg)
  m_vals <- grid$M
  if (is.null(m_vals) || !length(m_vals)) {
    m_vals <- seq(0.5, 1.5, 0.25)
  }
  pr_opts <- production_plot_opts(ctx$cfg, "profile")
  p <- plot_production_functions(
    prof_pe,
    M_filter = m_vals,
    title = "Production function (M x SR_surv_Beta profile grid)",
    show_trajectories = pr_opts$show_trajectories,
    xlim = pr_opts$xlim,
    ylim = pr_opts$ylim,
    trajectory_yield_scale = pr_opts$trajectory_yield_scale
  )
  save_diag_figure(ctx, p, "production-functions-profile.png", width = 9, height = 6)
}

task_plot_kobe_scenarios <- function(ctx) {
  pe <- ctx_cache(ctx, "pe")
  if (is.null(pe)) {
    ss3 <- ctx_cache(ctx, "ss3")
    pe <- if (!is.null(ss3)) ss3$pe else NULL
  }
  if (is.null(pe)) {
    stop("PE data missing; run pe.load first.", call. = FALSE)
  }
  kb_opts <- kobe_plot_opts(ctx$cfg, "scenarios")
  p <- plot_kobe_phase_functions(
    pe,
    title = "Kobe phase plot by scenario",
    show_trajectories = kb_opts$show_trajectories,
    xlim = kb_opts$xlim,
    ylim = kb_opts$ylim,
    show_equilibrium = kb_opts$show_equilibrium %||% TRUE,
    facet_scales = kb_opts$facet_scales %||% "fixed",
    trajectory_linewidth = kb_opts$trajectory_linewidth,
    trajectory_colour = kb_opts$trajectory_colour
  )
  save_diag_figure(
    ctx,
    p,
    "kobe-phase-scenarios.png",
    width = kb_opts$width %||% 9,
    height = kb_opts$height %||% 6
  )
}

task_plot_kobe_profile <- function(ctx) {
  if (!has_capability(ctx$cfg, "likelihood_profiles")) {
    message("[plot] skipped Kobe profile (likelihood_profiles disabled for study)")
    return(invisible(NULL))
  }
  prof <- ctx_cache(ctx, "prof") %||% collect_likelihood_profiles(ctx$cfg, ctx$profile_scenario)
  if (is.null(prof)) {
    stop("Profile data missing; run profiles.collect or profiles.import first.", call. = FALSE)
  }
  ctx_cache(ctx, "prof", prof)
  prof_pe <- load_profile_pe(
    ctx$cfg,
    prof,
    maxY = ctx$args$maxY,
    parallel = diag_parallel_curvess(ctx),
    workers = ctx$args$workers
  )
  if (is.null(prof_pe)) {
    stop("load_profile_pe returned NULL.", call. = FALSE)
  }
  grid <- profile_grid_from_cfg(ctx$cfg)
  m_vals <- grid$M
  if (is.null(m_vals) || !length(m_vals)) {
    m_vals <- seq(0.5, 1.5, 0.25)
  }
  kb_opts <- kobe_plot_opts(ctx$cfg, "profile")
  p <- plot_kobe_phase_functions(
    prof_pe,
    M_filter = m_vals,
    title = "Kobe phase plot (M x SR_surv_Beta profile grid)",
    show_trajectories = kb_opts$show_trajectories,
    xlim = kb_opts$xlim,
    ylim = kb_opts$ylim,
    show_equilibrium = kb_opts$show_equilibrium %||% TRUE,
    facet_scales = kb_opts$facet_scales %||% "fixed",
    trajectory_linewidth = kb_opts$trajectory_linewidth,
    trajectory_colour = kb_opts$trajectory_colour
  )
  save_diag_figure(
    ctx,
    p,
    "kobe-phase-profile.png",
    width = kb_opts$width %||% 9,
    height = kb_opts$height %||% 6
  )
}

task_plot_ll_components <- function(ctx) {
  if (!has_capability(ctx$cfg, "likelihood_profiles")) {
    message("[plot] skipped likelihood profile components (likelihood_profiles disabled for study)")
    return(invisible(NULL))
  }
  prof <- ctx_cache(ctx, "prof") %||% collect_likelihood_profiles(ctx$cfg, ctx$profile_scenario)
  if (is.null(prof)) {
    stop("Profile data missing; run profiles.collect or profiles.import first.", call. = FALSE)
  }
  ctx_cache(ctx, "prof", prof)
  p <- plot_likelihood_profile_components(prof)
  save_diag_figure(ctx, p, "ll-profile-components.png", width = 10, height = 12)
}

task_plot_pe_four_panel <- function(ctx) {
  pe <- ctx_cache(ctx, "pe")
  if (is.null(pe)) {
    ss3 <- ctx_cache(ctx, "ss3")
    pe <- if (!is.null(ss3)) ss3$pe else NULL
  }
  if (is.null(pe)) {
    stop("PE data missing; run pe.load first.", call. = FALSE)
  }
  p <- plot_pe_four_panel(pe)
  save_diag_figure(ctx, p, "pe-four-panel.png", width = 9, height = 12)
}

diag_task_registry <- function() {
  list(
    profiles.setup = list(fn = task_profiles_setup),
    profiles.run = list(fn = task_profiles_run),
    profiles.collect = list(fn = task_profiles_collect),
    profiles.import = list(fn = task_profiles_import),
    pe.load = list(fn = task_pe_load),
    plot.production_scenarios = list(fn = task_plot_production_scenarios),
    plot.production_profile = list(fn = task_plot_production_profile),
    plot.kobe_scenarios = list(fn = task_plot_kobe_scenarios),
    plot.kobe_profile = list(fn = task_plot_kobe_profile),
    plot.ll_components = list(fn = task_plot_ll_components),
    plot.pe_four_panel = list(fn = task_plot_pe_four_panel)
  )
}

run_diag_task <- function(ctx, task_id) {
  task_id <- trimws(task_id)
  reg <- diag_task_registry()
  if (!task_id %in% names(reg)) {
    stop(
      "Unknown task: ", task_id,
      "\nAvailable: ", paste(names(reg), collapse = ", "),
      call. = FALSE
    )
  }
  t0 <- proc.time()[["elapsed"]]
  out <- reg[[task_id]]$fn(ctx)
  ctx$results[[paste0("task.", task_id)]] <- list(
    ok = TRUE,
    seconds = proc.time()[["elapsed"]] - t0
  )
  invisible(out)
}

run_diag_tasks <- function(ctx, tasks = NULL) {
  tasks <- tasks %||% resolve_diag_task_list(ctx$args)
  tasks <- unique(trimws(tasks))
  message("Diagnostics tasks: ", paste(tasks, collapse = ", "))
  message(
    "SS base: ", ctx$cfg$ss_base,
    "\nScenarios: ", paste(ctx$cfg$scenarios, collapse = ", "),
    "\nProfile scenario: ", ctx$profile_scenario
  )
  for (task_id in tasks) {
    run_diag_task(ctx, task_id)
  }
  figs <- ctx$results[vapply(ctx$results, is.character, logical(1))]
  if (length(figs)) {
    message("\nFigures:")
    for (nm in names(figs)) {
      message("  ", figs[[nm]])
    }
  }
  invisible(ctx)
}

diag_task_help <- function() {
  cat(
    "Diagnostic task presets (use with --tasks):\n",
    "  first_pass         setup + run + collect profiles, load PE, three figures\n",
    "  first_pass_import  import profiles + load PE + three figures\n",
    "  figures            plot only (needs existing SS3/profile outputs)\n",
    "  figures_scenarios  PE by scenario: production, Kobe, four-panel\n",
    "  process_error      alias for figures_scenarios\n",
    "  profiles           profile setup, run, collect only\n",
    "\nIndividual tasks:\n",
    paste0("  ", paste(names(diag_task_registry()), collapse = "\n  "), "\n"),
    sep = ""
  )
}
