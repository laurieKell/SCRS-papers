library(r4ss)

profileDir="C:/active/scrs-papers/ss-ll/data/ICCAT_2017_Model_Run_3_SSv330.23.02/00_Model_run_files/profile"
hVec     =seq(0.6, 1.0, by = 0.05)

# Edit the starter file as needed (see above)
starter <- SS_readstarter(file.path(mydir, "starter.ss"))
starter[["ctlfile"]] <- "control_modified.ss"
starter[["prior_like"]] <- 1
SS_writestarter(starter, dir = mydir, overwrite = TRUE)

# Run the profile
SS_profile(
  dir = mydir,
  masterctlfile = "control.ss_new",
  newctlfile = "control_modified.ss",
  string = "steep",
  profilevec = h.vec
)
