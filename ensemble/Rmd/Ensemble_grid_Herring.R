#########################################################################################################################
#                                                                                                                       #
#   SS3 ENSEMBLE MODEL script                                                                                           #
#   Run all scenarios at once and plot results, do diagnostic, do Ensemble on final weighted grid                       #
#                                                                                                                       #
#   September 2021                                                                                                      #
#                                                                                                                       #
#   Authors:                                                                                                            #
#   Francesco Masnadi (CNR-IRBIM & UNIBO, Ancona)                                                                       #
#   Massimiliano Cardinale (SLU Aqua, Lysekil)                                                                          #
#   Edits by David Gilljam (SLU Aqua, ?regrund)                                                                   #
#                                                                                                                       #
#  This script is now tailored on ensemble grid of 4 runs with 5 survey,15 Length data slot and 10 Age data slot.       #
#                                                                                                                       #
#  - ATTENTION: Age slot have to be activate in the Diagnostic loop before running the script (now only Length active)  #
#  - Also "sspar(mfrow=c(X,X)" in Diagnostic loop must be set accordingly                                               #
#  - Retro analysis is now set to -1 year to speed up simulation, please set to real value (at least 3 years)           #
#  - Mohn Rho value bound set now for long-live species, please change for short-live species                           #
#  - In Ensemble part make sure you set correct Fref depending on your choice (es. MSY,Btgt,..)                         #
#  - ATTENTION: Weighting vector automatically produce is a pure mean of all diagnostic scores as they are,             #
#    if the user want to use different Weighting methods (ex. merge multiple similar diagnostics into a single score)   #
#    have to produce the vector externally and save it in the main dir as "weight_vector.csv"                           #
#                                                                                                                       #
#                                                                                                                       #
#########################################################################################################################
rm(list=ls())
library(r4ss)
library(ss3diags)
library(tidyverse)
library(readr)
library(parallel)
library(doParallel)
#registerDoParallel(10)
#devtools::install_github("jabbamodel/ss3diags", force = TRUE)
#remotes::install_github("r4ss/r4ss")
cl <- makeCluster(10, type="PSOCK")
registerDoParallel(cl)

## See Ensemble_grid_Forecast_Herring.R for the derivation of RATIOLIM AND TRIG
RATIOLIM  = 2.333 # ratio of Btarget on Blim: e.g. if target=SSB30% and lim=SSB15% --> 35/15=2.333
RATIOTRIG = 1.667 # ratio of 1 on trigger percentage: e.g. if trigger=60% of target --> 1/0.6=1.667

STARTER_FILE  = "starter.ss"
WTATAGE_FILE  = "wtatage.ss"
PAR_FILE      = "ss.par"
FORECAST_FILE = "forecast.ss"
CTL_FILE      = "HerringSD2532.ctl"
DAT_FILE      = "HerringSD2532.dat"

runs <- paste0("Run", c(1:3)) # 1:N Number of run     
nmodels <- 3

MC.CORES = 32
MC.SILENT = TRUE

MVLN_MC       = 5000     ## Set to 5000 for final forecast MVLN runs and figures. 500 is fine for testing.

##
## Some system specific setups (Linux vs Windows)
##
if (Sys.info()["sysname"] == "Windows") {

    SS_EXE   <- "~/Max/Stock_synthesis/ss3_3.21/ss_win.exe"
    main.dir <- "~/Max/Commitees/ICES/WKBENCH/Central Baltic herring/Ensemble"
    MC.CORES <- 3
    print("mc.cores set to 1, NO PARALLELISATION on Windows :( ")

} else {

    SS_EXE <- "~/Max/Executives_SS/ss_linux"
    main.dir <- "~/Max/WKBENCH 2023/Central Baltic herring/Ensemble"
}

dir  <- main.dir
setwd(dir)

###############################################################################
## Run of all grid of models together and make SS plots (Skip if already done!)
###############################################################################
##for (i in 1:length(runs)){

##print(dir.runN <- paste0(dir,"/",runs[i]))
##flush.console()
##    Sys.sleep(3)

## XXX DG change:  mclapply() shouldn't be called within the loop
dir.runs <- file.path(paste0(dir, "/", runs))
print(sprintf("RUNNING MODELS IN PARALLEL: %s [mc.cores = %d, mc.silent = %s]",
              dir.runs, MC.CORES, as.character(MC.SILENT)));
flush.console()
Sys.sleep(5)

##This is the linux code for running in parallel
#mclapply(file.path(paste0(dir, "/", runs)),
#         run,
#         extras = "-nox", skipfinished = F, show_in_console = TRUE, exe=SS_EXE,
#         mc.cores = MC.CORES, mc.silent = MC.SILENT)

#This is the windoWs code for running in parallel
foreach(i=seq(1:3), .packages="r4ss") %dopar% {
  cat("core:", i + 1, "\n")
  run(dir.runs[[i]],extras = "-nox", skipfinished = F, show_in_console = TRUE, exe=SS_EXE)
  }

print("Model runs DONE, now generating SS standard plots, one at a time."); flush.console()
Sys.sleep(3)

for (i in 1:length(runs)){

    tmp <- SS_output(dir = file.path(runs[i]))
    SS_plots(tmp, uncertainty=T, datplot=F, forecast=T, maxyr=2025, minbthresh=0.25,fitrange = F)
}

###########################################################
# build the empty diags table for weighting the models later
############################################################
results <- data.frame(Run=runs, Convergence=NA, Total_LL=NA, N_Params=NA,
                      Runs_test_cpue1=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_cpue2=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_cpue3=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_cpue4=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_cpue5=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len1=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len2=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len3=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len4=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len5=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len6=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len7=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len8=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len9=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len10=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len11=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len12=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len13=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len14=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_len15=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age1=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age2=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age3=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age4=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age5=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age6=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age7=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age8=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age9=factor(NA, levels = c("Passed","Failed")),
                      Runs_test_age10=factor(NA, levels = c("Passed","Failed")),
                      RMSE_Perc=NA, RMSE_Perc_1=NA,RMSE_Perc_2=NA,
                      Retro_Rho_SSB= NA, Forecast_Rho_SSB=NA,
                      Retro_Rho_F= NA, Forecast_Rho_F=NA,
                      MASE_cpue1=NA, MASE_cpue2=NA, MASE_cpue3=NA, MASE_cpue4=NA, MASE_cpue5=NA,
                      MASE_len1=NA, MASE_len2=NA, MASE_len3=NA, MASE_len4=NA, MASE_len5=NA,
                      MASE_len6=NA, MASE_len7=NA,MASE_len8=NA,MASE_len9=NA,MASE_len10=NA,
                      MASE_len11=NA,MASE_len12=NA,MASE_len13=NA,MASE_len14=NA,MASE_len15=NA,
                      MASE_len16=NA,MASE_len17=NA,MASE_len18=NA,MASE_len19=NA,MASE_len20=NA,
                      MASE_age1=NA,MASE_age2=NA,MASE_age3=NA,MASE_age4=NA,MASE_age5=NA,
                      MASE_age6=NA,MASE_age7=NA,MASE_age8=NA,MASE_age9=NA,MASE_age10=NA
                      )

############################################################
## Make files for retro & change starter file
## Needed to allow retro analysis to run on multiple cores
############################################################

retro_files = paste0("retro", c(0,-1,-2,-3,-4,-5))

retroBatch <- list()

for(i in 1:length(runs)){
    dir.runN <- paste0(dir,"/",runs[i])
    dirname.Retrospective <- paste0(dir.runN,'/Retrospective')
    dir.create(path=dirname.Retrospective, showWarnings = TRUE, recursive = TRUE)

    for(j in 1:length(retro_files)){
        print("Creating retro directories and copying the nececcary SS config/data files.")

        dir.runN.new <- paste0(dir.runN,"/Retrospective/", retro_files[j])
        print(dir.runN.new)

        ## UPDATE THE RETROSPECTIVE BATCH LIST, FOR mclapply() below
        retroBatch <- c(retroBatch, as.list(dir.runN.new))

        dir.create(path=dir.runN.new, showWarnings = T, recursive = TRUE)

                                        # copy the SS base files in every retro subfolder
        file.copy(paste(dir.runN, "starter.ss", sep="/"),
                  paste(dir.runN.new, "starter.ss", sep="/"), overwrite=TRUE)

                                        # copy other files
        file.copy(paste(dir.runN, "HerringSD2532.ctl", sep="/"),
                  paste(dir.runN.new, "HerringSD2532.ctl", sep="/"), overwrite=TRUE)
        file.copy(paste(dir.runN, "HerringSD2532.dat", sep="/"),
                  paste(dir.runN.new, "HerringSD2532.dat", sep="/"), overwrite=TRUE)
        file.copy(paste(dir.runN, "forecast.ss", sep="/"),
                  paste(dir.runN.new, "forecast.ss", sep="/"), overwrite=TRUE)
        file.copy(paste(dir.runN, "wtatage.ss", sep="/"),
                  paste(dir.runN.new, "wtatage.ss", sep="/"), overwrite=TRUE)
        file.copy(paste(dir.runN, "wtatage.ss", sep="/"),
                  paste(dir.runN.new, "wtatage.ss", sep="/"), overwrite=TRUE)
        file.copy(paste(dir.runN, "ss.par", sep="/"),
                  paste(dir.runN.new, "ss.par", sep="/"), overwrite=TRUE)

                                        # Edit "starter.ss"
        starter.file <- readLines(paste(dir.runN.new, "/starter.ss", sep=""))

        linen <- NULL
        linen <- grep("# retrospective year relative to end year", starter.file)
        val_retro <- c(0,-1,-2,-3,-4,-5)
        starter.file[linen] <- paste0(val_retro[j] , " # retrospective year relative to end year" )

        write(starter.file, paste(dir.runN.new, "/starter.ss", sep=""))
    }
}

##
## ---- DG TEST TO SPEED UP RETRO RUNS. Seems to work nicely :) -----
##
##  No need to keep the retro sim runs in the diagnostics loop below
##
print(sprintf("RUNNING RETRO MODELS IN PARALLEL [mc.cores = %d, mc.silent = %s, total retros: %d]",
              MC.CORES, as.character(MC.SILENT), length(retroBatch)))
flush.console()
Sys.sleep(5)

##This is the linux code for running in parallel
#mclapply(retroBatch,
#         run,
#         extras = "-nox", skipfinished = F, show_in_console = TRUE, exe=SS_EXE,
#         mc.cores = MC.CORES, mc.silent = MC.SILENT)

##This is the windows code for running in parallel
n=6*nmodels
foreach(i=seq(1:n), .packages="r4ss") %dopar% {
  cat("core:", i + 1, "\n")
  run(retroBatch[[i]],extras = "-nox", skipfinished = F, show_in_console = TRUE, exe=SS_EXE)
}
#### ---------------------------------------------------------------------

##########################################################
##*********************************************************
## Diagnostic loop
##********************************************************
##########################################################

for(i in runs){
  # load all scenarios as list
    dir <- main.dir
    setwd(dir)
    readme<-paste0("/",i)
    tmp <- SS_output(dir=paste0(dir,readme),covar=T)

  ##############################
  ## check for Convergence
  ##############################
  ## put in final diags table
  results$Convergence[results$Run== i] <- tmp$maximum_gradient_component
  results$Total_LL[results$Run== i] <-   tmp$likelihoods_used$values[1]
  results$N_Params[results$Run== i] <-   tmp$N_estimated_parameters
  ## Par_nearBound <-as.data.frame(ifelse(tmp$parameters$Status  == "OK" | tmp$parameters$Status  == "act", 1, 0) %>% na.omit())

    ##*********************************************************
    ## Basic Residual Diagnostics (joint-residual and run test)
    ##********************************************************
# dir.create("Plotdiags",showWarnings = F)
    dir.runN <- paste0(dir,"/",i)
    dir.diag <- paste0(dir.runN,'/Plotdiags')
    dir.create(path=dir.diag, showWarnings = T, recursive = T)

    ## Joint-Residual (JABBAres)
    ## Check conflict between indices and mean length or age
    sspar(mfrow=c(1,2),plot.cex = 0.8)   # change to 3 if also Age present
    jr.cpue <- SSplotJABBAres(tmp,subplots="cpue",add=T,col=sscol(3)[c(1,3,2)])
    #jr.len <- SSplotJABBAres(tmp,subplots="len",add=T,col=sscol(3)[c(1,3,2)])
    jr.age <-SSplotJABBAres(tmp,subplots="age",add=T,col=sscol(3)[c(1,3,2)])  # activate if age data is present!
    dev.print(jpeg,paste0(dir = dir.diag,"/JointResiduals_",i,".jpg"), width = 8, height = 3.5, res = 300, units = "in")
    ## put in final diags table
    results$RMSE_Perc[results$Run== i] <-   jr.cpue$RMSE.perc[jr.cpue$indices=="Combined"]
    #results$RMSE_Perc_1[results$Run== i] <-   jr.len$RMSE.perc[jr.len$indices=="Combined"]
    results$RMSE_Perc_1[results$Run== i] <-   jr.age$RMSE.perc[jr.age$indices=="Combined"]

    ## Check Runs Test
    sspar(mfrow=c(3,2),plot.cex = 0.8)    # change based on number of fleet and length/age data
    rt.cpue <- SSplotRunstest(tmp,subplots="cpue",add=T, mixing="two.sided")
    #rt.len <-SSplotRunstest(tmp,subplots="len",add=T, mixing="two.sided")
    rt.age <- SSplotRunstest(tmp,subplots="age",add=T,mixing="two.sided")   # activate if age data is present!
    dev.print(jpeg,paste0(dir = dir.diag,"/RunTestResidual_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")
    ## put in final diags table
    ##CPUE
    results$Runs_test_cpue1[results$Run== i] <-   rt.cpue$test[1]
    results$Runs_test_cpue2[results$Run== i] <-   rt.cpue$test[2]
    results$Runs_test_cpue3[results$Run== i] <-   rt.cpue$test[3]
    results$Runs_test_cpue4[results$Run== i] <-   rt.cpue$test[4]
    results$Runs_test_cpue5[results$Run== i] <-   rt.cpue$test[5]
    ## length
    #results$Runs_test_len1[results$Run== i] <-   rt.len$test[1]
    #results$Runs_test_len2[results$Run== i] <-   rt.len$test[2]
    #results$Runs_test_len3[results$Run== i] <-   rt.len$test[3]
    #results$Runs_test_len4[results$Run== i] <-   rt.len$test[4]
    #results$Runs_test_len5[results$Run== i] <-   rt.len$test[5]
    #results$Runs_test_len6[results$Run== i] <-   rt.len$test[6]
    #results$Runs_test_len7[results$Run== i] <-   rt.len$test[7]
    #results$Runs_test_len8[results$Run== i] <-   rt.len$test[8]
    #results$Runs_test_len9[results$Run== i] <-   rt.len$test[9]
    #results$Runs_test_len10[results$Run== i] <-   rt.len$test[10]
    #results$Runs_test_len11[results$Run== i] <-   rt.len$test[11]
    #results$Runs_test_len12[results$Run== i] <-   rt.len$test[12]
    #results$Runs_test_len13[results$Run== i] <-   rt.len$test[13]
    #results$Runs_test_len14[results$Run== i] <-   rt.len$test[14]
    #results$Runs_test_len15[results$Run== i] <-   rt.len$test[15]
  # Age
   results$Runs_test_age1[results$Run== i] <-   rt.age$test[1]
   results$Runs_test_age2[results$Run== i] <-   rt.age$test[2]
   results$Runs_test_age3[results$Run== i] <-   rt.age$test[3]
   results$Runs_test_age4[results$Run== i] <-   rt.age$test[4]
   results$Runs_test_age5[results$Run== i] <-   rt.age$test[5]
   results$Runs_test_age6[results$Run== i] <-   rt.age$test[6]
   results$Runs_test_age7[results$Run== i] <-   rt.age$test[7]
   results$Runs_test_age8[results$Run== i] <-   rt.age$test[8]
   results$Runs_test_age9[results$Run== i] <-   rt.age$test[9]
   results$Runs_test_age10[results$Run== i] <-   rt.age$test[10]

  ##############################
  # Retrospective analyses
##############################

    ##
    ## DG note: To utilise the full set of cores, we should do parallel runs of
    ##          all scenario runs' retro-0, then all runs' retro-1 etc.
    ##          Then we would run 27 simulations in parallel, 6 times,
    ##          instead of 6 simulations in parallel, 27 times.
    ##
    ## DG note2: Move the actual running of the retrospective sims out from this diagnostics
    ##           loop, to the previous loop where the directories and files are set up.
    ##           Then, we could run the diagnostics on finished runs, without having
    ##           to do the restrospective simulations.
    ##
    ##    FIXED! SEE RETRO BATCH RUNS BEFORE THE DIAGNOSTICS LOOP!
    ##

    ## Automatically running the retrospective analyses
    start.retro <- 0    #end year of model
    end.retro   <- 5    #number of years for retrospective

    ## Retro simulations ----------------------------------------------------------------------------
    ## run retro on 6 cores
    dir.runN <- paste0(dir,"/", i, "/Retrospective")
    if (0) {
        ## RETRO SIMULATIONS ARE NOW RUN BEFORE THIS DIAGNOSTICS LOOP!
        mclapply(file.path(paste0(dir.runN, "/", retro_files)),
                 run, extras = "-nohess", skipfinished = T, exe=SS_EXE, show_in_console = TRUE,mc.cores=MC.CORES)
    }
    ##-----------------------------------------------------------------------------------------------

    ## Read "SS_doRetro" output
    ##
    ## XXX DG: dirname.Retrospective is not updated in this for-loop, the last Run-dir will be
    ##         repeatedly used!
    ##
    ## DG Fix:
    retroModels <- SSgetoutput(dirvec = file.path(paste0(dir.runN),
    ## OLD CODE: retroModels <- SSgetoutput(dirvec = file.path(paste0(dirname.Retrospective),
                                                  paste("retro", start.retro:-end.retro, sep="")))

    ## Summarize output
    retroSummary <- r4ss::SSsummarize(retroModels)
    endyrvec <- retroSummary$endyrs + start.retro:-end.retro

    sspar(mfrow=c(1,2),plot.cex = 0.8)
    retro.ssb <- SSplotRetro(retroSummary,add=T,legendcex=0.8,tickEndYr=F,xylabs=F,legendloc = "bottomleft",uncertainty = T,showrho = T,forecast = T,labels="SSB (t)", endyrvec = c(endyrvec), subplots = c("SSB"))
    retro.f <- SSplotRetro(retroSummary,add=T,legendcex=0.8,tickEndYr=F,xylabs=F,legendloc = "bottomleft",uncertainty = T,showrho = T,forecast = T,labels="SSB (t)", endyrvec = c(endyrvec), subplots = c("F"))
    dev.print(jpeg,paste0(dir = dir.diag,"/Retro_",i,".jpg"), width = 12, height = 5, res = 300, units = "in")

    ## put in final diags table
    results$Retro_Rho_SSB[results$Run== i] <-   retro.ssb$Rho[retro.ssb$peel=="Combined"]
    results$Forecast_Rho_SSB[results$Run== i] <-   retro.ssb$ForecastRho[retro.ssb$peel=="Combined"]
    results$Retro_Rho_F[results$Run== i] <-   retro.f$Rho[retro.f$peel=="Combined"]
    results$Forecast_Rho_F[results$Run== i] <-   retro.f$ForecastRho[retro.f$peel=="Combined"]

  ##############################
  ## Hindcasting analyses
  ##############################
    ## Do Hindcast with Cross-Validation of CPUE observations
    sspar(mfrow=c(2,1),plot.cex = 0.9)     # change based on number of cpue or survey
    mase_cpue <- SSplotHCxval(retroSummary,add=T)
    mase_cpueB <- SSmase(retroSummary,quants = "cpue")
    dev.print(jpeg,paste0(dir = dir.diag,"/HCxval_CPUE_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")
    results$MASE_cpue1[results$Run== i] <-   mase_cpueB$MASE[1]
    results$MASE_cpue2[results$Run== i] <-   mase_cpueB$MASE[2]
    results$MASE_cpue3[results$Run== i] <-   mase_cpueB$MASE[3]
    results$MASE_cpue4[results$Run== i] <-   mase_cpueB$MASE[4]
    results$MASE_cpue5[results$Run== i] <-   mase_cpueB$MASE[5]

    ## Also Hindcast with Cross-Validation for mean length or age
    ## Use new converter fuction SSretroComps()
    hccomps = SSretroComps(retroModels)
    ## Specify subplots = "age" or "len" in SSplotHCxval
    sspar(mfrow=c(2,2),plot.cex = 0.7) # change based on number of fleet or survey per length or age data
    mase_len.plot <- SSplotHCxval(hccomps,add=T,subplots = "age",legendloc="topright",legend = FALSE, indexUncertainty = T,legendcex = 1)
    mase_len.plot2 <- SSplotHCxval(hccomps,add=T,subplots = "age",legendloc="topright",legend = FALSE, indexUncertainty = T,legendcex = 1, Season = 4)
#  mase_age <- SSplotHCxval(hccomps,add=T,subplots = "age",legendloc="topright",legend = FALSE, indexUncertainty = T,legendcex = 1, Season = 4) # activate if age data is present!
    dev.print(jpeg,paste0(dir = dir.diag,"/HCxval_age_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")

    ## put in final diags table
    #mase_len <- SSmase(hccomps,quants = "len")
    #mase_len2 <- SSmase(hccomps,quants = "len", Season = c(4))
    mase_age <- SSmase(hccomps,quants = "age")
    mase_age2 <- SSmase(hccomps,quants = "age")
    #results$MASE_len1[results$Run== i] <-   mase_len$MASE.adj[1]
    #results$MASE_len2[results$Run== i] <-   mase_len$MASE.adj[2]
    #results$MASE_len3[results$Run== i] <-   mase_len2$MASE.adj[3]
    #results$MASE_len4[results$Run== i] <-   mase_len$MASE.adj[4]
    #results$MASE_len5[results$Run== i] <-   mase_len$MASE.adj[5]
    #results$MASE_len6[results$Run== i] <-   mase_len$MASE.adj[6]
    #results$MASE_len7[results$Run== i] <-   mase_len$MASE.adj[7]
    #results$MASE_len8[results$Run== i] <-   mase_len$MASE.adj[8]
    #results$MASE_len9[results$Run== i] <-   mase_len$MASE.adj[9]
    #results$MASE_len10[results$Run== i] <-   mase_len$MASE.adj[10]
    #results$MASE_len11[results$Run== i] <-   mase_len$MASE.adj[11]
    #results$MASE_len12[results$Run== i] <-   mase_len$MASE.adj[12]
    #results$MASE_len13[results$Run== i] <-   mase_len$MASE.adj[13]
    #results$MASE_len14[results$Run== i] <-   mase_len$MASE.adj[14]
    #results$MASE_len15[results$Run== i] <-   mase_len$MASE.adj[15]
    #results$MASE_len16[results$Run== i] <-   mase_len$MASE.adj[16]
    #results$MASE_len17[results$Run== i] <-   mase_len$MASE.adj[17]
    #results$MASE_len18[results$Run== i] <-   mase_len$MASE.adj[18]
    #results$MASE_len19[results$Run== i] <-   mase_len$MASE.adj[19]
    #results$MASE_len20[results$Run== i] <-   mase_len$MASE.adj[20]
  # age
    results$MASE_age1[results$Run== i] <-   mase_age$MASE.adj[1]
    results$MASE_age2[results$Run== i] <-   mase_age$MASE.adj[2]
    results$MASE_age3[results$Run== i] <-   mase_age2$MASE.adj[3]
    results$MASE_age4[results$Run== i] <-   mase_age$MASE.adj[4]
    results$MASE_age5[results$Run== i] <-   mase_age$MASE.adj[5]
    results$MASE_age6[results$Run== i] <-   mase_age$MASE.adj[6]
    results$MASE_age7[results$Run== i] <-   mase_age$MASE.adj[7]
    results$MASE_age8[results$Run== i] <-   mase_age$MASE.adj[8]
    results$MASE_age9[results$Run== i] <-   mase_age$MASE.adj[9]
    results$MASE_age10[results$Run== i] <-   mase_age$MASE.adj[10]

} ## Diagnostics/retrospective loop end

##*********************************************************
## save the table of diagnostic with real value
##********************************************************
##diags_table <- as.data.frame(t(results) %>% na.omit())
tmp <- t(results)
diags_table <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]
write.csv(diags_table, file = paste0(main.dir,"/Diags_table.csv"))

#*********************************************************
# save the table of diagnostic for weigthing purpose (value 0 to 1)
#********************************************************
results2 <- results %>% dplyr::select(-Convergence,-Total_LL,-N_Params,-Run)

results2$Runs_test_cpue1 <-   ifelse(results2$Runs_test_cpue1 == 'Passed', 1, 0)
results2$Runs_test_cpue2 <-   ifelse(results2$Runs_test_cpue2 == 'Passed', 1, 0)
results2$Runs_test_cpue3 <-   ifelse(results2$Runs_test_cpue3 == 'Passed', 1, 0)
results2$Runs_test_cpue4 <-   ifelse(results2$Runs_test_cpue4 == 'Passed', 1, 0)
results2$Runs_test_cpue5 <-   ifelse(results2$Runs_test_cpue5 == 'Passed', 1, 0)

# length
#results2$Runs_test_len1 <-   ifelse(results2$Runs_test_len1 == 'Passed', 1, 0)
#results2$Runs_test_len2 <-   ifelse(results2$Runs_test_len2 == 'Passed', 1, 0)
#results2$Runs_test_len3 <-   ifelse(results2$Runs_test_len3 == 'Passed', 1, 0)
#results2$Runs_test_len4 <-   ifelse(results2$Runs_test_len4 == 'Passed', 1, 0)
#results2$Runs_test_len5 <-   ifelse(results2$Runs_test_len5 == 'Passed', 1, 0)
#results2$Runs_test_len6 <-   ifelse(results2$Runs_test_len6 == 'Passed', 1, 0)
#results2$Runs_test_len7 <-   ifelse(results2$Runs_test_len7 == 'Passed', 1, 0)
#results2$Runs_test_len8 <-   ifelse(results2$Runs_test_len8 == 'Passed', 1, 0)
#results2$Runs_test_len9 <-   ifelse(results2$Runs_test_len9 == 'Passed', 1, 0)
#results2$Runs_test_len10 <-   ifelse(results2$Runs_test_len10 == 'Passed', 1, 0)
#results2$Runs_test_len11 <-   ifelse(results2$Runs_test_len11 == 'Passed', 1, 0)
#results2$Runs_test_len12 <-   ifelse(results2$Runs_test_len12 == 'Passed', 1, 0)
#results2$Runs_test_len13 <-   ifelse(results2$Runs_test_len13 == 'Passed', 1, 0)
#results2$Runs_test_len14 <-   ifelse(results2$Runs_test_len14 == 'Passed', 1, 0)
#results2$Runs_test_len15 <-   ifelse(results2$Runs_test_len15 == 'Passed', 1, 0)

# Age
results2$Runs_test_age1 <-   ifelse(results2$Runs_test_age1 == 'Passed', 1, 0)
results2$Runs_test_age2 <-   ifelse(results2$Runs_test_age2 == 'Passed', 1, 0)
results2$Runs_test_age3 <-   ifelse(results2$Runs_test_age3 == 'Passed', 1, 0)
results2$Runs_test_age4 <-   ifelse(results2$Runs_test_age4 == 'Passed', 1, 0)
results2$Runs_test_age5 <-   ifelse(results2$Runs_test_age5 == 'Passed', 1, 0)
results2$Runs_test_age6 <-   ifelse(results2$Runs_test_age6 == 'Passed', 1, 0)
results2$Runs_test_age7 <-   ifelse(results2$Runs_test_age7 == 'Passed', 1, 0)
results2$Runs_test_age8 <-   ifelse(results2$Runs_test_age8 == 'Passed', 1, 0)
results2$Runs_test_age9 <-   ifelse(results2$Runs_test_age9 == 'Passed', 1, 0)
results2$Runs_test_age10 <-   ifelse(results2$Runs_test_age10 == 'Passed', 1, 0)

results2$RMSE_Perc <- ifelse(results2$RMSE_Perc < 30, 1, 0)
results2$RMSE_Perc_1 <- ifelse(results2$RMSE_Perc_1 < 30, 1, 0)

# Age
#  results2$RMSE_Perc_2 <- ifelse(results2$RMSE_Perc_2 < 30, 1, 0)
results2$Retro_Rho_SSB <- ifelse(results2$Retro_Rho_SSB > -0.15 & results2$Retro_Rho_SSB < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Forecast_Rho_SSB <- ifelse(results2$Forecast_Rho_SSB > -0.15 & results2$Forecast_Rho_SSB < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Retro_Rho_F <- ifelse(results2$Retro_Rho_F > -0.15 & results2$Retro_Rho_F < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Forecast_Rho_F <- ifelse(results2$Forecast_Rho_F > -0.15 & results2$Forecast_Rho_F < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$MASE_cpue1 <- ifelse(results2$MASE_cpue1 < 1, 1, 0)
results2$MASE_cpue2 <- ifelse(results2$MASE_cpue2 < 1, 1, 0)
results2$MASE_cpue3 <- ifelse(results2$MASE_cpue3 < 1, 1, 0)
results2$MASE_cpue4 <- ifelse(results2$MASE_cpue4 < 1, 1, 0)
results2$MASE_cpue5 <- ifelse(results2$MASE_cpue5 < 1, 1, 0)
#results2$MASE_len1 <- ifelse(results2$MASE_len1 < 1, 1, 0)
#results2$MASE_len2 <- ifelse(results2$MASE_len2 < 1, 1, 0)
#results2$MASE_len3 <- ifelse(results2$MASE_len3 < 1, 1, 0)
#results2$MASE_len4 <- ifelse(results2$MASE_len4 < 1, 1, 0)
#results2$MASE_len5 <- ifelse(results2$MASE_len5 < 1, 1, 0)
#results2$MASE_len6 <- ifelse(results2$MASE_len6 < 1, 1, 0)
#results2$MASE_len7 <- ifelse(results2$MASE_len7 < 1, 1, 0)
#results2$MASE_len8 <- ifelse(results2$MASE_len8 < 1, 1, 0)
#results2$MASE_len9 <- ifelse(results2$MASE_len9 < 1, 1, 0)
#results2$MASE_len10 <- ifelse(results2$MASE_len10 < 1, 1, 0)
#results2$MASE_len11 <- ifelse(results2$MASE_len11 < 1, 1, 0)
#results2$MASE_len12 <- ifelse(results2$MASE_len12 < 1, 1, 0)
#results2$MASE_len13 <- ifelse(results2$MASE_len13 < 1, 1, 0)
#results2$MASE_len14 <- ifelse(results2$MASE_len14 < 1, 1, 0)
#results2$MASE_len15 <- ifelse(results2$MASE_len15 < 1, 1, 0)
#results2$MASE_len16 <- ifelse(results2$MASE_len16 < 1, 1, 0)
#results2$MASE_len17 <- ifelse(results2$MASE_len17 < 1, 1, 0)
#results2$MASE_len18 <- ifelse(results2$MASE_len18 < 1, 1, 0)
#results2$MASE_len19 <- ifelse(results2$MASE_len19 < 1, 1, 0)
#results2$MASE_len20 <- ifelse(results2$MASE_len20 < 1, 1, 0)

# Age
 results2$MASE_age1 <- ifelse(results2$MASE_age1 < 1, 1, 0)
 results2$MASE_age2 <- ifelse(results2$MASE_age2 < 1, 1, 0)
 results2$MASE_age3 <- ifelse(results2$MASE_age3 < 1, 1, 0)
 results2$MASE_age4 <- ifelse(results2$MASE_age4 < 1, 1, 0)
 results2$MASE_age5 <- ifelse(results2$MASE_age5 < 1, 1, 0)
 results2$MASE_age6 <- ifelse(results2$MASE_age6 < 1, 1, 0)
 results2$MASE_age7 <- ifelse(results2$MASE_age7 < 1, 1, 0)
 results2$MASE_age8 <- ifelse(results2$MASE_age8 < 1, 1, 0)
 results2$MASE_age9 <- ifelse(results2$MASE_age9 < 1, 1, 0)
 results2$MASE_age10 <- ifelse(results2$MASE_age10 < 1, 1, 0)

weight_table <-  as.data.frame(t(results2) %>% na.omit())
##tmp <- t(results2)
##weight_table <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]
write.csv(weight_table, file = paste0(main.dir,"/Weight_table.csv"))

####################################################
#### weight_vector to be used in Ensemble
####################################################
weight_vector <- rowMeans(as.data.frame(lapply(as.data.frame(t(weight_table)), as.numeric)))
write.csv(t(weight_vector), file = paste0(main.dir,"/weight_vector.csv"),row.names=FALSE)

##****************************************************************
## Approximate uncertainty with MVLN (hessian)
##****************************************************************
setwd(main.dir)  # set dir to the initial one

#### weight_vector from diagnostic scores coming from "Ensemble_diags_runs.R" script
weight_vector <- unname(unlist(read_csv(paste0(main.dir,"/weight_vector.csv"))[1,]))

kbproj = NULL

# Compile MVLN posteriors by scenario run
for(i in 1:length(runs)){
    ## load all scenarios as list
    run = SS_output(dir=file.path(runs[i]))
    ## get MVLN mvn.temp for each scenario
    ## Make sure you set correct Fref depending on your choice
    ##mvn.temp = SSdeltaMVLN_DG(run, run = runs[i],
    mvn.temp = SSdeltaMVLN_CF(run, run = runs[i],
                           years = run$startyr:run$endyr+1,
                           ##years = run$startyr:run$endyr+1,
                           addprj = TRUE, mc = MVLN_MC, weight = weight_vector[i], Fref = "Btgt", plot=FALSE, verbose=TRUE)
    kbproj = rbind(kbproj,
                   data.frame(mvn.temp$kb, RR0 = mvn.temp$kb$Recr/mvn.temp$refpts[mvn.temp$refpts$RefPoint=="R0","value"], model=runs[i]))
    ## save labels once
    if(i==1) labels = mvn.temp$labels
}

#### save Ensemble r.data
save(kbproj, file="Ensemble_model_CF.rdata")
load(file="Ensemble_model.rdata")

##########################
## FINAL Kobe_plot ensemble
##########################
sspar(mfrow=c(1,1),plot.cex = 0.7)
SSplotKobe(kbproj, fill=T, joint=F, posterior="points", ci.levels = c(0.5, 0.95, 0.99), ylim = c(0,6.5),
           ylab=expression(F/F[trg]),
           xlab=expression(SSB/SSB[trg]), legendpos = "right", legend = TRUE)
dev.print(jpeg,paste0("Kobe_final_points_cv.jpg"), width = 12, height = 8, res = 300, units = "in")

##
## (1) Show trajectories one by one. One figure with all years, one with the last 20 years
## (2) Show trajectories All together. One figure with all years, one with the last 20 years
##
##years   <- c(run$startyr, run$endyr)
years   <- c(run$startyr, run$endyr+1)
indRuns <- kbproj$run
quants  <-  c("stock", "harvest", "SSB", "F", "Recr", "Catch")
labels  <- expression(SSB/SSB[trg], "F/F"[SB ~ 35], "SSB", "F", "Recruits", "Catch")

for (i in 1:2) {

    ## (1)
    sspar(mfrow=c(3, 2), plot.cex = 0.9)
    ##SSplotEnsemble(kbproj, add=T, legendcex = 0.3, legendloc="topright", legendncol = 3, ylabs = labels, xlim = years, lwd=1)
    ##SSplotEnsemble(kbproj, add=T, subplots = quants[1], legendcex = 0.3, legendloc="topright", legendncol = 3, ylabs = labels[1], xlim = years, lwd=1)
    SSplotEnsemble(kbproj, add=T, subplots = quants[1], legendcex = 0.3, legendloc="topright", legendncol = 3, ylabs = labels[1], xlim = years, lwd=1)
    abline(h = 1/RATIOLIM, lty = 1)   # Add SSB/BLim horisontal line
    SSplotEnsemble(kbproj, add=T, subplots = quants[2:6], legendcex = 0.3, legendloc="topright", legendncol = 3, ylabs = labels[2:6], xlim = years, lwd=1)

    dev.print(jpeg, sprintf("MLVN_Compare.jpg", paste(years, collapse="-")),
              width = 12, height = 8, res = 300, units = "in")

    ## (2)
    kbproj$run <- "All_runs"
    sspar(mfrow=c(3,2), plot.cex = 0.9)
    ##SSplotEnsemble(kbproj, models="all", add=T, subplots=quants[1], legendcex = 0.4,legendloc="topleft",legendncol = 2, ylabs=labels[1], xlim = years)
    SSplotEnsemble(kbproj, models="all", add=T, subplots=quants[1], legendcex = 0.4,legendloc="topleft",legendncol = 2, ylabs=labels[1], xlim = years)
    abline(h = 1/RATIOLIM, lty = 1)   # Add SSB/BLim horisontal line
    SSplotEnsemble(kbproj, models="all", add=T, subplots=quants[2:6], legendcex = 0.4,legendloc="topleft",legendncol = 2, ylabs=labels[2:6], xlim = years)
    dev.print(jpeg, sprintf("MLVN_All.jpg", paste(years, collapse="-")),
              width = 12, height = 8, res = 300, units = "in")

    ## Set years and reset for next loop
    kbproj$run <- indRuns
    years      <- c(run$startyr, run$endyr+1)
}

## get median + 90% CIs for SSB/SSBmsy across all models
prj_SBB_SBBtrg = aggregate(stock~year, kbproj , quantile, c(0.5,0.05,0.95))
prj_F_Btrg = aggregate(harvest~year, kbproj, quantile,c(0.5,0.05,0.95))
prj_SSB = aggregate(SSB~year, kbproj, quantile, c(0.5,0.05,0.95))
prj_F = aggregate(F~year,kbproj, quantile,c(0.5,0.05,0.95))
prj_Recr = aggregate(Recr~year, kbproj, quantile,c(0.5,0.05,0.95))
prj_RR0 = aggregate(RR0~year, kbproj, quantile,c(0.5,0.05,0.95))
write.csv(prj_SBB_SBBtrg, "Ensemble_SBB_SBBtrg.csv")
write.csv(prj_F_Btrg, "Ensemble_F_Btrg.csv")
write.csv(prj_SSB, "Ensemble_SBB.csv")
write.csv(prj_F, "Ensemble_F.csv")
write.csv(prj_Recr, "Ensemble_Recr.csv")
write.csv(prj_RR0, "Ensemble_RR0.csv")

save.image("ensemble_image.Rdata")

kbprojR0 <- kbproj[kbproj$year==years[[1]],]



