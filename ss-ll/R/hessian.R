# Install r4ss if needed
# install.packages("devtools")
# devtools::install_github("r4ss/r4ss")

library(r4ss)

# Read the main model output (for Hessian and parameter estimates)
ss3Dir="C:/active/scrs-papers/ss-ll/data/ICCAT_2017_Model_Run_3_SSv330.23.02/00_Model_run_files/1/1"
rep   =SS_output(dir=ss3Dir, verbose=FALSE)

# Read the profile output (assuming you used SS_profile or similar)
# If you used r4ss::SS_profile, it returns a list of outputs
profileDir ="C:/active/scrs-papers/ss-ll/data/ICCAT_2017_Model_Run_3_SSv330.23.02/00_Model_run_files/profile/output"
profile    =SSgetoutput(dirvec=list(profileDir))
profileSmry=SSsummarize(profile)

# Find the parameter of interest (e.g., "SR_LN(R0)")
paramName="SR_LN(R0)"

# Get the estimate and standard error from the base model
paramEst=rep$parameters[rep$parameters$Label==paramName, "Value"]
paramSe =rep$parameters[rep$parameters$Label==paramName, "StdDev"]

# Calculate the 95% CI from the Hessian
ciHessian=c(
  lower=param_est-1.96*param_se,
  upper=param_est+1.96*param_se)


# The profile summary gives likelihoods and parameter values
profileParam=profileSummary$pars[[paramName]]
profileLike =profileSummary$likelihoods[,"TOTAL"]
minLike     =min(profileLike,na.rm=TRUE)

# Calculate delta log-likelihood (relative to minimum)
deltaLike=profileLike - minLike
# Find parameter values within 1.92 log-likelihood units (approx 95% CI for 1 parameter)
profileCi=range(profileParam[deltaLike<1.92], na.rm=TRUE)


plot(
  profile_param, delta_like,
  type = "b", xlab = param_name, ylab = "Delta log-likelihood",
  main = paste("Likelihood profile vs. Hessian CI for", param_name)
)
abline(h = 1.92, col = "red", lty = 2)
abline(v = ci_hessian, col = "blue", lty = 3)
abline(v = param_est, col = "darkgreen", lty = 1)
abline(v = profile_ci, col = "purple", lty = 4)

legend("topright",
       legend = c("Hessian 95% CI", "Profile 95% CI", "MLE"),
       col = c("blue", "purple", "darkgreen"),
       lty = c(3, 4, 1)
)


cat("Hessian-based 95% CI:", round(ci_hessian, 3), "\n")
cat("Profile-based 95% CI:", round(profile_ci, 3), "\n")


If the blue (Hessian) and purple (profile) CIs are similar and centered on the MLE, the local and global uncertainty agree.

If the profile CI is much wider, asymmetric, or not centered on the MLE, the Hessian-based CI may be misleading. Trust the profile in such cases
