library(ggplot2)
library(ggpubr)
library(reshape2)

#' surv3Parm Stock-Recruitment Relationship
#'
#' @param S Spawning stock biomass (vector or scalar)
#' @param R0 unfished equilibrium recruitment
#' @param gamma parameter controlling shape, i.e curvature
#' @param beta Density-dependent parameter
#' @return Recruitment (R) for each S
surv3Parm<-function(S, R0, gamma, beta) {
  
  R=(R0*S)/(1+(S/beta)^gamma)
  
  return(R)}

dat=data.frame(S=seq(0,100)*5,
               R=surv3Parm(seq(0,100)*5, 250, 0.2, 1))
ggplot(dat)+geom_line(aes(S,R))
