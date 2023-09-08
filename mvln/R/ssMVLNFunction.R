# sskobe.R - DESC
# /home/mosqu003/sskobe.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

#' @examples
#' library(ss3diags)
#' 
#' data(sma)
#' sscor(sma$CoVar)
#'
#' hat(sma$derived_quants)
#'
#' ssmvln(sma$CoVar,sma$derived_quants)

## Returns stock & harvest as data.frame
sshat<-function(hat){
  names(hat)  =tolower(names(hat))
  
  y=rbind(hat[grep(paste("Recr","",  sep="_"),hat$label),],
          hat[grep(paste("Bratio","",sep="_"),hat$label),],
          hat[grep(paste("F",     "",sep="_"),hat$label),])
  y=y[substr(y$label,1,2)%in%c("F_","Br","Re"),]
  y=subset(y,!label%in%c("Recr_Initial","Recr_unfished","Recr_Virgin","Ret_Catch_MSY"))
  y=transform(y,cv2=(y$stddev/y$value)^2,
              var=log(1+(y$stddev/y$value)^2))[,-c(4,5)]
  
  names(y)=c("label","hat","stdLog","std","cv")
  y}

## Returns covariances as matrix
sscor<-function(covar){
  # UNIFY colnames
  setnames(covar, tolower(names(covar)))
  
  #
  flag <- unique(sort(c(
    grep(paste("Recr",  "",sep="_"),covar$label.i),
    grep(paste("Recr",  "",sep="_"),covar$label.j),
    grep(paste("Bratio","",sep="_"),covar$label.i),
    grep(paste("Bratio","",sep="_"),covar$label.j),
    grep(paste("F",     "",sep="_"),covar$label.i),
    grep(paste("F",     "",sep="_"),covar$label.j))))
  
  # GET correlations for all Fs, Bratios & Recrs
  cor <- covar[flag, c("label.i", "label.j", "corr")]
  
  flag <- substr(cor$label.i, 1, 2) %in% c("Re","F_", "Br") & 
    substr(cor$label.j, 1, 2) %in% c("Re","F_", "Br")
  cor <- cor[flag,]
  
  cor =rbind(cor,transform(cor,label.i=label.j,label.j=label.i))
  cor =subset(cor,!(is.na(label.i)|is.na(label.j)))
  
  cor =dcast(cor, label.i ~ label.j, value.var="corr")
  dmns=unname(unlist(cor[,1]))
  #  dimnames(cor)[[1]] <- dmns
  #  dmns=dimnames(cor)[[1]]
  
  cor=as.matrix(cor[,-1])
  cor[is.na(cor)]=0
  diag(cor)=1
  
  dimnames(cor)=list(dmns,dmns)
  cor}

# generates time series with serial correlations
ssmvln<-function(covar,hat,mc=5000,new=!FALSE){
  names(hat)  =tolower(names(hat))
  names(covar)=tolower(names(covar))
  
  #cov(x,y)
  cor=sscor(data.table(covar))
  
  if (mc<=1) return(cor)
  
  hat=data.table(hat)
  hat=subset(hat,label%in%dimnames(cor)[[1]])
  hat=sshat(hat)
  cor=cor[hat$label,hat$label]
  
  ## calculate Covariance matrix
  if(new)
    cvr=cor2cov(cor,hat$cv)
  else
    cvr=cor2cov(cor,hat$std^2)
  
  cvr[!is.finite(cvr)]=0
  
  mvlmu=log(hat$hat) 
  names(mvlmu)=dimnames(cor)[[1]]
  
  ## MC it all
  rtn=exp(mvtnorm::rmvnorm(mc, mean = mvlmu, sigma=cvr, method = c("svd")))
  
  ## return object
  rtn=data.table(rtn)
  names(rtn)=dimnames(cor)[[1]]
  
  dat=cbind(iter=seq(mc),rtn)
  dat=melt(dat,id="iter")
  dat=dat[, c("variable", "year"):= tstrsplit(variable, "_")]
  
  setnames(dat, c("variable", "value"), c("qname", "data"))
  return(dat)}


## FLQuants ####################################################################
## Reads in data.frame generates FLQuants ######################################

## Gets covariance between stock and harvest ###################################
covFLQ<-function(covar){
  dt=cbind(expand.grid(dimnames(covar)),value=c(covar))
  dt=subset(dt,substr(ac(Var1),8,nchar(ac(Var1)))==substr(ac(Var2),3,nchar(ac(Var2))))
  cov=cbind(dt,year=an(substr(ac(dt$Var2),3,nchar(ac(dt$Var2)))))
  as.FLQuant(data.frame(data=cov$value,year=cov$year))}

hatFLQ<-function(hat,data="hat"){
  r=subset(hat,substr(label,1,5)=="Recr_")
  r=transform(r,year=as.numeric(substr(label,6,nchar(label))))
  
  b=subset(hat,substr(label,1,7)=="Bratio_")
  b=transform(b,year=as.numeric(substr(label,8,nchar(label))))
  
  f=subset(hat,substr(label,1,2)=="F_")
  f=transform(f,year=as.numeric(substr(label,3,nchar(label))))
  
  rtn=rbind(cbind(qname="recruits",r),
            cbind(qname="stock",   b),
            cbind(qname="harvest", f))[,-2]
  rtn=cbind(rtn[,c("qname","year")],data=rtn[,data])
  
  as(rtn,"FLQuants")}

mvlnFLQ<-function(x){
  
  res=as(x, 'FLQuants')
  
  names(res)=c("recruits","stock","harvest")
  
  return(res)}
  
if(FALSE){
  library(calculus)
  library(compositions)
  
  mean =c(1.06, 0.76)
  sigma=matrix(c(1,-0.73,-0.73,1),2,2)
  
  fn<-function(x,y)
    dnorm.aplus(c(x,y),mean,sigma)
  integral(fn, bounds=list(x=c(1,0),y=c(Inf,1)))$value
  }