#EXPERIMENT 1A: NON-RANDOMIZED EXPERIMENT WITH X-SPECIFIC COMMON TREND, TIME-INVARIANT COVARIATES AND HOMOGENEOUS EFFECTS
#              PROPENSITY SCORE CORRECTLY SPECIFIED, OUTCOME REGRESSION CORRECTLY SPECIFIED


# clean the R environment
rm(list = ls())

repository="C:/Users/tommy/OneDrive/Desktop/tesi/DID simulation/MC simulation"
setwd(repository)
file.name="EXP_1A"

# install.packages("devtools")
library(devtools)
#install_github("clu0/diffindiff")
# devtools::install_github("pedrohcgs/DRDID")
library(diffindiff)

library(readstata13)
library(DRDID)
library(glmnet)
library(pacman)
library(data.table)
library(fixest)
library(stargazer)
library(dplyr)
library(magrittr)
library(sandwich)
library(miceadds)
library(ncpen)
library(GGally)
library(MASS)
library(tidyverse)
library(fabricatr)
library(JWileymisc)
library(rlearner)
library(EQL)
library(ggplot2)
library(readstata13)
library(miceadds)
library(randomForest)
library(party)

#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################

#create triplematching function for later use

tripleIPWRA <- function(id, t, d, X, Y, method='logit'){
  
  # triple DiD propensity score weighting following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X matrix of covariates
  # t is the time dummy
  # d is the treatment group dummy
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  # method selects the estimation method for the propensity score. The available
  #estimation methods are:"logit", "probit", "lasso" and "randomforest".
  
  
  
  
  # create four groups and prepare dataset
  T1=ifelse(t==1& d==1,1,0)
  T0=ifelse(t==0& d==1,1,0)
  C1=ifelse(t==1& d==0,1,0)
  C0=ifelse(t==0& d==0,1,0)
  
  data=as.data.frame(cbind(Y,id, t, d, C0, C1,T0, T1, X))
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=X[T1==1|C1==1,]
  
  #Estimate propensity score
  if (method=='logit'){
    mylogit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                   data = data1)
    data1$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                    data = data1)
    data1$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov1 <- poly(as.matrix(data1[,9:ncol(data1)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov1, data1$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov1, type='response')
    
  } else if (method=='randomforest'){
    data1rf=data1[,8:ncol(data1)] 
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X1))))
    fr.fit=cforest(T1~ .,data=data1rf, controls=mycontrols)
    data1$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data1$pscore[data1$pscore>1]=1
  data1$pscore[data1$pscore<0]=0
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=X[T1==1|C0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                   data = data2)
    data2$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                    data = data2)
    data2$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression
    # perform LASSO linear regression with 10-fold cross validation
    cov2 <- poly(as.matrix(data2[,9:ncol(data2)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov2, data2$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov2, type="response")
    
  } else if (method=='randomforest'){
    data2rf=data2[,8:ncol(data2)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X2))))
    fr.fit=cforest(T1~ .,data=data2rf, controls=mycontrols)
    data2$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data2$pscore[data2$pscore>1]=1
  data2$pscore[data2$pscore<0]=0
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=X[T1==1|T0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                   data = data3)
    data3$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                    data = data3)
    data3$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # Estimate propensity score with LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov3 <- poly(as.matrix(data3[,9:ncol(data3)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov3, data3$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov3, type="response")
    
  } else if (method=='randomforest'){
    
    #Estimate propensity score with random forest
    data3rf=data3[,8:ncol(data3)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X3))))
    fr.fit=cforest(T1~ .,data=data3rf, controls=mycontrols)
    data3$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  
  #Check that the p-value is in between 0 and 0
  data3$pscore[data3$pscore>1]=1
  data3$pscore[data3$pscore<0]=0
  
  
  #Merge propensity score
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  #trimming for pscore value close to 0 and to 1
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  
  #Preparing data for regression
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  Xmatrix=as.matrix(cs_X)
  
  #Computing weights
  cs_data$w_att=rep(1,cs_len)
  cs_data$w_att=ifelse(cs_data$T1==0, cs_data$pscore/(1-cs_data$pscore),1)
  
  
  #Regression with propensity score weights
  e_tripleIPWRA=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix
                                      +I(Xmatrix*t)+I(Xmatrix*d),weights=cs_data$w_att,cluster="d")
  
  return(e_tripleIPWRA)
  
}


###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Repeated Cross Section Data
###################################################################################
# Pre = T stands for pre-treatment period
# treat = F  stands for control group

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = F){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))
  
  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }
  
  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))
  
  # return fitted values
  return(list(out.reg = out.delta))
  
}

# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)
  
  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]
  
  if (!any(is.na(gam))) {
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw)
    
    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean)
    
    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}
#Loss function for estimation of the Bias reduced PS, based on Graham, Pinton and Egel (2012, 2016)

loss.ps.IPT <- function(gamma1, n, D, int.cov, iw){
  #Coefficients for quadratic extrapolation
  cn <- -(n - 1)
  bn <- -n + (n - 1) * log(n - 1)
  an <- -(n - 1) * (1 - log(n - 1) + 0.5 * (log(n - 1))^2)
  vstar <- log(n - 1)
  
  v <- gamma1 %*% t(int.cov)
  phi <- ifelse(v < vstar, - v - exp(v), an + bn * v + 0.5 * cn * (v^2))
  phi1 <- ifelse(v < vstar, - 1 - exp(v), bn + cn * v)
  phi2 <- ifelse(v < vstar, - exp(v), cn)
  
  #phi <- (v<vstar) * (- v - exp(v)) + (v>=vstar) * (an + bn * v + 0.5 * cn * (v^2))
  #phi1 <- (v<vstar) * (- 1 - exp(v)) + (v>=vstar) * (bn  + cn * v)
  #phi2 <- (v<vstar) * (- exp(v)) + (v>=vstar) * cn
  
  # Minus is because nlm minimizes functions, and we aim to maximize!
  res <- - sum(iw * (1 - D) * phi + v)
  
  attr(res, "gradient") <- - t(int.cov) %*% as.vector(iw * ((1-D) * phi1 + 1))
  attr(res, "hessian")  <-  - t(as.vector((1-D) * iw * phi2) * int.cov) %*% int.cov
  return(res)
}

###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  pslogit <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial",
                                         weights = i.weights))
  
  init.gamma <- suppressWarnings(stats::coef(pslogit))
  
  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))
  
  flag <- ifelse(pscore.cal$converged, 0, 1)
  
  gamma.cal <- try(pscore.cal$argument)
  
  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {
    
    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)
    
    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }
  
  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))
  
  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")
    
    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }
    
  }
  
  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))
  
}


lasso_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type =  "weighted", nboot = NULL,
                          inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates), degree=3, raw=TRUE)
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}


randforest_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                               boot = FALSE, boot.type =  "weighted", nboot = NULL,
                               inffunc = FALSE){
  #----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by RANDOM FOREST
  data=as.data.frame(cbind(D, int.cov))
  mycontrols <- cforest_unbiased(ntree=100, as.integer(sqrt(ncol(int.cov))))
  fr.fit=cforest(D~ .,data=data, controls=mycontrols)
  ps.fit=as.numeric(predict(fr.fit, type='prob'))
  
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using RANDOM FOREST.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  reg.cont.coeff.pre=cforest(y_C0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C0*ipw_C0)
  out.y.cont.pre=as.numeric(predict(reg.cont.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  reg.cont.coeff.post <- cforest(y_C1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C1*ipw_C1)
  out.y.cont.post=as.numeric(predict(reg.cont.coeff.post, newdata=data))
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  reg.treat.coeff.pre <- cforest(y_T0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T0)
  
  out.y.treat.pre=as.numeric(predict(reg.treat.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  reg.treat.coeff.post <- cforest(y_T1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T1)
  out.y.treat.post=as.numeric(predict(reg.treat.coeff.post, newdata=data))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}

triple_drdid_rc <- function(y, post, D,id, int.cov, i.weights = NULL,
                            boot = FALSE, boot.type =  "weighted", nboot = NULL,
                            inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore using triplematching tecnique
  T0=ifelse(dta_long$post==0 & dta_long$d==1, 1, 0)
  T1=ifelse(dta_long$post==1 & dta_long$d==1, 1, 0)
  C0=ifelse(dta_long$post==0 & dta_long$d==0, 1, 0)
  C1=ifelse(dta_long$post==1 & dta_long$d==0, 1, 0)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, post, d, covariates, y))
  
  len=nrow(data)
  
  
  #First propensity score matching between T1 and C1
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(covariates, T1==1 | C1==1)
  X1=as.matrix(X1)
  myprobit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                  data = data1)
  
  data1$pscore=myprobit$fitted.values
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(covariates, T1==1 | C0==1)
  X2=as.matrix(X2)
  myprobit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                  data = data2)
  
  data2$pscore=myprobit$fitted.values
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(covariates, T1==1 | T0==1)
  X3=as.matrix(X3)
  myprobit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                  data = data3)
  
  data3$pscore=myprobit$fitted.values
  
  
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  # create dataset with common support and trimming for pscore value close to 0 and to 1
  
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  cs_len=nrow(cs_data)
  
  
  cs_data$w_att=rep(1,cs_len)
  
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  D=cs_data$d
  y=cs_data$y
  cs_int.cov <- as.matrix(rep(1,nrow(cs_data)))
  int.cov=as.matrix(cbind(cs_int.cov,cs_data[,8:11]))
  post=cs_data$post
  pscore.ipt<- cs_data$pscore
  i.weights=i.weights[cs_data$id]
  ps.fit <- as.vector(pscore.ipt)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.post <-  as.vector(out.y.post$out.reg)
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
  
  
}

#create rep.row function for later use

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
# START OF THE SIMULATION
################################################################################

## Parameters and seed

set.seed(1)      # Seed
n = 1000         # Sample size
M = 500         # Number of experiments/iterations




# Results Storage 

did <- rep(0,M)
did2 <- rep(0,M)
did3 <- rep(0,M)
ipwdid<- rep(0,M)
ordid<- rep(0,M)
impdrdid<- rep(0,M)
lasso_impdrdid<- rep(0,M)
randfor_impdrdid<- rep(0,M)
trIPWRA<- rep(0,M)
lasso_trIPWRA<- rep(0,M)
randfor_trIPWRA<- rep(0,M)
trIPWRA2<- rep(0,M)
lasso_trIPWRA2<- rep(0,M)
randfor_trIPWRA2<- rep(0,M)
trweightRA<- rep(0,M)
trweightRA2<- rep(0,M)
did_time <- rep(0,M)
did2_time <- rep(0,M)
did3_time <- rep(0,M)
ipwdid_time<- rep(0,M)
ordid_time<- rep(0,M)
impdrdid_time<- rep(0,M)
lasso_impdrdid_time<- rep(0,M)
randfor_impdrdid_time<- rep(0,M)
trIPWRA_time<- rep(0,M)
lasso_trIPWRA_time<- rep(0,M)
randfor_trIPWRA_time<- rep(0,M)
trIPWRA2_time<- rep(0,M)
lasso_trIPWRA2_time<- rep(0,M)
randfor_trIPWRA2_time<- rep(0,M)
trweightRA_time<- rep(0,M)
trweightRA2_time<- rep(0,M)
ATT<- rep(0,M)

# START OF THE SIMULATION LOOP

for (i in 1:M){ #  M is the number of iterations
  
  
  # Sample size
  n <- 1000
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Proportion in each period
  lambda <- 0.5
  # NUmber of bootstrapped draws
  nboot <- 199
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
  
  # Create heterogenenous effects for the ATT, which is set approximately equal to zero
  index.unobs.het <- d * (index.lin)
  index.att <- 0
  
  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  #Gen realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)
  
  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects
  
  # Gen realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  
  # Generate "T"
  post <- as.numeric(stats::runif(n) <= lambda)
  
  # observed outcome
  y <- post * y1 + (1 - post) * y0
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a long data frame
  dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- dta_long[order(dta_long$id),]  
  
  #Standard standard TWFE
  start.time <- Sys.time()
  did_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did[i] <- did_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did2_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did2_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did2[i] <- did2_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did3_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post
               +x1*d+x2*d+x3*d+x4*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did3_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did3[i] <- did3_i$coefficients["post:d"]
  
  # Implement DID ipw Abadie(2005) with normalized weights
  start.time <- Sys.time()
  ipwdid_i=ipwdid(yname="y", tname = "post", idname = "id",  dname = "d",
                  xformla= ~ x1 + x2 + x3 + x4,
                  data = dta_long, panel = FALSE,
                  boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ipwdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ipwdid[i] <- ipwdid_i$ATT
  
  # Implement Outcome regression model following Heckman, Ichimura and Todd (1997) but with linear assumption of covariates
  
  start.time <- Sys.time()
  ordid_i=ordid(yname="y", tname = "post", idname = "id",  dname = "d",
                xformla= ~ x1 + x2 + x3 + x4,
                data = dta_long, panel = FALSE,
                boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ordid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ordid[i] <- ordid_i$ATT
  
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020):
  start.time <- Sys.time()
  impdrdid_i <- drdid(yname="y", tname = "post", idname = "id",  dname = "d",
                      xformla= ~ x1 + x2 + x3 + x4,
                      data = dta_long, panel = FALSE, estMethod = "imp")
  end.time <-Sys.time()
  tk <- end.time - start.time
  impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  impdrdid[i] <- impdrdid_i$ATT
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with LASSO:
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  boot = FALSE
  boot.type =  "weighted"
  nboot = NULL
  inffunc = FALSE
  
  start.time <- Sys.time()
  lasso_impdrdid[i]=lasso_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                   boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with RANDOM FOREST:
  start.time <- Sys.time()
  randfor_impdrdid[i]=randforest_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                          boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  
  # Implement logit tripleIPWRA building on Bludell(2004)
  
  X=(cbind(dta_long$x1,dta_long$x2,dta_long$x3,dta_long$x4))
  id=dta_long$id
  t=dta_long$post
  d=dta_long$d
  Y=dta_long$y
  
  start.time <- Sys.time()
  trIPWRA_i=tripleIPWRA(id, t, d, X, Y)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trIPWRA[i] <- coef(trIPWRA_i)["t:d"]
  
  # Implement LASSO tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  lasso_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='lasso')
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  lasso_trIPWRA[i] <- coef(lasso_trIPWRA_i)["t:d"]
  
  # Implement random forest tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  randfor_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='randomforest')
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  randfor_trIPWRA[i] <- coef(randfor_trIPWRA_i)["t:d"]
  
  # Implement logit tripleDRDiD building by on Sant'Anna and Zhao (2020) 
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  id=dta_long$id
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  
  
  start.time <- Sys.time()
  trweightRA_i=triple_drdid_rc(y=y, post=post, D=D, id, covariates, i.weights = NULL, boot = FALSE,
                               boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trweightRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trweightRA[i] <- trweightRA_i
  
  
  print(i)
}


#generate summary statistics for each estimator

did_mean=mean(did)
did_variance=mean((did-mean(did))^2)
did_rmse=sqrt(mean((did-ATT)^2))
did_bias=mean (mean(did)-ATT)
did_avtime=mean (did_time)

did2_mean=mean(did2)
did2_variance=mean((did2-mean(did2))^2)
did2_rmse=sqrt(mean((did2-ATT)^2))
did2_bias=mean (mean(did2)-ATT)
did2_avtime=mean (did2_time)

did3_mean=mean(did3)
did3_variance=mean((did3-mean(did3))^2)
did3_rmse=sqrt(mean((did3-ATT)^2))
did3_bias=mean (mean(did3)-ATT)
did3_avtime=mean (did3_time)

ipwdid_mean=mean(ipwdid)
ipwdid_variance=mean((ipwdid-mean(ipwdid))^2)
ipwdid_rmse=sqrt(mean((ipwdid-ATT)^2))
ipwdid_bias=mean (mean(ipwdid)-ATT)
ipwdid_avtime=mean (ipwdid_time)

ordid_mean=mean(ordid)
ordid_variance=mean((ordid-mean(ordid))^2)
ordid_rmse=sqrt(mean((ordid-ATT)^2))
ordid_bias=mean (mean(ordid)-ATT)
ordid_avtime=mean (ordid_time)

impdrdid_mean=mean(impdrdid)
impdrdid_variance=mean((impdrdid-mean(impdrdid))^2)
impdrdid_rmse=sqrt(mean((impdrdid-ATT)^2))
impdrdid_bias=mean (mean(impdrdid)-ATT)
impdrdid_avtime=mean (impdrdid_time)

lasso_impdrdid_mean=mean(lasso_impdrdid)
lasso_impdrdid_variance=mean((lasso_impdrdid-mean(lasso_impdrdid))^2)
lasso_impdrdid_rmse=sqrt(mean((lasso_impdrdid-ATT)^2))
lasso_impdrdid_bias=mean (mean(lasso_impdrdid)-ATT)
lasso_impdrdid_avtime=mean (lasso_impdrdid_time)

randfor_impdrdid_mean=mean(randfor_impdrdid)
randfor_impdrdid_variance=mean((randfor_impdrdid-mean(randfor_impdrdid))^2)
randfor_impdrdid_rmse=sqrt(mean((randfor_impdrdid-ATT)^2))
randfor_impdrdid_bias=mean (mean(randfor_impdrdid)-ATT)
randfor_impdrdid_avtime=mean (randfor_impdrdid_time)

trIPWRA_mean=mean(trIPWRA)
trIPWRA_variance=mean((trIPWRA-mean(trIPWRA))^2)
trIPWRA_rmse=sqrt(mean((trIPWRA-ATT)^2))
trIPWRA_bias=mean (mean(trIPWRA)-ATT)
trIPWRA_avtime=mean (trIPWRA_time)

lasso_trIPWRA_mean=mean(lasso_trIPWRA)
lasso_trIPWRA_variance=mean((lasso_trIPWRA-mean(lasso_trIPWRA))^2)
lasso_trIPWRA_rmse=sqrt(mean((lasso_trIPWRA-ATT)^2))
lasso_trIPWRA_bias=mean (mean(lasso_trIPWRA)-ATT)
lasso_trIPWRA_avtime=mean (lasso_trIPWRA_time)

randfor_trIPWRA_mean=mean(randfor_trIPWRA)
randfor_trIPWRA_variance=mean((randfor_trIPWRA-mean(randfor_trIPWRA))^2)
randfor_trIPWRA_rmse=sqrt(mean((randfor_trIPWRA-ATT)^2))
randfor_trIPWRA_bias=mean (mean(randfor_trIPWRA)-ATT)
randfor_trIPWRA_avtime=mean (randfor_trIPWRA_time)

trweightRA_mean=mean(trweightRA)
trweightRA_variance=mean((trweightRA-mean(trweightRA))^2)
trweightRA_rmse=sqrt(mean((trweightRA-ATT)^2))
trweightRA_bias=mean (mean(trweightRA)-ATT)
trweightRA_avtime=mean (trweightRA_time)


#construct a summary table

mean=cbind(did_mean,did2_mean,did3_mean, ipwdid_mean, ordid_mean, impdrdid_mean, lasso_impdrdid_mean,
           randfor_impdrdid_mean, trIPWRA_mean,lasso_trIPWRA_mean, randfor_trIPWRA_mean,
           trweightRA_mean)
bias=cbind(did_bias,did2_bias,did3_bias, ipwdid_bias, ordid_bias, impdrdid_bias, lasso_impdrdid_bias,
           randfor_impdrdid_bias, trIPWRA_bias,lasso_trIPWRA_bias, randfor_trIPWRA_bias,
           trweightRA_bias)
rmse=cbind(did_rmse,did2_rmse,did3_rmse, ipwdid_rmse, ordid_rmse, impdrdid_rmse, lasso_impdrdid_rmse,
           randfor_impdrdid_rmse, trIPWRA_rmse,lasso_trIPWRA_rmse, randfor_trIPWRA_rmse,
           trweightRA_rmse)
variance=cbind(did_variance,did2_variance,did3_variance, ipwdid_variance, ordid_variance, impdrdid_variance, lasso_impdrdid_variance,
               randfor_impdrdid_variance, trIPWRA_variance,lasso_trIPWRA_variance, randfor_trIPWRA_variance,
               trweightRA_variance)
time=cbind(did_avtime,did2_avtime,did3_avtime, ipwdid_avtime, ordid_avtime, impdrdid_avtime, lasso_impdrdid_avtime,
           randfor_impdrdid_avtime, trIPWRA_avtime,lasso_trIPWRA_avtime, randfor_trIPWRA_avtime,
           trweightRA_avtime)

tab <- matrix(1:48, ncol=4, byrow=TRUE)
tab[,1]=bias
tab[,2]=rmse
tab[,3]=variance
tab[,4]=time
tab=round(tab, digits = 3)


colnames(tab) <- c('Bias','RMSE', 'Variance','Time')
rownames(tab) <- c('TWFE','TWFE (T*X)','TWFE (T*X+D*X)','IPW', 'RA','DRDiD', 'LASSO DRDiD', 'RF DRDiD',
                   '3IPWRA', 'LASSO 3IPWRA', 'RF 3IPWRA','3WDRDiD')
latextable=stargazer(tab)
tab
#save table
write.table(tab, file=paste(file.name,'.txt',sep = ""))

#Create plot to visualize results
tabdata=as.data.frame(tab)
tabdata$name=row.names(tabdata)
ggplot(data=tabdata, aes(x=reorder(name, -abs(Bias)), y=abs(Bias))) + 
  geom_bar(stat = "identity")+
  scale_y_continuous(limits=c(0, 33))+
  coord_flip()+
  geom_point(aes(y=RMSE),
             stat="identity",
             position="dodge",
             alpha=1,
             shape=21,
             stroke = 1.8,
             size=3,
             colour = "white",
             fill = "black")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())
ggsave("EXP1A.png")

#commands to read table and transfer to latex
#stargazer(read.table(paste('EXP0B','.txt',sep = "")), summary=FALSE)

#EXPERIMENT 1B: NON-RANDOMIZED EXPERIMENT WITH X-SPECIFIC COMMON TREND, TIME-INVARIANT COVARIATES AND HOMOGENEOUS EFFECTS
#              PROPENSITY SCORE NOT CORRECTLY SPECIFIED, OUTCOME REGRESSION CORRECTLY SPECIFIED


# clean the R environment
rm(list = ls())

repository="C:/Users/tommy/OneDrive/Desktop/tesi/DID simulation/MC simulation"
setwd(repository)
file.name="EXP_1B"

# install.packages("devtools")
library(devtools)
#install_github("clu0/diffindiff")
# devtools::install_github("pedrohcgs/DRDID")
library(diffindiff)

library(readstata13)
library(DRDID)
library(glmnet)
library(pacman)
library(data.table)
library(fixest)
library(stargazer)
library(dplyr)
library(magrittr)
library(sandwich)
library(miceadds)
library(ncpen)
library(GGally)
library(MASS)
library(tidyverse)
library(fabricatr)
library(JWileymisc)
library(rlearner)
library(EQL)
library(ggplot2)
library(readstata13)
library(miceadds)
library(randomForest)
library(party)

#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################

#create triplematching function for later use

tripleIPWRA <- function(id, t, d, X, Y, method='logit'){
  
  # triple DiD propensity score weighting following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X matrix of covariates
  # t is the time dummy
  # d is the treatment group dummy
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  # method selects the estimation method for the propensity score. The available
  #estimation methods are:"logit", "probit", "lasso" and "randomforest".
  
  
  
  
  # create four groups and prepare dataset
  T1=ifelse(t==1& d==1,1,0)
  T0=ifelse(t==0& d==1,1,0)
  C1=ifelse(t==1& d==0,1,0)
  C0=ifelse(t==0& d==0,1,0)
  
  data=as.data.frame(cbind(Y,id, t, d, C0, C1,T0, T1, X))
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=X[T1==1|C1==1,]
  
  #Estimate propensity score
  if (method=='logit'){
    mylogit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                   data = data1)
    data1$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                    data = data1)
    data1$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov1 <- poly(as.matrix(data1[,9:ncol(data1)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov1, data1$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov1, type='response')
    
  } else if (method=='randomforest'){
    data1rf=data1[,8:ncol(data1)] 
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X1))))
    fr.fit=cforest(T1~ .,data=data1rf, controls=mycontrols)
    data1$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data1$pscore[data1$pscore>1]=1
  data1$pscore[data1$pscore<0]=0
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=X[T1==1|C0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                   data = data2)
    data2$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                    data = data2)
    data2$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression
    # perform LASSO linear regression with 10-fold cross validation
    cov2 <- poly(as.matrix(data2[,9:ncol(data2)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov2, data2$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov2, type="response")
    
  } else if (method=='randomforest'){
    data2rf=data2[,8:ncol(data2)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X2))))
    fr.fit=cforest(T1~ .,data=data2rf, controls=mycontrols)
    data2$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data2$pscore[data2$pscore>1]=1
  data2$pscore[data2$pscore<0]=0
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=X[T1==1|T0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                   data = data3)
    data3$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                    data = data3)
    data3$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # Estimate propensity score with LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov3 <- poly(as.matrix(data3[,9:ncol(data3)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov3, data3$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov3, type="response")
    
  } else if (method=='randomforest'){
    
    #Estimate propensity score with random forest
    data3rf=data3[,8:ncol(data3)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X3))))
    fr.fit=cforest(T1~ .,data=data3rf, controls=mycontrols)
    data3$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  
  #Check that the p-value is in between 0 and 0
  data3$pscore[data3$pscore>1]=1
  data3$pscore[data3$pscore<0]=0
  
  
  #Merge propensity score
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  #trimming for pscore value close to 0 and to 1
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  
  #Preparing data for regression
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  Xmatrix=as.matrix(cs_X)
  
  #Computing weights
  cs_data$w_att=rep(1,cs_len)
  cs_data$w_att=ifelse(cs_data$T1==0, cs_data$pscore/(1-cs_data$pscore),1)
  
  
  #Regression with propensity score weights
  e_tripleIPWRA=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix
                                      +I(Xmatrix*t)+I(Xmatrix*d),weights=cs_data$w_att,cluster="d")
  
  return(e_tripleIPWRA)
  
}


###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Repeated Cross Section Data
###################################################################################
# Pre = T stands for pre-treatment period
# treat = F  stands for control group

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = F){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))
  
  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }
  
  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))
  
  # return fitted values
  return(list(out.reg = out.delta))
  
}

# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)
  
  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]
  
  if (!any(is.na(gam))) {
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw)
    
    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean)
    
    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}
#Loss function for estimation of the Bias reduced PS, based on Graham, Pinton and Egel (2012, 2016)

loss.ps.IPT <- function(gamma1, n, D, int.cov, iw){
  #Coefficients for quadratic extrapolation
  cn <- -(n - 1)
  bn <- -n + (n - 1) * log(n - 1)
  an <- -(n - 1) * (1 - log(n - 1) + 0.5 * (log(n - 1))^2)
  vstar <- log(n - 1)
  
  v <- gamma1 %*% t(int.cov)
  phi <- ifelse(v < vstar, - v - exp(v), an + bn * v + 0.5 * cn * (v^2))
  phi1 <- ifelse(v < vstar, - 1 - exp(v), bn + cn * v)
  phi2 <- ifelse(v < vstar, - exp(v), cn)
  
  #phi <- (v<vstar) * (- v - exp(v)) + (v>=vstar) * (an + bn * v + 0.5 * cn * (v^2))
  #phi1 <- (v<vstar) * (- 1 - exp(v)) + (v>=vstar) * (bn  + cn * v)
  #phi2 <- (v<vstar) * (- exp(v)) + (v>=vstar) * cn
  
  # Minus is because nlm minimizes functions, and we aim to maximize!
  res <- - sum(iw * (1 - D) * phi + v)
  
  attr(res, "gradient") <- - t(int.cov) %*% as.vector(iw * ((1-D) * phi1 + 1))
  attr(res, "hessian")  <-  - t(as.vector((1-D) * iw * phi2) * int.cov) %*% int.cov
  return(res)
}

###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  pslogit <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial",
                                         weights = i.weights))
  
  init.gamma <- suppressWarnings(stats::coef(pslogit))
  
  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))
  
  flag <- ifelse(pscore.cal$converged, 0, 1)
  
  gamma.cal <- try(pscore.cal$argument)
  
  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {
    
    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)
    
    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }
  
  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))
  
  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")
    
    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }
    
  }
  
  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))
  
}


lasso_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type =  "weighted", nboot = NULL,
                          inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates), degree=3, raw=TRUE)
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}


randforest_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                               boot = FALSE, boot.type =  "weighted", nboot = NULL,
                               inffunc = FALSE){
  #----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by RANDOM FOREST
  data=as.data.frame(cbind(D, int.cov))
  mycontrols <- cforest_unbiased(ntree=100, as.integer(sqrt(ncol(int.cov))))
  fr.fit=cforest(D~ .,data=data, controls=mycontrols)
  ps.fit=as.numeric(predict(fr.fit, type='prob'))
  
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using RANDOM FOREST.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  reg.cont.coeff.pre=cforest(y_C0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C0*ipw_C0)
  out.y.cont.pre=as.numeric(predict(reg.cont.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  reg.cont.coeff.post <- cforest(y_C1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C1*ipw_C1)
  out.y.cont.post=as.numeric(predict(reg.cont.coeff.post, newdata=data))
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  reg.treat.coeff.pre <- cforest(y_T0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T0)
  
  out.y.treat.pre=as.numeric(predict(reg.treat.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  reg.treat.coeff.post <- cforest(y_T1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T1)
  out.y.treat.post=as.numeric(predict(reg.treat.coeff.post, newdata=data))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}

triple_drdid_rc <- function(y, post, D,id, int.cov, i.weights = NULL,
                            boot = FALSE, boot.type =  "weighted", nboot = NULL,
                            inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore using triplematching tecnique
  T0=ifelse(dta_long$post==0 & dta_long$d==1, 1, 0)
  T1=ifelse(dta_long$post==1 & dta_long$d==1, 1, 0)
  C0=ifelse(dta_long$post==0 & dta_long$d==0, 1, 0)
  C1=ifelse(dta_long$post==1 & dta_long$d==0, 1, 0)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, post, d, covariates, y))
  
  len=nrow(data)
  
  
  #First propensity score matching between T1 and C1
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(covariates, T1==1 | C1==1)
  X1=as.matrix(X1)
  myprobit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                  data = data1)
  
  data1$pscore=myprobit$fitted.values
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(covariates, T1==1 | C0==1)
  X2=as.matrix(X2)
  myprobit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                  data = data2)
  
  data2$pscore=myprobit$fitted.values
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(covariates, T1==1 | T0==1)
  X3=as.matrix(X3)
  myprobit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                  data = data3)
  
  data3$pscore=myprobit$fitted.values
  
  
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  # create dataset with common support and trimming for pscore value close to 0 and to 1
  
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  cs_len=nrow(cs_data)
  
  
  cs_data$w_att=rep(1,cs_len)
  
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  D=cs_data$d
  y=cs_data$y
  cs_int.cov <- as.matrix(rep(1,nrow(cs_data)))
  int.cov=as.matrix(cbind(cs_int.cov,cs_data[,8:11]))
  post=cs_data$post
  pscore.ipt<- cs_data$pscore
  i.weights=i.weights[cs_data$id]
  ps.fit <- as.vector(pscore.ipt)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.post <-  as.vector(out.y.post$out.reg)
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
  
  
}

#create rep.row function for later use

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
# START OF THE SIMULATION
################################################################################

## Parameters and seed

set.seed(1)      # Seed
n = 1000         # Sample size
M = 500         # Number of experiments/iterations




# Results Storage 

did <- rep(0,M)
did2 <- rep(0,M)
did3 <- rep(0,M)
ipwdid<- rep(0,M)
ordid<- rep(0,M)
impdrdid<- rep(0,M)
lasso_impdrdid<- rep(0,M)
randfor_impdrdid<- rep(0,M)
trIPWRA<- rep(0,M)
lasso_trIPWRA<- rep(0,M)
randfor_trIPWRA<- rep(0,M)
trIPWRA2<- rep(0,M)
lasso_trIPWRA2<- rep(0,M)
randfor_trIPWRA2<- rep(0,M)
trweightRA<- rep(0,M)
trweightRA2<- rep(0,M)
did_time <- rep(0,M)
did2_time <- rep(0,M)
did3_time <- rep(0,M)
ipwdid_time<- rep(0,M)
ordid_time<- rep(0,M)
impdrdid_time<- rep(0,M)
lasso_impdrdid_time<- rep(0,M)
randfor_impdrdid_time<- rep(0,M)
trIPWRA_time<- rep(0,M)
lasso_trIPWRA_time<- rep(0,M)
randfor_trIPWRA_time<- rep(0,M)
trIPWRA2_time<- rep(0,M)
lasso_trIPWRA2_time<- rep(0,M)
randfor_trIPWRA2_time<- rep(0,M)
trweightRA_time<- rep(0,M)
trweightRA2_time<- rep(0,M)
ATT<- rep(0,M)

# START OF THE SIMULATION LOOP

for (i in 1:M){ #  M is the number of iterations
  
  
  # Sample size
  n <- 1000
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Proportion in each period
  lambda <- 0.5
  # NUmber of bootstrapped draws
  nboot <- 199
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
  
  # Create heterogenenous effects for the ATT, which is set approximately equal to zero
  index.unobs.het <- d * (index.lin)
  index.att <- 0
  
  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  #Gen realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)
  
  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects
  
  # Gen realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  
  # Generate "T"
  post <- as.numeric(stats::runif(n) <= lambda)
  
  # observed outcome
  y <- post * y1 + (1 - post) * y0
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a long data frame
  dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- dta_long[order(dta_long$id),]  
  
  #Standard standard TWFE
  start.time <- Sys.time()
  did_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did[i] <- did_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did2_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did2_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did2[i] <- did2_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did3_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post
               +x1*d+x2*d+x3*d+x4*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did3_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did3[i] <- did3_i$coefficients["post:d"]
  
  # Implement DID ipw Abadie(2005) with normalized weights
  start.time <- Sys.time()
  ipwdid_i=ipwdid(yname="y", tname = "post", idname = "id",  dname = "d",
                  xformla= ~ x1 + x2 + x3 + x4,
                  data = dta_long, panel = FALSE,
                  boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ipwdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ipwdid[i] <- ipwdid_i$ATT
  
  # Implement Outcome regression model following Heckman, Ichimura and Todd (1997) but with linear assumption of covariates
  
  start.time <- Sys.time()
  ordid_i=ordid(yname="y", tname = "post", idname = "id",  dname = "d",
                xformla= ~ x1 + x2 + x3 + x4,
                data = dta_long, panel = FALSE,
                boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ordid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ordid[i] <- ordid_i$ATT
  
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020):
  start.time <- Sys.time()
  impdrdid_i <- drdid(yname="y", tname = "post", idname = "id",  dname = "d",
                      xformla= ~ x1 + x2 + x3 + x4,
                      data = dta_long, panel = FALSE, estMethod = "imp")
  end.time <-Sys.time()
  tk <- end.time - start.time
  impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  impdrdid[i] <- impdrdid_i$ATT
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with LASSO:
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  boot = FALSE
  boot.type =  "weighted"
  nboot = NULL
  inffunc = FALSE
  
  start.time <- Sys.time()
  lasso_impdrdid[i]=lasso_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                   boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with RANDOM FOREST:
  start.time <- Sys.time()
  randfor_impdrdid[i]=randforest_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                          boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  
  # Implement logit tripleIPWRA building on Bludell(2004)
  
  X=(cbind(dta_long$x1,dta_long$x2,dta_long$x3,dta_long$x4))
  id=dta_long$id
  t=dta_long$post
  d=dta_long$d
  Y=dta_long$y
  
  start.time <- Sys.time()
  trIPWRA_i=tripleIPWRA(id, t, d, X, Y)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trIPWRA[i] <- coef(trIPWRA_i)["t:d"]
  
  # Implement LASSO tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  lasso_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='lasso')
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  lasso_trIPWRA[i] <- coef(lasso_trIPWRA_i)["t:d"]
  
  # Implement random forest tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  randfor_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='randomforest')
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  randfor_trIPWRA[i] <- coef(randfor_trIPWRA_i)["t:d"]
  
  # Implement logit tripleDRDiD building by on Sant'Anna and Zhao (2020) 
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  id=dta_long$id
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  
  
  start.time <- Sys.time()
  trweightRA_i=triple_drdid_rc(y=y, post=post, D=D, id, covariates, i.weights = NULL, boot = FALSE,
                               boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trweightRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trweightRA[i] <- trweightRA_i
  
  
  print(i)
}


#generate summary statistics for each estimator

did_mean=mean(did)
did_variance=mean((did-mean(did))^2)
did_rmse=sqrt(mean((did-ATT)^2))
did_bias=mean (mean(did)-ATT)
did_avtime=mean (did_time)

did2_mean=mean(did2)
did2_variance=mean((did2-mean(did2))^2)
did2_rmse=sqrt(mean((did2-ATT)^2))
did2_bias=mean (mean(did2)-ATT)
did2_avtime=mean (did2_time)

did3_mean=mean(did3)
did3_variance=mean((did3-mean(did3))^2)
did3_rmse=sqrt(mean((did3-ATT)^2))
did3_bias=mean (mean(did3)-ATT)
did3_avtime=mean (did3_time)

ipwdid_mean=mean(ipwdid)
ipwdid_variance=mean((ipwdid-mean(ipwdid))^2)
ipwdid_rmse=sqrt(mean((ipwdid-ATT)^2))
ipwdid_bias=mean (mean(ipwdid)-ATT)
ipwdid_avtime=mean (ipwdid_time)

ordid_mean=mean(ordid)
ordid_variance=mean((ordid-mean(ordid))^2)
ordid_rmse=sqrt(mean((ordid-ATT)^2))
ordid_bias=mean (mean(ordid)-ATT)
ordid_avtime=mean (ordid_time)

impdrdid_mean=mean(impdrdid)
impdrdid_variance=mean((impdrdid-mean(impdrdid))^2)
impdrdid_rmse=sqrt(mean((impdrdid-ATT)^2))
impdrdid_bias=mean (mean(impdrdid)-ATT)
impdrdid_avtime=mean (impdrdid_time)

lasso_impdrdid_mean=mean(lasso_impdrdid)
lasso_impdrdid_variance=mean((lasso_impdrdid-mean(lasso_impdrdid))^2)
lasso_impdrdid_rmse=sqrt(mean((lasso_impdrdid-ATT)^2))
lasso_impdrdid_bias=mean (mean(lasso_impdrdid)-ATT)
lasso_impdrdid_avtime=mean (lasso_impdrdid_time)

randfor_impdrdid_mean=mean(randfor_impdrdid)
randfor_impdrdid_variance=mean((randfor_impdrdid-mean(randfor_impdrdid))^2)
randfor_impdrdid_rmse=sqrt(mean((randfor_impdrdid-ATT)^2))
randfor_impdrdid_bias=mean (mean(randfor_impdrdid)-ATT)
randfor_impdrdid_avtime=mean (randfor_impdrdid_time)

trIPWRA_mean=mean(trIPWRA)
trIPWRA_variance=mean((trIPWRA-mean(trIPWRA))^2)
trIPWRA_rmse=sqrt(mean((trIPWRA-ATT)^2))
trIPWRA_bias=mean (mean(trIPWRA)-ATT)
trIPWRA_avtime=mean (trIPWRA_time)

lasso_trIPWRA_mean=mean(lasso_trIPWRA)
lasso_trIPWRA_variance=mean((lasso_trIPWRA-mean(lasso_trIPWRA))^2)
lasso_trIPWRA_rmse=sqrt(mean((lasso_trIPWRA-ATT)^2))
lasso_trIPWRA_bias=mean (mean(lasso_trIPWRA)-ATT)
lasso_trIPWRA_avtime=mean (lasso_trIPWRA_time)

randfor_trIPWRA_mean=mean(randfor_trIPWRA)
randfor_trIPWRA_variance=mean((randfor_trIPWRA-mean(randfor_trIPWRA))^2)
randfor_trIPWRA_rmse=sqrt(mean((randfor_trIPWRA-ATT)^2))
randfor_trIPWRA_bias=mean (mean(randfor_trIPWRA)-ATT)
randfor_trIPWRA_avtime=mean (randfor_trIPWRA_time)

trweightRA_mean=mean(trweightRA)
trweightRA_variance=mean((trweightRA-mean(trweightRA))^2)
trweightRA_rmse=sqrt(mean((trweightRA-ATT)^2))
trweightRA_bias=mean (mean(trweightRA)-ATT)
trweightRA_avtime=mean (trweightRA_time)


#construct a summary table

mean=cbind(did_mean,did2_mean,did3_mean, ipwdid_mean, ordid_mean, impdrdid_mean, lasso_impdrdid_mean,
           randfor_impdrdid_mean, trIPWRA_mean,lasso_trIPWRA_mean, randfor_trIPWRA_mean,
           trweightRA_mean)
bias=cbind(did_bias,did2_bias,did3_bias, ipwdid_bias, ordid_bias, impdrdid_bias, lasso_impdrdid_bias,
           randfor_impdrdid_bias, trIPWRA_bias,lasso_trIPWRA_bias, randfor_trIPWRA_bias,
           trweightRA_bias)
rmse=cbind(did_rmse,did2_rmse,did3_rmse, ipwdid_rmse, ordid_rmse, impdrdid_rmse, lasso_impdrdid_rmse,
           randfor_impdrdid_rmse, trIPWRA_rmse,lasso_trIPWRA_rmse, randfor_trIPWRA_rmse,
           trweightRA_rmse)
variance=cbind(did_variance,did2_variance,did3_variance, ipwdid_variance, ordid_variance, impdrdid_variance, lasso_impdrdid_variance,
               randfor_impdrdid_variance, trIPWRA_variance,lasso_trIPWRA_variance, randfor_trIPWRA_variance,
               trweightRA_variance)
time=cbind(did_avtime,did2_avtime,did3_avtime, ipwdid_avtime, ordid_avtime, impdrdid_avtime, lasso_impdrdid_avtime,
           randfor_impdrdid_avtime, trIPWRA_avtime,lasso_trIPWRA_avtime, randfor_trIPWRA_avtime,
           trweightRA_avtime)

tab <- matrix(1:48, ncol=4, byrow=TRUE)
tab[,1]=bias
tab[,2]=rmse
tab[,3]=variance
tab[,4]=time
tab=round(tab, digits = 3)


colnames(tab) <- c('Bias','RMSE', 'Variance','Time')
rownames(tab) <- c('TWFE','TWFE (T*X)','TWFE (T*X+D*X)','IPW', 'RA','DRDiD', 'LASSO DRDiD', 'RF DRDiD',
                   '3IPWRA', 'LASSO 3IPWRA', 'RF 3IPWRA','3WDRDiD')
latextable=stargazer(tab)
tab
#save table
write.table(tab, file=paste(file.name,'.txt',sep = ""))

#Create plot to visualize results
tabdata=as.data.frame(tab)
tabdata$name=row.names(tabdata)
ggplot(data=tabdata, aes(x=reorder(name, -abs(Bias)), y=abs(Bias))) + 
  geom_bar(stat = "identity")+
  scale_y_continuous(limits=c(0, 33))+
  coord_flip()+
  geom_point(aes(y=RMSE),
             stat="identity",
             position="dodge",
             alpha=1,
             shape=21,
             stroke = 1.8,
             size=3,
             colour = "white",
             fill = "black")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())
ggsave("EXP1B.png")
#commands to read table and transfer to latex
#stargazer(read.table(paste('EXP0B','.txt',sep = "")), summary=FALSE)

#EXPERIMENT 1C: NON-RANDOMIZED EXPERIMENT WITH X-SPECIFIC COMMON TREND, TIME-INVARIANT COVARIATES AND HOMOGENEOUS EFFECTS
#              PROPENSITY SCORE CORRECTLY SPECIFIED, OUTCOME REGRESSION NOT CORRECTLY SPECIFIED


# clean the R environment
rm(list = ls())

repository="C:/Users/tommy/OneDrive/Desktop/tesi/DID simulation/MC simulation"
setwd(repository)
file.name="EXP_1C"

# install.packages("devtools")
library(devtools)
#install_github("clu0/diffindiff")
# devtools::install_github("pedrohcgs/DRDID")
library(diffindiff)

library(readstata13)
library(DRDID)
library(glmnet)
library(pacman)
library(data.table)
library(fixest)
library(stargazer)
library(dplyr)
library(magrittr)
library(sandwich)
library(miceadds)
library(ncpen)
library(GGally)
library(MASS)
library(tidyverse)
library(fabricatr)
library(JWileymisc)
library(rlearner)
library(EQL)
library(ggplot2)
library(readstata13)
library(miceadds)
library(randomForest)
library(party)

#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################

#create triplematching function for later use

tripleIPWRA <- function(id, t, d, X, Y, method='logit'){
  
  # triple DiD propensity score weighting following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X matrix of covariates
  # t is the time dummy
  # d is the treatment group dummy
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  # method selects the estimation method for the propensity score. The available
  #estimation methods are:"logit", "probit", "lasso" and "randomforest".
  
  
  
  
  # create four groups and prepare dataset
  T1=ifelse(t==1& d==1,1,0)
  T0=ifelse(t==0& d==1,1,0)
  C1=ifelse(t==1& d==0,1,0)
  C0=ifelse(t==0& d==0,1,0)
  
  data=as.data.frame(cbind(Y,id, t, d, C0, C1,T0, T1, X))
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=X[T1==1|C1==1,]
  
  #Estimate propensity score
  if (method=='logit'){
    mylogit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                   data = data1)
    data1$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                    data = data1)
    data1$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov1 <- poly(as.matrix(data1[,9:ncol(data1)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov1, data1$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov1, type='response')
    
  } else if (method=='randomforest'){
    data1rf=data1[,8:ncol(data1)] 
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X1))))
    fr.fit=cforest(T1~ .,data=data1rf, controls=mycontrols)
    data1$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data1$pscore[data1$pscore>1]=1
  data1$pscore[data1$pscore<0]=0
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=X[T1==1|C0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                   data = data2)
    data2$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                    data = data2)
    data2$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression
    # perform LASSO linear regression with 10-fold cross validation
    cov2 <- poly(as.matrix(data2[,9:ncol(data2)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov2, data2$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov2, type="response")
    
  } else if (method=='randomforest'){
    data2rf=data2[,8:ncol(data2)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X2))))
    fr.fit=cforest(T1~ .,data=data2rf, controls=mycontrols)
    data2$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data2$pscore[data2$pscore>1]=1
  data2$pscore[data2$pscore<0]=0
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=X[T1==1|T0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                   data = data3)
    data3$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                    data = data3)
    data3$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # Estimate propensity score with LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov3 <- poly(as.matrix(data3[,9:ncol(data3)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov3, data3$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov3, type="response")
    
  } else if (method=='randomforest'){
    
    #Estimate propensity score with random forest
    data3rf=data3[,8:ncol(data3)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X3))))
    fr.fit=cforest(T1~ .,data=data3rf, controls=mycontrols)
    data3$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  
  #Check that the p-value is in between 0 and 0
  data3$pscore[data3$pscore>1]=1
  data3$pscore[data3$pscore<0]=0
  
  
  #Merge propensity score
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  #trimming for pscore value close to 0 and to 1
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  
  #Preparing data for regression
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  Xmatrix=as.matrix(cs_X)
  
  #Computing weights
  cs_data$w_att=rep(1,cs_len)
  cs_data$w_att=ifelse(cs_data$T1==0, cs_data$pscore/(1-cs_data$pscore),1)
  
  
  #Regression with propensity score weights
  e_tripleIPWRA=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix
                                      +I(Xmatrix*t)+I(Xmatrix*d),weights=cs_data$w_att,cluster="d")
  
  return(e_tripleIPWRA)
  
}


###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Repeated Cross Section Data
###################################################################################
# Pre = T stands for pre-treatment period
# treat = F  stands for control group

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = F){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))
  
  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }
  
  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))
  
  # return fitted values
  return(list(out.reg = out.delta))
  
}

# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)
  
  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]
  
  if (!any(is.na(gam))) {
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw)
    
    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean)
    
    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}
#Loss function for estimation of the Bias reduced PS, based on Graham, Pinton and Egel (2012, 2016)

loss.ps.IPT <- function(gamma1, n, D, int.cov, iw){
  #Coefficients for quadratic extrapolation
  cn <- -(n - 1)
  bn <- -n + (n - 1) * log(n - 1)
  an <- -(n - 1) * (1 - log(n - 1) + 0.5 * (log(n - 1))^2)
  vstar <- log(n - 1)
  
  v <- gamma1 %*% t(int.cov)
  phi <- ifelse(v < vstar, - v - exp(v), an + bn * v + 0.5 * cn * (v^2))
  phi1 <- ifelse(v < vstar, - 1 - exp(v), bn + cn * v)
  phi2 <- ifelse(v < vstar, - exp(v), cn)
  
  #phi <- (v<vstar) * (- v - exp(v)) + (v>=vstar) * (an + bn * v + 0.5 * cn * (v^2))
  #phi1 <- (v<vstar) * (- 1 - exp(v)) + (v>=vstar) * (bn  + cn * v)
  #phi2 <- (v<vstar) * (- exp(v)) + (v>=vstar) * cn
  
  # Minus is because nlm minimizes functions, and we aim to maximize!
  res <- - sum(iw * (1 - D) * phi + v)
  
  attr(res, "gradient") <- - t(int.cov) %*% as.vector(iw * ((1-D) * phi1 + 1))
  attr(res, "hessian")  <-  - t(as.vector((1-D) * iw * phi2) * int.cov) %*% int.cov
  return(res)
}

###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  pslogit <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial",
                                         weights = i.weights))
  
  init.gamma <- suppressWarnings(stats::coef(pslogit))
  
  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))
  
  flag <- ifelse(pscore.cal$converged, 0, 1)
  
  gamma.cal <- try(pscore.cal$argument)
  
  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {
    
    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)
    
    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }
  
  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))
  
  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")
    
    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }
    
  }
  
  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))
  
}


lasso_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type =  "weighted", nboot = NULL,
                          inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates), degree=3, raw=TRUE)
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}


randforest_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                               boot = FALSE, boot.type =  "weighted", nboot = NULL,
                               inffunc = FALSE){
  #----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by RANDOM FOREST
  data=as.data.frame(cbind(D, int.cov))
  mycontrols <- cforest_unbiased(ntree=100, as.integer(sqrt(ncol(int.cov))))
  fr.fit=cforest(D~ .,data=data, controls=mycontrols)
  ps.fit=as.numeric(predict(fr.fit, type='prob'))
  
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using RANDOM FOREST.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  reg.cont.coeff.pre=cforest(y_C0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C0*ipw_C0)
  out.y.cont.pre=as.numeric(predict(reg.cont.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  reg.cont.coeff.post <- cforest(y_C1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C1*ipw_C1)
  out.y.cont.post=as.numeric(predict(reg.cont.coeff.post, newdata=data))
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  reg.treat.coeff.pre <- cforest(y_T0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T0)
  
  out.y.treat.pre=as.numeric(predict(reg.treat.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  reg.treat.coeff.post <- cforest(y_T1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T1)
  out.y.treat.post=as.numeric(predict(reg.treat.coeff.post, newdata=data))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}

triple_drdid_rc <- function(y, post, D,id, int.cov, i.weights = NULL,
                            boot = FALSE, boot.type =  "weighted", nboot = NULL,
                            inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore using triplematching tecnique
  T0=ifelse(dta_long$post==0 & dta_long$d==1, 1, 0)
  T1=ifelse(dta_long$post==1 & dta_long$d==1, 1, 0)
  C0=ifelse(dta_long$post==0 & dta_long$d==0, 1, 0)
  C1=ifelse(dta_long$post==1 & dta_long$d==0, 1, 0)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, post, d, covariates, y))
  
  len=nrow(data)
  
  
  #First propensity score matching between T1 and C1
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(covariates, T1==1 | C1==1)
  X1=as.matrix(X1)
  myprobit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                  data = data1)
  
  data1$pscore=myprobit$fitted.values
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(covariates, T1==1 | C0==1)
  X2=as.matrix(X2)
  myprobit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                  data = data2)
  
  data2$pscore=myprobit$fitted.values
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(covariates, T1==1 | T0==1)
  X3=as.matrix(X3)
  myprobit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                  data = data3)
  
  data3$pscore=myprobit$fitted.values
  
  
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  # create dataset with common support and trimming for pscore value close to 0 and to 1
  
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  cs_len=nrow(cs_data)
  
  
  cs_data$w_att=rep(1,cs_len)
  
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  D=cs_data$d
  y=cs_data$y
  cs_int.cov <- as.matrix(rep(1,nrow(cs_data)))
  int.cov=as.matrix(cbind(cs_int.cov,cs_data[,8:11]))
  post=cs_data$post
  pscore.ipt<- cs_data$pscore
  i.weights=i.weights[cs_data$id]
  ps.fit <- as.vector(pscore.ipt)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.post <-  as.vector(out.y.post$out.reg)
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
  
  
}

#create rep.row function for later use

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
# START OF THE SIMULATION
################################################################################

## Parameters and seed

set.seed(1)      # Seed
n = 1000         # Sample size
M = 500         # Number of experiments/iterations




# Results Storage 

did <- rep(0,M)
did2 <- rep(0,M)
did3 <- rep(0,M)
ipwdid<- rep(0,M)
ordid<- rep(0,M)
impdrdid<- rep(0,M)
lasso_impdrdid<- rep(0,M)
randfor_impdrdid<- rep(0,M)
trIPWRA<- rep(0,M)
lasso_trIPWRA<- rep(0,M)
randfor_trIPWRA<- rep(0,M)
trIPWRA2<- rep(0,M)
lasso_trIPWRA2<- rep(0,M)
randfor_trIPWRA2<- rep(0,M)
trweightRA<- rep(0,M)
trweightRA2<- rep(0,M)
did_time <- rep(0,M)
did2_time <- rep(0,M)
did3_time <- rep(0,M)
ipwdid_time<- rep(0,M)
ordid_time<- rep(0,M)
impdrdid_time<- rep(0,M)
lasso_impdrdid_time<- rep(0,M)
randfor_impdrdid_time<- rep(0,M)
trIPWRA_time<- rep(0,M)
lasso_trIPWRA_time<- rep(0,M)
randfor_trIPWRA_time<- rep(0,M)
trIPWRA2_time<- rep(0,M)
lasso_trIPWRA2_time<- rep(0,M)
randfor_trIPWRA2_time<- rep(0,M)
trweightRA_time<- rep(0,M)
trweightRA2_time<- rep(0,M)
ATT<- rep(0,M)

# START OF THE SIMULATION LOOP

for (i in 1:M){ #  M is the number of iterations
  
  
  # Sample size
  n <- 1000
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Proportion in each period
  lambda <- 0.5
  # NUmber of bootstrapped draws
  nboot <- 199
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  
  # Create heterogenenous effects for the ATT, which is set approximately equal to zero
  index.unobs.het <- d * (index.lin)
  index.att <- 0
  
  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  #Gen realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)
  
  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects
  
  # Gen realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  
  # Generate "T"
  post <- as.numeric(stats::runif(n) <= lambda)
  
  # observed outcome
  y <- post * y1 + (1 - post) * y0
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a long data frame
  dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- dta_long[order(dta_long$id),]  
  
  #Standard standard TWFE
  start.time <- Sys.time()
  did_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did[i] <- did_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did2_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did2_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did2[i] <- did2_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did3_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post
               +x1*d+x2*d+x3*d+x4*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did3_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did3[i] <- did3_i$coefficients["post:d"]
  
  # Implement DID ipw Abadie(2005) with normalized weights
  start.time <- Sys.time()
  ipwdid_i=ipwdid(yname="y", tname = "post", idname = "id",  dname = "d",
                  xformla= ~ x1 + x2 + x3 + x4,
                  data = dta_long, panel = FALSE,
                  boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ipwdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ipwdid[i] <- ipwdid_i$ATT
  
  # Implement Outcome regression model following Heckman, Ichimura and Todd (1997) but with linear assumption of covariates
  
  start.time <- Sys.time()
  ordid_i=ordid(yname="y", tname = "post", idname = "id",  dname = "d",
                xformla= ~ x1 + x2 + x3 + x4,
                data = dta_long, panel = FALSE,
                boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ordid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ordid[i] <- ordid_i$ATT
  
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020):
  start.time <- Sys.time()
  impdrdid_i <- drdid(yname="y", tname = "post", idname = "id",  dname = "d",
                      xformla= ~ x1 + x2 + x3 + x4,
                      data = dta_long, panel = FALSE, estMethod = "imp")
  end.time <-Sys.time()
  tk <- end.time - start.time
  impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  impdrdid[i] <- impdrdid_i$ATT
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with LASSO:
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  boot = FALSE
  boot.type =  "weighted"
  nboot = NULL
  inffunc = FALSE
  
  start.time <- Sys.time()
  lasso_impdrdid[i]=lasso_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                   boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with RANDOM FOREST:
  start.time <- Sys.time()
  randfor_impdrdid[i]=randforest_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                          boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  
  # Implement logit tripleIPWRA building on Bludell(2004)
  
  X=(cbind(dta_long$x1,dta_long$x2,dta_long$x3,dta_long$x4))
  id=dta_long$id
  t=dta_long$post
  d=dta_long$d
  Y=dta_long$y
  
  start.time <- Sys.time()
  trIPWRA_i=tripleIPWRA(id, t, d, X, Y)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trIPWRA[i] <- coef(trIPWRA_i)["t:d"]
  
  # Implement LASSO tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  lasso_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='lasso')
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  lasso_trIPWRA[i] <- coef(lasso_trIPWRA_i)["t:d"]
  
  # Implement random forest tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  randfor_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='randomforest')
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  randfor_trIPWRA[i] <- coef(randfor_trIPWRA_i)["t:d"]
  
  # Implement logit tripleDRDiD building by on Sant'Anna and Zhao (2020) 
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  id=dta_long$id
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  
  
  start.time <- Sys.time()
  trweightRA_i=triple_drdid_rc(y=y, post=post, D=D, id, covariates, i.weights = NULL, boot = FALSE,
                               boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trweightRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trweightRA[i] <- trweightRA_i
  
  
  print(i)
}


#generate summary statistics for each estimator

did_mean=mean(did)
did_variance=mean((did-mean(did))^2)
did_rmse=sqrt(mean((did-ATT)^2))
did_bias=mean (mean(did)-ATT)
did_avtime=mean (did_time)

did2_mean=mean(did2)
did2_variance=mean((did2-mean(did2))^2)
did2_rmse=sqrt(mean((did2-ATT)^2))
did2_bias=mean (mean(did2)-ATT)
did2_avtime=mean (did2_time)

did3_mean=mean(did3)
did3_variance=mean((did3-mean(did3))^2)
did3_rmse=sqrt(mean((did3-ATT)^2))
did3_bias=mean (mean(did3)-ATT)
did3_avtime=mean (did3_time)

ipwdid_mean=mean(ipwdid)
ipwdid_variance=mean((ipwdid-mean(ipwdid))^2)
ipwdid_rmse=sqrt(mean((ipwdid-ATT)^2))
ipwdid_bias=mean (mean(ipwdid)-ATT)
ipwdid_avtime=mean (ipwdid_time)

ordid_mean=mean(ordid)
ordid_variance=mean((ordid-mean(ordid))^2)
ordid_rmse=sqrt(mean((ordid-ATT)^2))
ordid_bias=mean (mean(ordid)-ATT)
ordid_avtime=mean (ordid_time)

impdrdid_mean=mean(impdrdid)
impdrdid_variance=mean((impdrdid-mean(impdrdid))^2)
impdrdid_rmse=sqrt(mean((impdrdid-ATT)^2))
impdrdid_bias=mean (mean(impdrdid)-ATT)
impdrdid_avtime=mean (impdrdid_time)

lasso_impdrdid_mean=mean(lasso_impdrdid)
lasso_impdrdid_variance=mean((lasso_impdrdid-mean(lasso_impdrdid))^2)
lasso_impdrdid_rmse=sqrt(mean((lasso_impdrdid-ATT)^2))
lasso_impdrdid_bias=mean (mean(lasso_impdrdid)-ATT)
lasso_impdrdid_avtime=mean (lasso_impdrdid_time)

randfor_impdrdid_mean=mean(randfor_impdrdid)
randfor_impdrdid_variance=mean((randfor_impdrdid-mean(randfor_impdrdid))^2)
randfor_impdrdid_rmse=sqrt(mean((randfor_impdrdid-ATT)^2))
randfor_impdrdid_bias=mean (mean(randfor_impdrdid)-ATT)
randfor_impdrdid_avtime=mean (randfor_impdrdid_time)

trIPWRA_mean=mean(trIPWRA)
trIPWRA_variance=mean((trIPWRA-mean(trIPWRA))^2)
trIPWRA_rmse=sqrt(mean((trIPWRA-ATT)^2))
trIPWRA_bias=mean (mean(trIPWRA)-ATT)
trIPWRA_avtime=mean (trIPWRA_time)

lasso_trIPWRA_mean=mean(lasso_trIPWRA)
lasso_trIPWRA_variance=mean((lasso_trIPWRA-mean(lasso_trIPWRA))^2)
lasso_trIPWRA_rmse=sqrt(mean((lasso_trIPWRA-ATT)^2))
lasso_trIPWRA_bias=mean (mean(lasso_trIPWRA)-ATT)
lasso_trIPWRA_avtime=mean (lasso_trIPWRA_time)

randfor_trIPWRA_mean=mean(randfor_trIPWRA)
randfor_trIPWRA_variance=mean((randfor_trIPWRA-mean(randfor_trIPWRA))^2)
randfor_trIPWRA_rmse=sqrt(mean((randfor_trIPWRA-ATT)^2))
randfor_trIPWRA_bias=mean (mean(randfor_trIPWRA)-ATT)
randfor_trIPWRA_avtime=mean (randfor_trIPWRA_time)

trweightRA_mean=mean(trweightRA)
trweightRA_variance=mean((trweightRA-mean(trweightRA))^2)
trweightRA_rmse=sqrt(mean((trweightRA-ATT)^2))
trweightRA_bias=mean (mean(trweightRA)-ATT)
trweightRA_avtime=mean (trweightRA_time)


#construct a summary table

mean=cbind(did_mean,did2_mean,did3_mean, ipwdid_mean, ordid_mean, impdrdid_mean, lasso_impdrdid_mean,
           randfor_impdrdid_mean, trIPWRA_mean,lasso_trIPWRA_mean, randfor_trIPWRA_mean,
           trweightRA_mean)
bias=cbind(did_bias,did2_bias,did3_bias, ipwdid_bias, ordid_bias, impdrdid_bias, lasso_impdrdid_bias,
           randfor_impdrdid_bias, trIPWRA_bias,lasso_trIPWRA_bias, randfor_trIPWRA_bias,
           trweightRA_bias)
rmse=cbind(did_rmse,did2_rmse,did3_rmse, ipwdid_rmse, ordid_rmse, impdrdid_rmse, lasso_impdrdid_rmse,
           randfor_impdrdid_rmse, trIPWRA_rmse,lasso_trIPWRA_rmse, randfor_trIPWRA_rmse,
           trweightRA_rmse)
variance=cbind(did_variance,did2_variance,did3_variance, ipwdid_variance, ordid_variance, impdrdid_variance, lasso_impdrdid_variance,
               randfor_impdrdid_variance, trIPWRA_variance,lasso_trIPWRA_variance, randfor_trIPWRA_variance,
               trweightRA_variance)
time=cbind(did_avtime,did2_avtime,did3_avtime, ipwdid_avtime, ordid_avtime, impdrdid_avtime, lasso_impdrdid_avtime,
           randfor_impdrdid_avtime, trIPWRA_avtime,lasso_trIPWRA_avtime, randfor_trIPWRA_avtime,
           trweightRA_avtime)

tab <- matrix(1:48, ncol=4, byrow=TRUE)
tab[,1]=bias
tab[,2]=rmse
tab[,3]=variance
tab[,4]=time
tab=round(tab, digits = 3)


colnames(tab) <- c('Bias','RMSE', 'Variance','Time')
rownames(tab) <- c('TWFE','TWFE (T*X)','TWFE (T*X+D*X)','IPW', 'RA','DRDiD', 'LASSO DRDiD', 'RF DRDiD',
                   '3IPWRA', 'LASSO 3IPWRA', 'RF 3IPWRA','3WDRDiD')
latextable=stargazer(tab)
tab

#save table
write.table(tab, file=paste(file.name,'.txt',sep = ""))

#Create plot to visualize results
tabdata=as.data.frame(tab)
tabdata$name=row.names(tabdata)
ggplot(data=tabdata, aes(x=reorder(name, -abs(Bias)), y=abs(Bias))) + 
  geom_bar(stat = "identity")+
  scale_y_continuous(limits=c(0, 33))+
  coord_flip()+
  geom_point(aes(y=RMSE),
             stat="identity",
             position="dodge",
             alpha=1,
             shape=21,
             stroke = 1.8,
             size=3,
             colour = "white",
             fill = "black")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())
ggsave("EXP1C.png")

#commands to read table and transfer to latex
#stargazer(read.table(paste('EXP0B','.txt',sep = "")), summary=FALSE)

#EXPERIMENT 1D: NON-RANDOMIZED EXPERIMENT WITH X-SPECIFIC COMMON TREND, TIME-INVARIANT COVARIATES AND HOMOGENEOUS EFFECTS
#              PROPENSITY SCORE NOT CORRECTLY SPECIFIED, OUTCOME REGRESSION NOT CORRECTLY SPECIFIED


# clean the R environment
rm(list = ls())

repository="C:/Users/tommy/OneDrive/Desktop/tesi/DID simulation/MC simulation"
setwd(repository)
file.name="EXP_1D"

# install.packages("devtools")
library(devtools)
#install_github("clu0/diffindiff")
# devtools::install_github("pedrohcgs/DRDID")
library(diffindiff)

library(readstata13)
library(DRDID)
library(glmnet)
library(pacman)
library(data.table)
library(fixest)
library(stargazer)
library(dplyr)
library(magrittr)
library(sandwich)
library(miceadds)
library(ncpen)
library(GGally)
library(MASS)
library(tidyverse)
library(fabricatr)
library(JWileymisc)
library(rlearner)
library(EQL)
library(ggplot2)
library(readstata13)
library(miceadds)
library(randomForest)
library(party)

#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################

#create triplematching function for later use

tripleIPWRA <- function(id, t, d, X, Y, method='logit'){
  
  # triple DiD propensity score weighting following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X matrix of covariates
  # t is the time dummy
  # d is the treatment group dummy
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  # method selects the estimation method for the propensity score. The available
  #estimation methods are:"logit", "probit", "lasso" and "randomforest".
  
  
  
  
  # create four groups and prepare dataset
  T1=ifelse(t==1& d==1,1,0)
  T0=ifelse(t==0& d==1,1,0)
  C1=ifelse(t==1& d==0,1,0)
  C0=ifelse(t==0& d==0,1,0)
  
  data=as.data.frame(cbind(Y,id, t, d, C0, C1,T0, T1, X))
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=X[T1==1|C1==1,]
  
  #Estimate propensity score
  if (method=='logit'){
    mylogit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                   data = data1)
    data1$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                    data = data1)
    data1$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov1 <- poly(as.matrix(data1[,9:ncol(data1)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov1, data1$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov1, type='response')
    
  } else if (method=='randomforest'){
    data1rf=data1[,8:ncol(data1)] 
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X1))))
    fr.fit=cforest(T1~ .,data=data1rf, controls=mycontrols)
    data1$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data1$pscore[data1$pscore>1]=1
  data1$pscore[data1$pscore<0]=0
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=X[T1==1|C0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                   data = data2)
    data2$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                    data = data2)
    data2$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression
    # perform LASSO linear regression with 10-fold cross validation
    cov2 <- poly(as.matrix(data2[,9:ncol(data2)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov2, data2$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov2, type="response")
    
  } else if (method=='randomforest'){
    data2rf=data2[,8:ncol(data2)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X2))))
    fr.fit=cforest(T1~ .,data=data2rf, controls=mycontrols)
    data2$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data2$pscore[data2$pscore>1]=1
  data2$pscore[data2$pscore<0]=0
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=X[T1==1|T0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                   data = data3)
    data3$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                    data = data3)
    data3$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # Estimate propensity score with LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov3 <- poly(as.matrix(data3[,9:ncol(data3)]), degree=3, raw=TRUE)
    
    lasso.fit <- cv.glmnet(cov3, data3$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov3, type="response")
    
  } else if (method=='randomforest'){
    
    #Estimate propensity score with random forest
    data3rf=data3[,8:ncol(data3)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X3))))
    fr.fit=cforest(T1~ .,data=data3rf, controls=mycontrols)
    data3$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  
  #Check that the p-value is in between 0 and 0
  data3$pscore[data3$pscore>1]=1
  data3$pscore[data3$pscore<0]=0
  
  
  #Merge propensity score
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  #trimming for pscore value close to 0 and to 1
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  
  #Preparing data for regression
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  Xmatrix=as.matrix(cs_X)
  
  #Computing weights
  cs_data$w_att=rep(1,cs_len)
  cs_data$w_att=ifelse(cs_data$T1==0, cs_data$pscore/(1-cs_data$pscore),1)
  
  
  #Regression with propensity score weights
  e_tripleIPWRA=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix
                                      +I(Xmatrix*t)+I(Xmatrix*d),weights=cs_data$w_att,cluster="d")
  
  return(e_tripleIPWRA)
  
}


###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Repeated Cross Section Data
###################################################################################
# Pre = T stands for pre-treatment period
# treat = F  stands for control group

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = F){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))
  
  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }
  
  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))
  
  # return fitted values
  return(list(out.reg = out.delta))
  
}

# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)
  
  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]
  
  if (!any(is.na(gam))) {
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw)
    
    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean)
    
    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}
#Loss function for estimation of the Bias reduced PS, based on Graham, Pinton and Egel (2012, 2016)

loss.ps.IPT <- function(gamma1, n, D, int.cov, iw){
  #Coefficients for quadratic extrapolation
  cn <- -(n - 1)
  bn <- -n + (n - 1) * log(n - 1)
  an <- -(n - 1) * (1 - log(n - 1) + 0.5 * (log(n - 1))^2)
  vstar <- log(n - 1)
  
  v <- gamma1 %*% t(int.cov)
  phi <- ifelse(v < vstar, - v - exp(v), an + bn * v + 0.5 * cn * (v^2))
  phi1 <- ifelse(v < vstar, - 1 - exp(v), bn + cn * v)
  phi2 <- ifelse(v < vstar, - exp(v), cn)
  
  #phi <- (v<vstar) * (- v - exp(v)) + (v>=vstar) * (an + bn * v + 0.5 * cn * (v^2))
  #phi1 <- (v<vstar) * (- 1 - exp(v)) + (v>=vstar) * (bn  + cn * v)
  #phi2 <- (v<vstar) * (- exp(v)) + (v>=vstar) * cn
  
  # Minus is because nlm minimizes functions, and we aim to maximize!
  res <- - sum(iw * (1 - D) * phi + v)
  
  attr(res, "gradient") <- - t(int.cov) %*% as.vector(iw * ((1-D) * phi1 + 1))
  attr(res, "hessian")  <-  - t(as.vector((1-D) * iw * phi2) * int.cov) %*% int.cov
  return(res)
}

###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  pslogit <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial",
                                         weights = i.weights))
  
  init.gamma <- suppressWarnings(stats::coef(pslogit))
  
  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))
  
  flag <- ifelse(pscore.cal$converged, 0, 1)
  
  gamma.cal <- try(pscore.cal$argument)
  
  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {
    
    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)
    
    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }
  
  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))
  
  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")
    
    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }
    
  }
  
  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))
  
}


lasso_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type =  "weighted", nboot = NULL,
                          inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates), degree=3, raw=TRUE)
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}


randforest_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                               boot = FALSE, boot.type =  "weighted", nboot = NULL,
                               inffunc = FALSE){
  #----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by RANDOM FOREST
  data=as.data.frame(cbind(D, int.cov))
  mycontrols <- cforest_unbiased(ntree=100, as.integer(sqrt(ncol(int.cov))))
  fr.fit=cforest(D~ .,data=data, controls=mycontrols)
  ps.fit=as.numeric(predict(fr.fit, type='prob'))
  
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using RANDOM FOREST.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  reg.cont.coeff.pre=cforest(y_C0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C0*ipw_C0)
  out.y.cont.pre=as.numeric(predict(reg.cont.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  reg.cont.coeff.post <- cforest(y_C1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_C1*ipw_C1)
  out.y.cont.post=as.numeric(predict(reg.cont.coeff.post, newdata=data))
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  reg.treat.coeff.pre <- cforest(y_T0 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T0)
  
  out.y.treat.pre=as.numeric(predict(reg.treat.coeff.pre, newdata=data))
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  reg.treat.coeff.post <- cforest(y_T1 ~ .,data=dataor, controls=mycontrols, weights =i.weights_T1)
  out.y.treat.post=as.numeric(predict(reg.treat.coeff.post, newdata=data))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att)
}

triple_drdid_rc <- function(y, post, D,id, int.cov, i.weights = NULL,
                            boot = FALSE, boot.type =  "weighted", nboot = NULL,
                            inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore using triplematching tecnique
  T0=ifelse(dta_long$post==0 & dta_long$d==1, 1, 0)
  T1=ifelse(dta_long$post==1 & dta_long$d==1, 1, 0)
  C0=ifelse(dta_long$post==0 & dta_long$d==0, 1, 0)
  C1=ifelse(dta_long$post==1 & dta_long$d==0, 1, 0)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, post, d, covariates, y))
  
  len=nrow(data)
  
  
  #First propensity score matching between T1 and C1
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(covariates, T1==1 | C1==1)
  X1=as.matrix(X1)
  myprobit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                  data = data1)
  
  data1$pscore=myprobit$fitted.values
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(covariates, T1==1 | C0==1)
  X2=as.matrix(X2)
  myprobit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                  data = data2)
  
  data2$pscore=myprobit$fitted.values
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(covariates, T1==1 | T0==1)
  X3=as.matrix(X3)
  myprobit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                  data = data3)
  
  data3$pscore=myprobit$fitted.values
  
  
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  # create dataset with common support and trimming for pscore value close to 0 and to 1
  
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  cs_len=nrow(cs_data)
  
  
  cs_data$w_att=rep(1,cs_len)
  
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  D=cs_data$d
  y=cs_data$y
  cs_int.cov <- as.matrix(rep(1,nrow(cs_data)))
  int.cov=as.matrix(cbind(cs_int.cov,cs_data[,8:11]))
  post=cs_data$post
  pscore.ipt<- cs_data$pscore
  i.weights=i.weights[cs_data$id]
  ps.fit <- as.vector(pscore.ipt)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.post <-  as.vector(out.y.post$out.reg)
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
  
  
}

#create rep.row function for later use

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
# START OF THE SIMULATION
################################################################################

## Parameters and seed

set.seed(1)      # Seed
n = 1000         # Sample size
M = 500         # Number of experiments/iterations




# Results Storage 

did <- rep(0,M)
did2 <- rep(0,M)
did3 <- rep(0,M)
ipwdid<- rep(0,M)
ordid<- rep(0,M)
impdrdid<- rep(0,M)
lasso_impdrdid<- rep(0,M)
randfor_impdrdid<- rep(0,M)
trIPWRA<- rep(0,M)
lasso_trIPWRA<- rep(0,M)
randfor_trIPWRA<- rep(0,M)
trIPWRA2<- rep(0,M)
lasso_trIPWRA2<- rep(0,M)
randfor_trIPWRA2<- rep(0,M)
trweightRA<- rep(0,M)
trweightRA2<- rep(0,M)
did_time <- rep(0,M)
did2_time <- rep(0,M)
did3_time <- rep(0,M)
ipwdid_time<- rep(0,M)
ordid_time<- rep(0,M)
impdrdid_time<- rep(0,M)
lasso_impdrdid_time<- rep(0,M)
randfor_impdrdid_time<- rep(0,M)
trIPWRA_time<- rep(0,M)
lasso_trIPWRA_time<- rep(0,M)
randfor_trIPWRA_time<- rep(0,M)
trIPWRA2_time<- rep(0,M)
lasso_trIPWRA2_time<- rep(0,M)
randfor_trIPWRA2_time<- rep(0,M)
trweightRA_time<- rep(0,M)
trweightRA2_time<- rep(0,M)
ATT<- rep(0,M)

# START OF THE SIMULATION LOOP

for (i in 1:M){ #  M is the number of iterations
  
  
  # Sample size
  n <- 1000
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Proportion in each period
  lambda <- 0.5
  # NUmber of bootstrapped draws
  nboot <- 199
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  
  # Create heterogenenous effects for the ATT, which is set approximately equal to zero
  index.unobs.het <- d * (index.lin)
  index.att <- 0
  
  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  #Gen realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)
  
  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects
  
  # Gen realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  
  # Generate "T"
  post <- as.numeric(stats::runif(n) <= lambda)
  
  # observed outcome
  y <- post * y1 + (1 - post) * y0
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a long data frame
  dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- dta_long[order(dta_long$id),]  
  
  #Standard standard TWFE
  start.time <- Sys.time()
  did_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did[i] <- did_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did2_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did2_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did2[i] <- did2_i$coefficients["post:d"]
  
  #TWFE with time interactions
  start.time <- Sys.time()
  did3_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d
               +x1*post+x2*post+x3*post+x4*post
               +x1*d+x2*d+x3*d+x4*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  did3_time[i]=round(as.numeric(substr(tk,0,6)),5)
  did3[i] <- did3_i$coefficients["post:d"]
  
  # Implement DID ipw Abadie(2005) with normalized weights
  start.time <- Sys.time()
  ipwdid_i=ipwdid(yname="y", tname = "post", idname = "id",  dname = "d",
                  xformla= ~ x1 + x2 + x3 + x4,
                  data = dta_long, panel = FALSE,
                  boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ipwdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ipwdid[i] <- ipwdid_i$ATT
  
  # Implement Outcome regression model following Heckman, Ichimura and Todd (1997) but with linear assumption of covariates
  
  start.time <- Sys.time()
  ordid_i=ordid(yname="y", tname = "post", idname = "id",  dname = "d",
                xformla= ~ x1 + x2 + x3 + x4,
                data = dta_long, panel = FALSE,
                boot = FALSE, nboot = 199)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ordid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ordid[i] <- ordid_i$ATT
  
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020):
  start.time <- Sys.time()
  impdrdid_i <- drdid(yname="y", tname = "post", idname = "id",  dname = "d",
                      xformla= ~ x1 + x2 + x3 + x4,
                      data = dta_long, panel = FALSE, estMethod = "imp")
  end.time <-Sys.time()
  tk <- end.time - start.time
  impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  impdrdid[i] <- impdrdid_i$ATT
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with LASSO:
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  boot = FALSE
  boot.type =  "weighted"
  nboot = NULL
  inffunc = FALSE
  
  start.time <- Sys.time()
  lasso_impdrdid[i]=lasso_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                   boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020) with RANDOM FOREST:
  start.time <- Sys.time()
  randfor_impdrdid[i]=randforest_drdid_rc(y=y, post=post, D=D, covariates, i.weights = NULL, boot = FALSE,
                                          boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_impdrdid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  
  
  # Implement logit tripleIPWRA building on Bludell(2004)
  
  X=(cbind(dta_long$x1,dta_long$x2,dta_long$x3,dta_long$x4))
  id=dta_long$id
  t=dta_long$post
  d=dta_long$d
  Y=dta_long$y
  
  start.time <- Sys.time()
  trIPWRA_i=tripleIPWRA(id, t, d, X, Y)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trIPWRA[i] <- coef(trIPWRA_i)["t:d"]
  
  # Implement LASSO tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  lasso_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='lasso')
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  lasso_trIPWRA[i] <- coef(lasso_trIPWRA_i)["t:d"]
  
  # Implement random forest tripleIPWRA building on Bludell(2004)
  
  start.time <- Sys.time()
  randfor_trIPWRA_i=tripleIPWRA(id, t, d, X, Y, method='randomforest')
  end.time <-Sys.time()
  tk <- end.time - start.time
  randfor_trIPWRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  randfor_trIPWRA[i] <- coef(randfor_trIPWRA_i)["t:d"]
  
  # Implement logit tripleDRDiD building by on Sant'Anna and Zhao (2020) 
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  id=dta_long$id
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  
  
  start.time <- Sys.time()
  trweightRA_i=triple_drdid_rc(y=y, post=post, D=D, id, covariates, i.weights = NULL, boot = FALSE,
                               boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  trweightRA_time[i]=round(as.numeric(substr(tk,0,6)),5)
  trweightRA[i] <- trweightRA_i
  
  
  print(i)
}


#generate summary statistics for each estimator

did_mean=mean(did)
did_variance=mean((did-mean(did))^2)
did_rmse=sqrt(mean((did-ATT)^2))
did_bias=mean (mean(did)-ATT)
did_avtime=mean (did_time)

did2_mean=mean(did2)
did2_variance=mean((did2-mean(did2))^2)
did2_rmse=sqrt(mean((did2-ATT)^2))
did2_bias=mean (mean(did2)-ATT)
did2_avtime=mean (did2_time)

did3_mean=mean(did3)
did3_variance=mean((did3-mean(did3))^2)
did3_rmse=sqrt(mean((did3-ATT)^2))
did3_bias=mean (mean(did3)-ATT)
did3_avtime=mean (did3_time)

ipwdid_mean=mean(ipwdid)
ipwdid_variance=mean((ipwdid-mean(ipwdid))^2)
ipwdid_rmse=sqrt(mean((ipwdid-ATT)^2))
ipwdid_bias=mean (mean(ipwdid)-ATT)
ipwdid_avtime=mean (ipwdid_time)

ordid_mean=mean(ordid)
ordid_variance=mean((ordid-mean(ordid))^2)
ordid_rmse=sqrt(mean((ordid-ATT)^2))
ordid_bias=mean (mean(ordid)-ATT)
ordid_avtime=mean (ordid_time)

impdrdid_mean=mean(impdrdid)
impdrdid_variance=mean((impdrdid-mean(impdrdid))^2)
impdrdid_rmse=sqrt(mean((impdrdid-ATT)^2))
impdrdid_bias=mean (mean(impdrdid)-ATT)
impdrdid_avtime=mean (impdrdid_time)

lasso_impdrdid_mean=mean(lasso_impdrdid)
lasso_impdrdid_variance=mean((lasso_impdrdid-mean(lasso_impdrdid))^2)
lasso_impdrdid_rmse=sqrt(mean((lasso_impdrdid-ATT)^2))
lasso_impdrdid_bias=mean (mean(lasso_impdrdid)-ATT)
lasso_impdrdid_avtime=mean (lasso_impdrdid_time)

randfor_impdrdid_mean=mean(randfor_impdrdid)
randfor_impdrdid_variance=mean((randfor_impdrdid-mean(randfor_impdrdid))^2)
randfor_impdrdid_rmse=sqrt(mean((randfor_impdrdid-ATT)^2))
randfor_impdrdid_bias=mean (mean(randfor_impdrdid)-ATT)
randfor_impdrdid_avtime=mean (randfor_impdrdid_time)

trIPWRA_mean=mean(trIPWRA)
trIPWRA_variance=mean((trIPWRA-mean(trIPWRA))^2)
trIPWRA_rmse=sqrt(mean((trIPWRA-ATT)^2))
trIPWRA_bias=mean (mean(trIPWRA)-ATT)
trIPWRA_avtime=mean (trIPWRA_time)

lasso_trIPWRA_mean=mean(lasso_trIPWRA)
lasso_trIPWRA_variance=mean((lasso_trIPWRA-mean(lasso_trIPWRA))^2)
lasso_trIPWRA_rmse=sqrt(mean((lasso_trIPWRA-ATT)^2))
lasso_trIPWRA_bias=mean (mean(lasso_trIPWRA)-ATT)
lasso_trIPWRA_avtime=mean (lasso_trIPWRA_time)

randfor_trIPWRA_mean=mean(randfor_trIPWRA)
randfor_trIPWRA_variance=mean((randfor_trIPWRA-mean(randfor_trIPWRA))^2)
randfor_trIPWRA_rmse=sqrt(mean((randfor_trIPWRA-ATT)^2))
randfor_trIPWRA_bias=mean (mean(randfor_trIPWRA)-ATT)
randfor_trIPWRA_avtime=mean (randfor_trIPWRA_time)

trweightRA_mean=mean(trweightRA)
trweightRA_variance=mean((trweightRA-mean(trweightRA))^2)
trweightRA_rmse=sqrt(mean((trweightRA-ATT)^2))
trweightRA_bias=mean (mean(trweightRA)-ATT)
trweightRA_avtime=mean (trweightRA_time)


#construct a summary table

mean=cbind(did_mean,did2_mean,did3_mean, ipwdid_mean, ordid_mean, impdrdid_mean, lasso_impdrdid_mean,
           randfor_impdrdid_mean, trIPWRA_mean,lasso_trIPWRA_mean, randfor_trIPWRA_mean,
           trweightRA_mean)
bias=cbind(did_bias,did2_bias,did3_bias, ipwdid_bias, ordid_bias, impdrdid_bias, lasso_impdrdid_bias,
           randfor_impdrdid_bias, trIPWRA_bias,lasso_trIPWRA_bias, randfor_trIPWRA_bias,
           trweightRA_bias)
rmse=cbind(did_rmse,did2_rmse,did3_rmse, ipwdid_rmse, ordid_rmse, impdrdid_rmse, lasso_impdrdid_rmse,
           randfor_impdrdid_rmse, trIPWRA_rmse,lasso_trIPWRA_rmse, randfor_trIPWRA_rmse,
           trweightRA_rmse)
variance=cbind(did_variance,did2_variance,did3_variance, ipwdid_variance, ordid_variance, impdrdid_variance, lasso_impdrdid_variance,
               randfor_impdrdid_variance, trIPWRA_variance,lasso_trIPWRA_variance, randfor_trIPWRA_variance,
               trweightRA_variance)
time=cbind(did_avtime,did2_avtime,did3_avtime, ipwdid_avtime, ordid_avtime, impdrdid_avtime, lasso_impdrdid_avtime,
           randfor_impdrdid_avtime, trIPWRA_avtime,lasso_trIPWRA_avtime, randfor_trIPWRA_avtime,
           trweightRA_avtime)

tab <- matrix(1:48, ncol=4, byrow=TRUE)
tab[,1]=bias
tab[,2]=rmse
tab[,3]=variance
tab[,4]=time
tab=round(tab, digits = 3)


colnames(tab) <- c('Bias','RMSE', 'Variance','Time')
rownames(tab) <- c('TWFE','TWFE (T*X)','TWFE (T*X+D*X)','IPW', 'RA','DRDiD', 'LASSO DRDiD', 'RF DRDiD',
                   '3IPWRA', 'LASSO 3IPWRA', 'RF 3IPWRA','3WDRDiD')
latextable=stargazer(tab)
tab
#save table
write.table(tab, file=paste(file.name,'.txt',sep = ""))

#Create plot to visualize results
tabdata=as.data.frame(tab)
tabdata$name=row.names(tabdata)
ggplot(data=tabdata, aes(x=reorder(name, -abs(Bias)), y=abs(Bias))) + 
  geom_bar(stat = "identity")+
  scale_y_continuous(limits=c(0, 33))+
  coord_flip()+
  geom_point(aes(y=RMSE),
             stat="identity",
             position="dodge",
             alpha=1,
             shape=21,
             stroke = 1.8,
             size=3,
             colour = "white",
             fill = "black")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())
ggsave("EXP1D.png")

#commands to read table and transfer to latex
#stargazer(read.table(paste('EXP0B','.txt',sep = "")), summary=FALSE)