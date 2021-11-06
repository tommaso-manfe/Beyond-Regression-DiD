#EXPERIMENT 7: NON-RANDOMIZED EXPERIMENT, NO STATIONARITY OF COVARIATES AND OF THE UNOBSERVED CHARACTERISTIC
# NON-LINEAR DATA GENERATING PROCESS, HETEROGENOUS EFFECTS

# clean the R environment
rm(list = ls())


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

#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################

# create Chang (2020) DMLDiD formula for later use

DMLDiD=function(Y,D,p,T){
  N=length(Y)
  B=100
  set.seed(123)
  random=sample(1:1000,B)
  
  thetabar=c(0)
  for (l in 1:B){
    k=2
    samplesplit=function(k,N){
      c1=1:N
      smp_size <- floor((1/k) * length(c1))
      
      ## set the seed to make your partition reproducible
      set.seed(random[l])
      train_ind <- sample(seq_len(length(c1)), size = smp_size)
      
      k1 <- c1[train_ind]
      k2 <- c1[-train_ind]
      return(rbind(k1,k2))
    }
    K=samplesplit(k,N)
    
    thetaDML=c(0)
    
    for (q in 1:k){
      ##Trimming
      set.seed(333)
      CV=cv.glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1)
      fit=glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1,lambda=CV$lambda.1se)
      beta1hat=fit$beta
      beta1hat <- as.numeric(as.character(beta1hat))
      
      ghat=1/(1+exp(-p[K[q,],]%*%beta1hat))
      
      index1=K[q,][which(ghat<0.95 & ghat>0.05)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=p[-K[q,],]
      XX=XX[index,]
      
      
      
      set.seed(333)
      CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
      fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
      beta2hat=fit$beta
      beta2hat <- as.numeric(as.character(beta2hat))
      
      ellhat2=p[index1,]%*%beta2hat
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1])
      s=s[which(s<abs(min(s)))]
      
      thetaDML[q]=mean(s)
      
    }
    
    thetabar[l]=mean(thetaDML)
    
    
  }
  finaltheta=mean(thetabar)
  finaltheta
  
  ##Variance
  var=c(0)
  for (m in 1:B){
    k=2
    samplesplit=function(k,N){
      c1=1:N
      smp_size <- floor((1/k) * length(c1))
      
      ## set the seed to make your partition reproducible
      set.seed(random[m])
      train_ind <- sample(seq_len(length(c1)), size = smp_size)
      
      k1 <- c1[train_ind]
      k2 <- c1[-train_ind]
      return(rbind(k1,k2))
    }
    K=samplesplit(k,N)
    
    varDML=c(0)
    for (q in 1:k){
      ##Trimming
      set.seed(333)
      CV=cv.glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1)
      fit=glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1,lambda=CV$lambda.1se)
      beta1hat=fit$beta
      beta1hat <- as.numeric(as.character(beta1hat))
      
      ghat=1/(1+exp(-p[K[q,],]%*%beta1hat))
      
      index1=K[q,][which(ghat<0.97 & ghat>0.03)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=p[-K[q,],]
      XX=XX[index,]
      
      set.seed(333)
      CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
      fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
      beta2hat=fit$beta
      beta2hat <- as.numeric(as.character(beta2hat))
      
      ellhat2=p[index1,]%*%beta2hat
      
      
      
      G=-(1-2*lambda)*finaltheta/(lambda*(1-lambda))-mean(Y[index1]*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1]))
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/mean(D[index1])/(lambda*(1-lambda))-D[index1]*finaltheta/mean(D[index1])+G*(T[index1]-lambda)
      s=s[which(s<abs(min(s)))]
      
      varDML[q]=mean(s^2)
    }
    
    var[m]=mean(varDML)
  }
  
  sd=sqrt(mean(var))/sqrt(N)
  sd
  return(c(finaltheta,sd))
}


#create triplematching function for later use

triplematching <- function(id, T0, T1, C0, C1, t, d, X, Y){
  
  #triple DiD propensity score matching following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X must be specified as a dataframe!
  # t is the time effect
  # d is the treatment group effect
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  
  library(ggplot2)
  library(readstata13)
  library(miceadds)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, t, d, X, Y))
  
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(X, T1==1 | C1==1)
  X1=as.matrix(X1)
  myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                  data = data1)
  
  data1$pscore=myprobit$fitted.values
  
  # propensity scores graphs treated and untreated
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data1, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT1=subset(data1, T1==1)
  dataC1=subset(data1, C1==1)
  
  minT1=min(dataT1$pscore)   #treated pscore range
  maxT1=max(dataT1$pscore)
  
  minC1=min(dataC1$pscore)   #untreated pscore range
  maxC1=max(dataC1$pscore)
  
  data1$support= rep(0,nrow(data1))
  
  for (k in 1:K){
    if (data1$pscore[k]>=minT1 & data1$pscore[k]<=maxC1){
      data1$support[k]=1
    }
  }
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(X, T1==1 | C0==1)
  X2=as.matrix(X2)
  myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                  data = data2)
  
  data2$pscore=myprobit$fitted.values
  
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data2, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT2=subset(data2, T1==1)
  dataC2=subset(data2, C0==1)
  
  minT2=min(dataT2$pscore)   #treated pscore range
  maxT2=max(dataT2$pscore)
  
  minC2=min(dataC2$pscore)   #untreated pscore range
  maxC2=max(dataC2$pscore)
  
  data2$support= rep(0,nrow(data2))
  
  for (k in 1:K){
    if (data2$pscore[k]>=minT2 & data2$pscore[k]<=maxC2){
      data2$support[k]=1
    }
  }
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(X, T1==1 | T0==1)
  X3=as.matrix(X3)
  myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                  data = data3)
  
  data3$pscore=myprobit$fitted.values
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data3, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT3=subset(data3, T1==1)
  dataC3=subset(data3, T0==1)
  
  minT3=min(dataT3$pscore)   #treated pscore range
  maxT3=max(dataT3$pscore)
  
  minC3=min(dataC3$pscore)   #untreated pscore range
  maxC3=max(dataC3$pscore)
  
  data3$support= rep(0,nrow(data3))
  
  for (k in 1:K){
    if (data3$pscore[k]>=minT3 & data3$pscore[k]<=maxC3){
      data3$support[k]=1
    }
  }
  
  data$pscore=rep(0,len)
  data$commonsupport=rep(0,len)
  data$commonsupport[data1$id]=data1$support[1:length(data1$id)]
  data$pscore[data1$id]=data1$pscore[1:length(data1$id)]
  data$commonsupport[data2$id]=data2$support[1:length(data2$id)]
  data$pscore[data2$id]=data2$pscore[1:length(data2$id)]
  data$commonsupport[data3$id]=data3$support[1:length(data3$id)]
  data$pscore[data3$id]=data3$pscore[1:length(data3$id)]
  
  # create dataset with common support and trimming for pscore value close to 0 and to 1
  
  cs_data=subset(data, commonsupport==1)
  cs_data=subset(cs_data, pscore>=0.02 & pscore<=0.98)
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  
  
  cs_data$w_att=rep(1,cs_len)
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  
  
  Xmatrix=as.matrix(cs_X)
  triplematching=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix,weights=cs_data$w_att,cluster="d")
  
  # version with non-clustered standard errors would be 
  # triplematching= lm(Y ~ t*d+t+d+Xmatrix, data = cs_data, weights=w_att) 
  
  return(triplematching)
}

triplematchinglasso <- function(id, T0, T1, C0, C1, t, d, X, Y){
  
  #triple DiD propensity score matching following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X must be specified as a dataframe!
  # t is the time effect
  # d is the treatment group effect
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  #
  
  library(ggplot2)
  library(readstata13)
  library(miceadds)
  
  # would work also using matricex X=as.matrix(cbind(X))
  
  
  data=as.data.frame(cbind(id,T0, T1, C0, C1, t, d, X, Y))
  
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  
  
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=subset(X, T1==1 | C1==1)
  X1=as.matrix(X1)
  
  # LASSO linear regression 
  
  # Split data into training and testing datasets.
  # 0.7 of the data will be used for Training and 0.3 of the
  # data will be used for Testing.
  
  train_rows <- sample(1:K, .7*K)
  x.train <- X1[train_rows, ]
  x.test <- X1[-train_rows, ]
  
  y.train <- data1$T1[train_rows]
  y.test <- data1$T1[-train_rows]
  
  # perform LASSO linear regression with 10-fold cross validation
  
  lasso.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                         alpha=1, family="gaussian")
  
  data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=X1)
  
  for (k in 1:K){
    if (data1$pscore[k]<0) {
      data1$pscore[k]=0
    } else if (data1$pscore[k]>1) {
      data1$pscore[k]=1
    }
  }
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data1, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT1=subset(data1, T1==1)
  dataC1=subset(data1, C1==1)
  
  minT1=min(dataT1$pscore)   #treated pscore range
  maxT1=max(dataT1$pscore)
  
  minC1=min(dataC1$pscore)   #untreated pscore range
  maxC1=max(dataC1$pscore)
  
  data1$support= rep(0,nrow(data1))
  
  for (k in 1:K){
    if (data1$pscore[k]>=minT1 & data1$pscore[k]<=maxC1){
      data1$support[k]=1
    }
  }
  
  #Second propensity score matching between T1 and C0
  
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=subset(X, T1==1 | C0==1)
  X2=as.matrix(X2)
  
  # LASSO linear regression 
  
  # Split data into training and testing datasets.
  # 0.7 of the data will be used for Training and 0.3 of the
  # data will be used for Testing.
  
  train_rows <- sample(1:K, .7*K)
  x.train <- X2[train_rows, ]
  x.test <- X2[-train_rows, ]
  
  y.train <- data2$T1[train_rows]
  y.test <- data2$T1[-train_rows]
  
  # perform LASSO linear regression with 10-fold cross validation
  
  lasso.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                         alpha=1, family="gaussian")
  
  data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=X2)
  
  for (k in 1:K){
    if (data2$pscore[k]<0) {
      data2$pscore[k]=0
    } else if (data2$pscore[k]>1) {
      data2$pscore[k]=1
    }
  }
  
  
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data2, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT2=subset(data2, T1==1)
  dataC2=subset(data2, C0==1)
  
  minT2=min(dataT2$pscore)   #treated pscore range
  maxT2=max(dataT2$pscore)
  
  minC2=min(dataC2$pscore)   #untreated pscore range
  maxC2=max(dataC2$pscore)
  
  data2$support= rep(0,nrow(data2))
  
  for (k in 1:K){
    if (data2$pscore[k]>=minT2 & data2$pscore[k]<=maxC2){
      data2$support[k]=1
    }
  }
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=subset(X, T1==1 | T0==1)
  X3=as.matrix(X3)
  
  # LASSO linear regression 
  
  # Split data into training and testing datasets.
  # 0.7 of the data will be used for Training and 0.3 of the
  # data will be used for Testing.
  
  train_rows <- sample(1:K, .7*K)
  x.train <- X3[train_rows, ]
  x.test <- X3[-train_rows, ]
  
  y.train <- data3$T1[train_rows]
  y.test <- data3$T1[-train_rows]
  
  # perform LASSO linear regression with 10-fold cross validation
  
  lasso.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                         alpha=1, family="gaussian")
  
  data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=X3)
  
  for (k in 1:K){
    if (data3$pscore[k]<0) {
      data3$pscore[k]=0
    } else if (data3$pscore[k]>1) {
      data3$pscore[k]=1
    }
  }
  
  
  #labs <- paste("Propensity score", c("treated", "untreated"))
  #prs_df %>%
  #  mutate(T1 = ifelse(T1 == 1, labs[1], labs[2])) %>
  #  ggplot(data3, aes(x=pscore)) +
  #  geom_histogram(color = "white") +
  #  facet_wrap(~T1) +
  #  xlab("Propensity score") +
  #  theme_bw()
  
  #common support
  
  dataT3=subset(data3, T1==1)
  dataC3=subset(data3, T0==1)
  
  minT3=min(dataT3$pscore)   #treated pscore range
  maxT3=max(dataT3$pscore)
  
  minC3=min(dataC3$pscore)   #untreated pscore range
  maxC3=max(dataC3$pscore)
  
  data3$support= rep(0,nrow(data3))
  
  for (k in 1:K){
    if (data3$pscore[k]>=minT3 & data3$pscore[k]<=maxC3){
      data3$support[k]=1
    }
  }
  
  data$pscore=rep(0,len)
  data$commonsupport=rep(0,len)
  data$commonsupport[data1$id]=data1$support[1:length(data1$id)]
  data$pscore[data1$id]=data1$pscore[1:length(data1$id)]
  data$commonsupport[data2$id]=data2$support[1:length(data2$id)]
  data$pscore[data2$id]=data2$pscore[1:length(data2$id)]
  data$commonsupport[data3$id]=data3$support[1:length(data3$id)]
  data$pscore[data3$id]=data3$pscore[1:length(data3$id)]
  
  
  cs_data=subset(data, commonsupport==1)
  cs_data=subset(cs_data, pscore>=0.02 & pscore<=0.98)
  cs_len=nrow(cs_data)
  lenX=ncol(X)
  cs_X=as.data.frame(cbind(X, id))
  cs_X=cs_X[cs_data$id,1:lenX]
  
  
  cs_data$w_att=rep(1,cs_len)
  
  for (w in 1:cs_len){
    if (cs_data$T1[w]==0) {
      cs_data$w_att[w]=cs_data$pscore[w]/(1-cs_data$pscore[w])
    }
  }
  
  #trimming for pscore value close to 0 and to 1
  
  Xmatrix=as.matrix(cs_X)
  triplematchinglasso=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix,weights=cs_data$w_att,cluster="d")
  # version with non-clustered standard errors would be 
  # triplematching= lm(Y ~ t*d+t+d+Xmatrix, data = cs_data, weights=w_att) 
  return(triplematchinglasso)
}


#' @title Minimax balancing weights to utilise Lu, Nie and Wager (2019)
#' @description function that finds the minimax balancing weights, using the method in
#' Lu, Nie, Wager (2019)
#'
#' @param X the input features
#' @param Ti the time variable (0 or 1)
#' @param Si the state variable (0 or 1)
#' @param zeta tuning parameter, in [0,1]
#' @param solver cvx opmitization solver, choose from "ECOS" or "SCS"
#'
#' @return a vector gamma of the same length as Ti
#'
#' @export


minimax = function(X, Ti, Si, zeta=0.5, solver = c("ECOS", "SCS"), verbose = FALSE) {
  solver = match.arg(solver)
  nobs = nrow(X)
  pobs = ncol(X)
  gg = CVXR::Variable(nobs + 4)
  objective = (1 - zeta) * sum(gg[1:nobs]^2) + zeta * sum(gg[nobs + 1:4]^2)
  contraints = list(
    sum(gg[1:nobs]) == 0,
    t(X) %*% gg[1:nobs] <= gg[nobs + 1],
    -t(X) %*% gg[1:nobs] <= gg[nobs + 1],
    t(X) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
    -t(X) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
    t(X) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
    -t(X) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
    sum(Ti * Si * gg[1:nobs]) == 1,
    t(X) %*% (Ti * Si * gg[1:nobs]) <= colMeans(X) + gg[nobs + 4],
    - t(X) %*% (Ti * Si * gg[1:nobs]) <= - colMeans(X) + gg[nobs + 4]
  )
  cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
  cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
  result = cvx.output$getValue(gg)
  gamma = nobs * result[1:nobs]
  return(gamma)
}

#' @title generate hermite basis to utilise Lu, Nie and Wager (2019)
#' @description  Generates hermite basis for a given order, which we use for
#' non-parametric regression
#'
#' @param X an n by p matrix containing the input features
#' @param order order of hermite polynomials that we will use
#' @return a n by q matrix, which is a basis expansion of the input X
#' @export
generate.basis = function(X, order=3) {
  H = lapply(1:ncol(X), function(j) {
    sapply(1:order, function(k) hermite(X[,j], k, prob = TRUE) / sqrt(factorial(k)))
  })
  polys = lapply(1:order, function(r) {
    partitions = combn(r + ncol(X) -1, ncol(X) - 1,
                       function(vec) c(vec, r + ncol(X)) - c(0, vec) - 1)
    elems = sapply(1:ncol(partitions), function(iter) {
      part = partitions[,iter]
      idx = which(part > 0)
      elem = H[[idx[1]]][,part[idx[1]]]
      if (length(idx) > 1) {
        for (id in idx[-1]) {
          elem = elem * H[[id]][,part[id]]
        }
      }
      elem
    })
    scale(elems) / sqrt(ncol(elems)) / r
  })
  Reduce(cbind, polys)
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
M = 100           # Number of experiments/iterations




# Results Storage 

ols <- rep(0,M)
did <- rep(0,M)
ipwdid<- rep(0,M)
ordid<- rep(0,M)
impdrdid<- rep(0,M)
DiD_AMLE<- rep(0,M)
triplematch<- rep(0,M)
ATT<- rep(0,M)

# START OF THE SIMULATION LOOP

for (i in 1:M){ #  M is the number of iterations
  
  
  # create four heterogenous groups C0, C1, T0, T1
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  
  # Group 1: C0
  
  C0_i=rep(1, n/4)               #group dummy
  C1_i=rep(0, n/4)
  T0_i=rep(0, n/4)
  T1_i=rep(0, n/4)
  
  # create a multivariate distribution of the observed covariates:
  
  # 1.create a symmetric matrix of correlations among variables
  
  V<-matrix(c(
    1,.152,.096,.043,.109,
    .152,1,.400,-.016,.297,
    .096,.400,1,.092,.382,
    .043,-.016,.092,1,.103,
    .109,.297,.382,.103,1),
    5,5)
  
  # 2. create a vector of standard deviations for variables
  sigma<-c(0.5, 1, 4, 0.5, 1.5)
  
  # 3. Convert correlation matrix into covariance matrix Sigma
  Sigma<-cor2cov(V,sigma)
  
  # 4. vector of mean values
  mu<-c(1,3,10,2,5)
  
  # 5. generate multivariate distribution
  xdistr=mvrnorm(n=n/4,mu,Sigma,5,5)
  
  # 6. take the vectors of the covariates
  x1_i=xdistr[,1]
  x2_i=xdistr[,2]
  x3_i=xdistr[,3]
  x4_i=xdistr[,4]
  x5_i=xdistr[,5]
  
  # create time, treatment group and id parameter
  
  t_i= rep(0, n/4)            # time parameter
  d_i= rep(0, n/4)            # treatment group parameter
  id_i=1:(n/4)                # id specifier
  
  # put all the variables in a subset
  
  C0_data=as.data.frame(cbind(id_i,C0_i,C1_i, T0_i, T1_i, d_i, t_i, x1_i, x2_i, x3_i, x4_i,x5_i))
  
  # Group 2: C1
  
  C0_i=rep(0, n/4)               #group dummy
  C1_i=rep(1, n/4)
  T0_i=rep(0, n/4)
  T1_i=rep(0, n/4)
  
  # 1. create a symmetric matrix of correlations among variables
  V<-matrix(c(
    1,.152,.096,.043,.109,
    .152,1,.400,-.016,.297,
    .096,.400,1,.092,.382,
    .043,-.016,.092,1,.103,
    .109,.297,.382,.103,1),
    5,5)
  
  # 2. create a vector of standard deviations for variables
  sigma<-c(0.5, 1, 4, 0.5, 1.5)
  
  # 3. Convert correlation matrix into covariance matrix Sigma
  Sigma<-cor2cov(V,sigma)
  
  # 4. vector of mean values
  mu<-c(1.2,3.3,10.8,2.2,5.3)
  
  # 5. generate distribution
  xdistr=mvrnorm(n=n/4,mu,Sigma,5,5)
  
  # 6. take the covariate vectors 
  x1_i=xdistr[,1]
  x2_i=xdistr[,2]
  x3_i=xdistr[,3]
  x4_i=xdistr[,4]
  x5_i=xdistr[,5]
  
  
  # create time, treatment group and id parameter
  
  t_i= rep(1, n/4)            # time sampling
  d_i= rep(0, n/4)            # treatment group sampling
  id_i=(n/4+1):(n/2)          # id specifier
  
  # put all the variables in a subset
  
  C1_data=as.data.frame(cbind(id_i,C0_i,C1_i, T0_i, T1_i, d_i, t_i, x1_i, x2_i, x3_i, x4_i,x5_i))
  
  # Group 3: T0
  
  C0_i=rep(0, n/4)               #group dummy
  C1_i=rep(0, n/4)
  T0_i=rep(1, n/4)
  T1_i=rep(0, n/4)
  
  # create a multivariate distribution of the observed covariates:
  
  # 1.create a symmetric matrix of correlations among variables
  V<-matrix(c(
    1,.152,.096,.043,.109,
    .152,1,.400,-.016,.297,
    .096,.400,1,.092,.382,
    .043,-.016,.092,1,.103,
    .109,.297,.382,.103,1),
    5,5)
  #
  # 2.create a vector of standard deviations for variables
  sigma<-c(0.5, 1, 4, 0.5, 1.5)
  
  # 3.Convert correlation matrix into covariance matrix Sigma
  Sigma<-cor2cov(V,sigma)
  
  # 4.vector of mean values
  mu<-c(1.3,3.5,11.4,2.3,5.3)
  
  # generate distribution
  xdistr=mvrnorm(n=n/4,mu,Sigma,5,5)
  
  # take the vectors of the covariates
  x1_i=xdistr[,1]
  x2_i=xdistr[,2]
  x3_i=xdistr[,3]
  x4_i=xdistr[,4]
  x5_i=xdistr[,5]
  
  
  # create time, treatment group and id parameter
  
  t_i= rep(0, n/4)            # time sampling
  d_i= rep(1, n/4)            # treatment group sampling
  id_i=(n/2+1):(n/4*3)        # id specifier
  
  
  # put all the variables together
  
  T0_data=as.data.frame(cbind(id_i,C0_i,C1_i, T0_i, T1_i, d_i, t_i, x1_i, x2_i, x3_i, x4_i,x5_i))
  
  
  # Group 4: T1
  
  C0_i=rep(0, n/4)               #group dummy
  C1_i=rep(0, n/4)
  T0_i=rep(0, n/4)
  T1_i=rep(1, n/4)
  
  # create a multivariate distribution of the observed covariates:
  
  # 1.create a symmetric matrix of correlations among variables
  V<-matrix(c(
    1,.152,.096,.043,.109,
    .152,1,.400,-.016,.297,
    .096,.400,1,.092,.382,
    .043,-.016,.092,1,.103,
    .109,.297,.382,.103,1),
    5,5)
  
  # 2.create a vector of standard deviations for variables
  sigma<-c(0.5, 1, 4, 0.5, 1.5)
  
  # 3.Convert correlation matrix into covariance matrix Sigma
  Sigma<-cor2cov(V,sigma)
  
  # 4.vector of mean values
  mu<-c(1.5,3.7,12.4,2.5,5.6)
  
  # 5.generate distribution
  xdistr=mvrnorm(n=n/4,mu,Sigma,5,5)
  
  # 6.take the vectors of the covariates
  x1_i=xdistr[,1]
  x2_i=xdistr[,2]
  x3_i=xdistr[,3]
  x4_i=xdistr[,4]
  x5_i=xdistr[,5]
  
  
  # create time, treatment group and id parameter
  
  t_i= rep(1, n/4)            # time sampling
  d_i= rep(1, n/4)            # treatment group sampling
  id_i=(n/4*3+1):n            # id specifier
  
  #put all data together 
  
  T1_data=as.data.frame(cbind(id_i,C0_i,C1_i, T0_i, T1_i, d_i, t_i, x1_i, x2_i, x3_i, x4_i,x5_i))
  
  
  # append the four groups and redefinite variable for the whole sample
  
  data_i=rbind(C0_data,C1_data, T0_data, T1_data)
  
  x1_i=data_i$x1_i
  x2_i=data_i$x2_i
  x3_i=data_i$x3_i
  x4_i=data_i$x4_i
  x5_i=data_i$x5_i
  C0_i=data_i$C0_i
  C1_i=data_i$C1_i
  T0_i=data_i$T0_i
  T1_i=data_i$T1_i
  t_i=data_i$t_i
  d_i=data_i$d_i
  id_i=data_i$id_i
  
  # create a time-invariant unobservable characteristic
  # I create this unobservable correlated with the treatment to create bias when not using first differences.
  # To do that, I estimate a propensity score.
  # The estimated propensity score is highly correlated with the true propensity score,
  # thus I generate a random variable that is correlated with the estimated propensity score
  # and, therefore, also with the true one, and which is time invariant within the
  # treated and untreated group
  
  #propensity score estimation 
  
  est_pscore <- glm(T1_i ~ x1_i+x2_i+x3_i+x4_i+x5_i, family = binomial(link = "probit"), 
                    data = data_i)
  
  data_i$est_pscore=est_pscore$fitted.values
  
  # generate unobserved random variable which is correlated with the propensity score
  data_i$xunobs_i=rep(0,n)
  
  # for c0
  data_i$xunobs_i[1:(n/4)] = correlate(given = data_i$est_pscore[1:(n/4)], rho = 0.7,
                                draw_count, mean = 6)
  # for c1
  data_i$xunobs_i[(n/4+1):(n/2)] = correlate(given = data_i$est_pscore[(n/4+1):(n/2)], rho = 0.7,
                                       draw_count, mean = 7)
  
  # for T0
  data_i$xunobs_i[(n/2+1):(n/4*3)] = correlate(given = data_i$est_pscore[(n/2+1):(n/4*3)], rho = 0.7,
                                  draw_count, mean = 9)
  # for T1
  data_i$xunobs_i[(n/4*3+1):n] = correlate(given = data_i$est_pscore[(n/4*3+1):n], rho = 0.7,
                                         draw_count, mean = 10)
  
  xunobs_i=data_i$xunobs_i
  
  #drop est_pscore
  data_i=data_i[,-13]
  
  #create interaction matrix between covariates
  
  x_i=data_i[8:12]
  x_i[6:15]=as.data.frame(interact.data(data_i[8:12], base.cols = NULL, exclude.pair = NULL))
  x_i=as.matrix(x_i)
  
  #create the function determining potential and observed outcomes
  
  data_i$noise00_i = rnorm(n, mean = 0, sd = 5)   # random error
  data_i$noise01_i = rnorm(n, mean = 0, sd = 5)
  data_i$noise10_i = rnorm(n, mean = 0, sd = 5)
  data_i$noise11_i = rnorm(n, mean = 0, sd = 5)
  
  data_i$Y00_i=rep(0,n)                           #Ydt potential outcomes
  data_i$Y01_i=rep(0,n)
  data_i$Y10_i=rep(0,n)
  data_i$Y11_i=rep(0,n)
  
  # covariates matrix multiplied by effects
  # multiply the value of the covariate by an arbitrary effect
  
  effects=matrix(c(6, 4, 1.5, 2, 2, 0.35, 0.175, 0.205, 0.195, 0.085, 0.190 ,0.075, 0.075, 0.035, 0.125 ),
                 nrow = 1, ncol=15)
  effects=rep.row(effects, n)
  xe_i=effects*x_i
  
  #parameters 
  intercept = 3                         # Intercept 
  time_eff = 0.5*x1_i+x2_i+0.8*x4_i     # time effect
  treat_group= 0.4*x3_i                 # treatment group effect
  D=0.2*(x1_i+x2_i+0.3*x5_i)^2          # individual treatment effect 
  ATT[i]=mean(D[(n/4*3+1):n])           # ATT (here we are assuming heterogenous effects)
  
  
  # create potential outcomes
  # Y00 is the potential outcome with no-treatment at time 0
  # Y01 is the potential outcome with no-treatment at time 1
  # Y00 is the potential outcome with treatment at time 0
  # Y00 is the potential outcome with treatment at time 1
  
  for (b in 1:n){
    data_i$Y00_i[b] = intercept +2*sin(sum(xe_i[b,1:2]))+0.1*sum(xe_i[b,3:15])^2+xunobs_i[b]+data_i$noise00_i[b]
    data_i$Y01_i[b] = intercept +2*sin(sum(xe_i[b,1:2]))+0.1*sum(xe_i[b,3:15])^2+xunobs_i[b]+time_eff[b]+data_i$noise01_i[b]
    data_i$Y10_i[b] = intercept +2*sin(sum(xe_i[b,1:2]))+0.1*sum(xe_i[b,3:15])^2+xunobs_i[b]+treat_group[b]+data_i$noise10_i[b]
    data_i$Y11_i[b] = intercept +2*sin(sum(xe_i[b,1:2]))+0.1*sum(xe_i[b,3:15])^2+xunobs_i[b]+time_eff[b]+treat_group[b]+D[b]+data_i$noise11_i[b]
  }
  
  
  
  # create the observed outcome 
  
  y_i=rep(0, n)
  
  for (s in 1:n){
    if (data_i$t_i[s]==0 & data_i$d_i[s]==0) {
      data_i$y_i[s]=data_i$Y00_i[s]
    } else if (data_i$t_i[s]==1 & data_i$d_i[s]==0) {
      data_i$y_i[s]=data_i$Y01_i[s]
    } else if (data_i$t_i[s]==0 & data_i$d_i[s]==1) {
      data_i$y_i[s]=data_i$Y10_i[s]
    } else {
      data_i$y_i[s]=data_i$Y11_i[s]
    }
  }
  
  y_i=data_i$y_i
  
  #### the simulated dataset is now ready for the estimation ###
  
  #Naive OLS
  
  ols_i <- lm(y_i ~ x_i+d_i, data=data_i)
  ols[i] <- ols_i$coefficients["d_i"]
  
  #Standard difference in difference model
  
  did_i <- lm(y_i ~ x_i+t_i+d_i+t_i*d_i, data=data_i)
  
  # clustered standard errors?, but we should not be concerned about variance in our simulation 
  # mod1 <- miceadds::lm.cluster( data=data_i, formula=y_i ~ x_i+t_i+d_i+t_i*d_i,cluster="d_i")
  # Extract slope coefficient and save
  
  did[i] <- did_i$coefficients["t_i:d_i"]
  
  # Implement DID ipw Abadie(2005) with normalized weights
  
  ipwdid_i=ipwdid(yname="y_i", tname = "t_i", idname = "id_i",  dname = "d_i",
                  xformla= ~ x1_i+x2_i+x3_i+x4_i +x5_i+x1_i*x2_i+x1_i*x3_i+ x1_i*x4_i+ x1_i*x5_i+
                    x2_i*x3_i+ x2_i*x4_i+ x2_i*x5_i+ x3_i*x4_i+ x3_i*x5_i+ x4_i*x5_i,
                  data = data_i, panel = FALSE,
                  boot = TRUE, nboot = 199)
  ipwdid[i] <- ipwdid_i$ATT
  
  # Implement Outcome regression model following Heckman, Ichimura and Todd (1997) but with linear assumption of covariates
  
  ordid_i=ordid(yname="y_i", tname = "t_i", idname = "id_i", dname = "d_i",
                xformla= ~ x1_i+x2_i+x3_i+x4_i +x5_i+x1_i*x2_i+x1_i*x3_i+ x1_i*x4_i+ x1_i*x5_i+
                  x2_i*x3_i+ x2_i*x4_i+ x2_i*x5_i+ x3_i*x4_i+ x3_i*x5_i+ x4_i*x5_i,
                data = data_i, panel = FALSE,
                boot = TRUE, nboot = 199)
  ordid[i] <- ordid_i$ATT
  
  # Implement improved locally efficient DR DID (Sant'Anna and Zhao,2020):
  
  impdrdid_i <- drdid(yname = "y_i", tname = "t_i", idname = "id_i", dname = "d_i",
                      xformla= ~ x1_i+x2_i+x3_i+x4_i +x5_i+x1_i*x2_i+x1_i*x3_i+ x1_i*x4_i+ x1_i*x5_i+
                        x2_i*x3_i+ x2_i*x4_i+ x2_i*x5_i+ x3_i*x4_i+ x3_i*x5_i+ x4_i*x5_i,
                      data = data_i, panel = FALSE, estMethod = "imp")
  impdrdid[i] <- impdrdid_i$ATT
  
  
  # Implement Lu, Nie and Wager (2019)
  
  p=ncol(x_i)
  X = x_i[,1:5]
  tau = D
  f = time_eff
  h = treat_group
  Si = data_i$d_i
  Ti = data_i$t_i
  Y = data_i$y_i
  
  # generate a basis
  make_matrix = function(x) stats::model.matrix(~.-1, x)
  X = data.frame(X) %>% make_matrix
  order = 3
  Basis = generate.basis(X,order)
  
  gamma.minimax = minimax(Basis, Ti, Si)
  DiD_AMLE_i = DiD(X = Basis,
                   Y = Y,
                   Ti = Ti,
                   Si = Si,
                   constant_eff = "non_constant",
                   gamma = gamma.minimax)
  
  DiD_AMLE[i] <- DiD_AMLE_i$"TAU_hat"
  
  
  # Implement triple matching following Bludell(2004), propensity score estimated by probit regression
  
  X=as.data.frame(x_i)
  T0=data_i$T0_i
  T1=data_i$T1_i
  C0=data_i$C0_i
  C1=data_i$C1_i
  id=data_i$id_i
  t=data_i[,7]
  d=data_i[,6]
  Y=data_i$y_i
  
  
  triplematch_i=triplematching(id, T0, T1, C0, C1, t, d, X, Y)
  triplematch[i] <- coef(triplematch_i)["t:d"]
  
}


#generate summary statistics for each estimator

ols_mean=mean(ols)
ols_variance=mean((ols-mean(ols))^2)
ols_mse=mean((ols-ATT)^2)
ols_bias=mean (mean(ols)-ATT)

did_mean=mean(did)
did_variance=mean((did-mean(did))^2)
did_mse=mean((did-ATT)^2)
did_bias=mean (mean(did)-ATT)

ipwdid_mean=mean(ipwdid)
ipwdid_variance=mean((ipwdid-mean(ipwdid))^2)
ipwdid_mse=mean((ipwdid-ATT)^2)
ipwdid_bias=mean (mean(ipwdid)-ATT)

ordid_mean=mean(ordid)
ordid_variance=mean((ordid-mean(ordid))^2)
ordid_mse=mean((ordid-ATT)^2)
ordid_bias=mean (mean(ordid)-ATT)

impdrdid_mean=mean(impdrdid)
impdrdid_variance=mean((impdrdid-mean(impdrdid))^2)
impdrdid_mse=mean((impdrdid-ATT)^2)
impdrdid_bias=mean (mean(impdrdid)-ATT)

triplematching_mean=mean(triplematch)
triplematching_variance=mean((triplematch-mean(triplematch))^2)
triplematching_mse=mean((triplematch-ATT)^2)
triplematching_bias=mean (mean(triplematch)-ATT)

DiD_AMLE_mean=mean(DiD_AMLE)
DiD_AMLE_variance=mean((DiD_AMLE-mean(DiD_AMLE))^2)
DiD_AMLE_mse=mean((DiD_AMLE-ATT)^2)
DiD_AMLE_bias=mean (mean(DiD_AMLE)-ATT)

#construct a summary table

mean=cbind(ols_mean, did_mean, ipwdid_mean, ordid_mean, impdrdid_mean, triplematching_mean, DiD_AMLE_mean)
bias=cbind(ols_bias,did_bias, ipwdid_bias, ordid_bias, impdrdid_bias, triplematching_bias, DiD_AMLE_bias)
mse=cbind(ols_mse,did_mse, ipwdid_mse, ordid_mse, impdrdid_mse, triplematching_mse, DiD_AMLE_mse)
variance=cbind(ols_variance, did_variance, ipwdid_variance, ordid_variance, impdrdid_variance, triplematching_variance, DiD_AMLE_variance)

tab <- matrix(1:21, ncol=3, byrow=TRUE)
tab[,1]=bias
tab[,2]=mse
tab[,3]=variance



colnames(tab) <- c('bias','MSE', 'variance')
rownames(tab) <- c('ols','did', 'ipwdid','ordid', 'impdrdid', 'triplematching', 'DiD_AMLE')
tab <- as.table(tab)
tab