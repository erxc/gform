################################################################################
rm(list = ls())
options(warn=-1)

# balancing score-based g-formula

library(geepack)
library(MASS)
library(ResourceSelection)
library(ltmle)
library(SuperLearner)
library(dplyr)
library(pillar)
library(data.table)
library(splines)
library(BART)
library(Hmisc)

source("datagen_dd_sim.R")

set.seed(0919)
seeds <- floor(runif(1000)*10^8)

bart_bs_gform_func <- function(m){
  
  set.seed(seeds[m])
  
  rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
  logit <- function(x) (log(x) - log(1-x))
  
  # set up and parameters
  n <- 500 # number of participants
  K <- 3 # number of observation periods          
  
  # exposure model
  ## 40% censoring 
  # eta0 <- -2; eta1 <- -1; eta2 <- 0.75; eta3 <- -0.5; eta4 <- 0.2
  ## 20% censoring
  eta0 <- -3; eta1 <- -1; eta2 <- 0.75; eta3 <- -0.5; eta4 <- 0.2
  
  # covariate model
  beta0 <- 0; beta1 <- -2; beta2 <- 0.2; beta3 <- 1; beta4 <- 0.3
  # continuous covariate sd
  sigma <- 0.1
  
  # outcome model
  theta0 <- 0; theta1 <- -2; theta2 <- 1.25; theta3 <- 1.5; theta4 <- -1.4
  
  # generate data for all 500 patients
  df <- lapply(
    as.list(1:n),
    FUN = function(ind) {
      datagen_dd_nonlinear(
        ind,
        K = K,
        beta0 = beta0,
        beta1 = beta1,
        beta2 = beta2,
        beta3 = beta3,
        beta4 = beta4,
        theta0 = theta0,
        theta1 = theta1,
        theta2 = theta2,
        theta3 = theta3,
        theta4 = theta4,
        eta0 = eta0,
        eta1 = eta1,
        eta2 = eta2,
        eta3 = eta3,
        eta4 = eta4,
        sigma = sigma
      )
    }
  )
  
  # combine and create lag one observations - this is because we have assumed first order Markovian
  dffull <- rbindlist(df)
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_C") := shift(C, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L1") := shift(L1, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L2") := shift(L2, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L3") := shift(L3, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L4") := shift(L4, 1, NA, type='lag'), by=id]
  
  # Fit propensity score model and censoring score model
  nt <- 200
  # this is for the generative model of A 
  ## this is to find the X's for the generative model of A
  X_amodel <- dffull %>%
    select(lag_A, L1, L2, L3, L4) %>%
    as.data.frame()
  
  ## this is to find the Y's for the generative model of A
  Y_amodel <- dffull$A
  
  ## this is to fit the generative model of A
  invisible(capture.output(bart_amodel <-
    lbart(X_amodel, Y_amodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # this is for the generative model of C 
  ## this is to find the X's for the generative model of C
  X_cmodel <- dffull %>%
    select(A, L1, L2, L3, L4) %>% 
    as.data.frame()
  
  ## this is to find the Y's for the generative model of C
  Y_cmodel <- dffull$C  
  
  ## this is to fit the generative model of C
  invisible(capture.output(bart_cmodel <-
    lbart(X_cmodel, Y_cmodel, nskip = 1000, ndpost = 10000, ntree = nt, 
         nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # estimate propensity scores
  ## estimate zeta
  dffull$pred_obs <- bart_amodel$prob.train.mean
  ## estimate omega
  dffull$pred_obsc <- 1 - bart_cmodel$prob.train.mean
  ## overall ps
  dffull$E <- (dffull$pred_obs*dffull$A + (1-dffull$pred_obs)*(1-dffull$A)) * 
    dffull$pred_obsc  
  ## overall logit ps
  dffull$logE <- logit(dffull$E)                                                         
  # create lags
  dffull[, paste("lag_E") := shift(E, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_logE") := shift(logE, 1, NA, type='lag'), by=id]
  
  ## pooled ps model (these are flexible models, they are agnostic to the true models)
  X_emodel <- dffull[!is.na(dffull$lag_logE),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_logE, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  
  Y_emodel <- dffull[!is.na(dffull$lag_logE),]$logE
  
  invisible(capture.output(bart_emodel <- 
    wbart(X_emodel, Y_emodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptestmean = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # this is for the generative model of L2
  ## this is to find the X's for the generative model of L2
  ### it requires logE and L2 to be non-NA
  X_covmodel2 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_logE, lag_L2, logE, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)

  ## this is to find the Y's for the generative model of L2
  Y_covmodel2 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2),]$L2

  ## this is to fit the generative model of L2
  invisible(capture.output(bart_covmodel2 <-
    wbart(X_covmodel2, Y_covmodel2, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # the generative model of L3
  X_covmodel3 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_logE, lag_L2, lag_L3, logE, L2, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  Y_covmodel3 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3),]$L3
  
  invisible(capture.output(bart_covmodel3 <-
    wbart(X_covmodel3, Y_covmodel3, nskip = 1000, ndpost = 10000, ntree = nt,
         nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # the generative model of L4
  X_covmodel4 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3) & !is.na(dffull$lag_L4),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_logE, lag_L2, lag_L3, lag_L4, logE, L2, L3, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  ## this is to find the Y's for the generative model of L4
  Y_covmodel4 <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3) & !is.na(dffull$lag_L4),]$L4
  
  ## this is to fit the generative model of L4
  invisible(capture.output(bart_covmodel4 <-
    wbart(X_covmodel4, Y_covmodel4, nskip = 1000, ndpost = 10000, ntree = nt,
         nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # pooled Y model
  dffull_exclude_NA_Y <- dffull[(!is.na(dffull$Y)),]
  Y_ymodel <- dffull_exclude_NA_Y$Y
  X_ymodel <- dffull_exclude_NA_Y %>%
    mutate(t0 = as.factor(t0)) %>%
    select(A, logE, L2, L3, L4, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  invisible(capture.output(bart_ymodel <-
    lbart(X_ymodel, Y_ymodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  sigma_E <- mean(bart_emodel$sigma)
  sigma_L2 <- mean(bart_covmodel2$sigma)
  sigma_L3 <- mean(bart_covmodel3$sigma)
  sigma_L4 <- mean(bart_covmodel4$sigma)
  
  ## Monte Carlo simulation
  N <- 3000
  ids <- as.list(1:N)
  result_bart <- matrix(NA, nrow = N, ncol = K) # storing results in all five periods
  
  for (i in 1:N){
    logE0 <- L20 <- L30 <- L40 <- A0 <- Y1 <-
      logE1 <- L21 <- L31 <- L41 <- A1 <- Y2 <-
      logE2 <- L22 <- L32 <- L42 <- A2 <- Y3 <- as.numeric(rep(NA, N))
    
    # baseline
    ## bootstrapping the distribution of covariates
    dfbaseline <- dffull[dffull$t0 == 0,]
    idboot <- sample(1:nrow(dfbaseline), N, replace=T)
    uY <- rep(TRUE, N)
    logE0 <- dfbaseline$logE[idboot]
    L20 <- dfbaseline$L2[idboot]
    L30 <- dfbaseline$L3[idboot]
    L40 <- dfbaseline$L4[idboot]
    A0 <- as.numeric((L20*L30 + cos(L40)) > 1)
    
    ## samples from BART posteriors of the generative model of Y1, probabilities from lbart
    invisible(capture.output(bart_predict_Y1_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A0, logE = logE0, L2 = L20, L3 = L30, L4 = L40,
                   t00 = 1, t01 = 0, t02 = 0))$prob.test.mean))
    ## generate Y1
    Y1 <- rbinom(N, 1, bart_predict_Y1_prob)
    ## remaining individuals
    uY <- uY & !Y1 
    
    # time 1
    ## generate logE1
    invisible(capture.output(bart_predict_logE1_mean <-
      predict(bart_emodel, newdata =
        data.frame(lag_A = A0, lag_logE = logE0, t00 = 1, t01 = 0))[1,]))
    logE1[uY] <- rnorm(sum(uY), mean = bart_predict_logE1_mean[uY], sd = sigma_E)
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t1_L2_mean <-
      predict(bart_covmodel2, newdata = 
                data.frame(lag_A = A0, lag_logE = logE0, lag_L2 = L20, logE = logE1,
                           t00 = 1, t01 = 0))[1,]))
    L21[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t1_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t1_L3_mean <-
      predict(bart_covmodel3, newdata = 
                data.frame(lag_A = A0, lag_logE = logE0, lag_L2 = L20, 
                           lag_L3 = L30, logE = logE1, L2 = L21, 
                           t00 = 1, t01 = 0))[1,]))
    L31[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t1_L3_mean[uY], sd = sigma_L3)
    
    ## generate L4
    invisible(capture.output(bart_covmodel4_predict_t1_L4_mean <-
      predict(bart_covmodel4, newdata = 
                data.frame(lag_A = A0, lag_logE = logE0, lag_L2 = L20, 
                           lag_L3 = L30, lag_L4 = L40, logE = logE1, 
                           L2 = L21, L3 = L31, t00 = 1, t01 = 0))[1,]))
    L41[uY] <- rnorm(sum(uY), mean = bart_covmodel4_predict_t1_L4_mean[uY], sd = sigma_L4)
    
    ## generate A1
    A1[uY] <- as.numeric(((L21[uY]*L31[uY] + cos(L41[uY])) > 1) | A0[uY])

    ## generate Y2
    invisible(capture.output(bart_predict_Y2_prob <-
      predict(bart_ymodel, newdata = 
        data.frame(A = A1, logE = logE1, L2 = L21, L3 = L31, L4 = L41,
                   t00 = 0, t01 = 1, t02 = 0))$prob.test.mean))    
    Y2[uY] <- rbinom(sum(uY), 1, bart_predict_Y2_prob[uY])
    ## remaining individuals
    uY <- uY & !Y2 
    
    # time 2
    ## generate logE2
    invisible(capture.output(bart_predict_logE2_mean <-
      predict(bart_emodel, newdata =
                data.frame(lag_A = A1, lag_logE = logE1, t00 = 0, t01 = 1))[1,]))
    logE2[uY] <- rnorm(sum(uY), mean = bart_predict_logE2_mean[uY], sd = sigma_E)
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t2_L2_mean <-
      predict(bart_covmodel2, newdata = 
               data.frame(lag_A = A1, lag_logE = logE1, lag_L2 = L21, logE = logE2, 
                          t00 = 0, t01 = 1))[1,]))
    L22[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t2_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t2_L3_mean <-
      predict(bart_covmodel3, newdata = 
               data.frame(lag_A = A1, lag_logE = logE1, lag_L2 = L21, 
                          lag_L3 = L31, logE = logE2, L2 = L22, 
                          t00 = 0, t01 = 1))[1,]))
    L32[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t2_L3_mean[uY], sd = sigma_L3)
    
    ## generate L4
    invisible(capture.output(bart_covmodel4_predict_t2_L4_mean <-
      predict(bart_covmodel4, newdata = 
               data.frame(lag_A = A1, lag_logE = logE1, lag_L2 = L21, 
                          lag_L3 = L31, lag_L4 = L41, logE = logE2, 
                          L2 = L22, L3 = L32, t00 = 0, t01 = 1))[1,]))
    L42[uY] <- rnorm(sum(uY), mean = bart_covmodel4_predict_t2_L4_mean[uY], sd = sigma_L4)
    
    ## generate A2
    A2[uY] <- as.numeric(((L22[uY]*L32[uY] + cos(L42[uY])) > 1) | A1[uY])
    
    ## generate Y3
    invisible(capture.output(bart_predict_Y3_prob <-
      predict(bart_ymodel, newdata = 
        data.frame(A = A2, logE = logE2, L2 = L22, L3 = L32, L4 = L42,
                   t00 = 0, t01 = 0, t02 = 1))$prob.test.mean))    
    Y3[uY] <- rbinom(sum(uY), 1, bart_predict_Y3_prob[uY])
    
    # compute estimates
    tmpdata <- data.frame(Y1, Y2, Y3)
    
    tmpdata$Y2 <- ifelse(tmpdata$Y1 == 1, 1, tmpdata$Y2)
    tmpdata$Y3 <- ifelse(!is.na(tmpdata$Y2) & tmpdata$Y2==1, 1, tmpdata$Y3)
    
    result_bart[i,] <- colMeans(tmpdata)
    
    if (i %% 100 == 0) {
      cat("\r", paste("MC sample", i, sep=" "))
      flush.console()
    }

  }
  output <- colMeans(result_bart)
  
  return(output)
}

bart_bs_gform_func(1)
  

