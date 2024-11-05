################################################################################
rm(list = ls())
options(warn=-1)

# propensity score-based g-formula

library(geepack)
library(MASS)
library(ResourceSelection)
library(ltmle)
library(SuperLearner)
library(dplyr)
library(data.table)
library(splines)
library(BART)
library(Hmisc)

source("datagen_dd.R")

bart_cov_bs_gform_func <- function(m){
  
  set.seed(seeds[m])
  
  rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
  logit <- function(x) (log(x) - log(1-x))
  
  # set up and parameters
  n <- 500 # number of participants
  K <- 5 # number of observation periods          
  
  # exposure model
  eta0 <- -3; eta1 <- -1; eta2 <- 0.75; eta3 <- -0.5; eta4 <- 0 # 20% censoring
  
  # covariate model
  beta0 <- 0; beta1 <- -2; beta2 <- 0.2; beta3 <- 1
  # continuous covariate sd
  sigma <- 0.1 
  
  # outcome model
  theta0 <- -2; theta1 <- -3; theta2 <- 1; theta3 <- -6; theta4 <- 0
  
  # generate data for all 1000 patients
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
  
  # sum(unlist(lapply(df, function(x) 1*(sum(x$A) > 0))))
  
  # combine and create lag one observations - this is because we have assumed first order Markovian
  dffull <- rbindlist(df)
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_C") := shift(C, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L1") := shift(L1, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L2") := shift(L2, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L3") := shift(L3, 1, NA, type='lag'), by=id]
  
  # Fit propensity score model and censoring score model
  nt <- 200
  # this is for the generative model of A 
  ## this is to find the X's for the generative model of A
  X_amodel <- dffull %>%
    select(lag_A, L1, L2, L3) %>%
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
    select(A, L1, L2, L3) %>% 
    as.data.frame() 
  
  ## this is to find the Y's for the generative model of C
  Y_cmodel <- dffull$C  
  
  ## this is to fit the generative model of C
  invisible(capture.output(bart_cmodel <-
    lbart(X_cmodel, Y_cmodel, nskip = 1000, ndpost = 10000, ntree = nt, 
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # estimate propensity scores for static, treat-every-other-period strategy
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
  
  # this is for the generative model of L1
  ## this is to find the X's for the generative model of L1
  ### it requires L1 to be non-NA b/c L1 depends on L1
  X_covmodel1 <- dffull[!is.na(dffull$lag_L1),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  ## this is to find the Y's for the generative model of L1
  ### it requires L1 to be non-NA b/c L1 depends on L1
  Y_covmodel1 <- dffull[!is.na(dffull$lag_L1),]$L1
  
  ## this is to fit the generative model of L1
  invisible(capture.output(bart_covmodel1 <- 
    lbart(X_covmodel1, Y_covmodel1, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # this is for the generative model of L2 
  ## this is to find the X's for the generative model of L2
  ### it requires both L1 and L2 to be non-NA b/c L2 depends on both L1 and L2
  X_covmodel2 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, lag_L2, L1, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  ## this is to find the Y's for the generative model of L2
  ### it requires L2 to be non-NA b/c L2 depends on L2
  Y_covmodel2 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2),]$L2
  
  ## this is to fit the generative model of L2
  invisible(capture.output(bart_covmodel2 <- 
    wbart(X_covmodel2, Y_covmodel2, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptestmean = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # this is for the generative model of L3
  ## this is to find the X's for the generative model of L3
  ### it requires L1, L2 and L3 to be non-NA b/c L3 depends on L1, L2 and L3
  X_covmodel3 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2)  & !is.na(dffull$lag_L3),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, lag_L2, lag_L3, L1, L2, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  ## this is to find the Y's for the generative model of L3
  ### it requires L1, L2 and L3 to be non-NA b/c L3 depends on L1, L2 and L3
  Y_covmodel3 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2) & !is.na(dffull$lag_L3),]$L3
  
  ## this is to fit the generative model of L3
  invisible(capture.output(bart_covmodel3 <- 
    wbart(X_covmodel3, Y_covmodel3, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptestmean = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  ## pooled ps model (these are flexible models, they are agnostic to the true models)
  X_emodel <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L1) & 
                       !is.na(dffull$lag_L2) & !is.na(dffull$lag_L3),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, lag_L2, lag_L3, lag_logE, L1, L2, L3, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  Y_emodel <- dffull[!is.na(dffull$lag_logE) & !is.na(dffull$lag_L1) & 
                       !is.na(dffull$lag_L2) & !is.na(dffull$lag_L3),]$logE
  
  invisible(capture.output(bart_emodel <- 
    wbart(X_emodel, Y_emodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptestmean = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
 
  # pooled Y model
  dffull_exclude_NA_Y <- dffull[(!is.na(dffull$Y)),]
  Y_ymodel <- dffull_exclude_NA_Y$Y
  X_ymodel <- dffull_exclude_NA_Y %>%
    mutate(t0 = as.factor(t0)) %>%
    select(A, L1, L2, L3, logE, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  invisible(capture.output(bart_ymodel <-
    lbart(X_ymodel, Y_ymodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  sigma_L2 <- mean(bart_covmodel2$sigma)
  sigma_L3 <- mean(bart_covmodel3$sigma)
  sigma_E <- mean(bart_emodel$sigma)
  
  ## Monte Carlo simulation
  N <- 3000
  ids <- as.list(1:N)
  result_bart <- matrix(NA, nrow = N, ncol = K) # storing results in all five periods
  
  for (i in 1:N){
    L10 <- L20 <- L30 <- logE0 <- A0 <- Y1 <-
      L11 <- L21 <- L31 <- logE1 <- A1 <- Y2 <-
      L12 <- L22 <- L32 <- logE2 <- A2 <- Y3 <-
      L13 <- L23 <- L33 <- logE3 <- A3 <- Y4 <-
      L14 <- L24 <- L34 <- logE4 <- A4 <- Y5 <- as.numeric(rep(NA, N))
    
    # baseline
    ## bootstrapping the distribution of covariates
    dfbaseline <- dffull[dffull$t0 == 0,]
    idboot <- sample(1:nrow(dfbaseline), N, replace=T)
    uY <- rep(TRUE, N)
    L10 <- dfbaseline$L1[idboot]
    L20 <- dfbaseline$L2[idboot]
    L30 <- dfbaseline$L3[idboot]
    logE0 <- dfbaseline$logE[idboot]
    A0 <- as.numeric(L20 > 0.2)
    
    ## samples from BART posteriors of the generative model of Y1, probabilities from lbart
    invisible(capture.output(bart_predict_Y1_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A0, L1 = L10, L2 = L20, L3 = L30, logE = logE0, t00 = 1, t01 = 0, t02 = 0,
                   t03 = 0, t04 = 0))$prob.test.mean))
    ## generate Y1
    Y1 <- rbinom(N, 1, bart_predict_Y1_prob)
    ## remaining individuals
    uY <- uY & !Y1 
    
    # time 1
    ## generate L1
    invisible(capture.output(bart_predict_t1_L1_prob <-
      predict(bart_covmodel1, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, t01 = 1, t02 = 0, t03 = 0, t04 = 0))$prob.test.mean))
    L11[uY] <- rbinom(sum(uY), 1, bart_predict_t1_L1_prob[uY])
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t1_L2_mean <-
      predict(bart_covmodel2, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, L1 = L11,  
                   t01 = 1, t02 = 0, t03 = 0, t04 = 0))[1,]))
    L21[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t1_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t1_L3_mean <-
      predict(bart_covmodel3, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, lag_L3 = L30, L1 = L11, L2 = L21,
                   t01 = 1, t02 = 0, t03 = 0, t04 = 0))[1,]))
    L31[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t1_L3_mean[uY], sd = sigma_L3)
    
    ## generate logE
    invisible(capture.output(bart_predict_logE1_mean <-
      predict(bart_emodel, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, lag_L3 = L30, lag_logE = logE0, 
                   L1 = L11, L2 = L21, L3 = L31, t01 = 1, t02 = 0, t03 = 0, t04 = 0))[1,]))
    logE1[uY] <- rnorm(sum(uY), mean = bart_predict_logE1_mean[uY], sd = sigma_E)
    
    A1[uY] <- as.numeric((L21[uY] > 0.2) | A0[uY])
    
    invisible(capture.output(bart_predict_Y2_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A1, L1 = L11, L2 = L21, L3 = L31, logE = logE1, t00 = 0, t01 = 1, t02 = 0,
                   t03 = 0, t04 = 0))$prob.test.mean))
    ## generate Y2
    Y2[uY] <- rbinom(sum(uY), 1, bart_predict_Y2_prob[uY])
    ## remaining individuals
    uY <- uY & !Y2
    
    # time 2
    ## generate L1
    invisible(capture.output(bart_predict_t2_L1_prob <-
      predict(bart_covmodel1, newdata =
              data.frame(lag_A = A1, lag_L1 = L11, t01 = 0, t02 = 1, t03 = 0, t04 = 0))$prob.test.mean))
    L12[uY] <- rbinom(sum(uY), 1, bart_predict_t2_L1_prob[uY])
 
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t2_L2_mean <-
      predict(bart_covmodel2, newdata =
              data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, L1 = L12,  
                         t01 = 0, t02 = 1, t03 = 0, t04 = 0))[1,]))
    L22[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t2_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t2_L3_mean <-
      predict(bart_covmodel3, newdata =
              data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, lag_L3 = L31, L1 = L12, L2 = L22,
                         t01 = 0, t02 = 1, t03 = 0, t04 = 0))[1,]))
    L32[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t2_L3_mean[uY], sd = sigma_L3)
    
    ## generate logE
    invisible(capture.output(bart_predict_logE2_mean <-
      predict(bart_emodel, newdata =
        data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, lag_L3 = L31, lag_logE = logE1, 
                   L1 = L12, L2 = L22, L3 = L32, t01 = 0, t02 = 1, t03 = 0, t04 = 0))[1,]))
    logE2[uY] <- rnorm(sum(uY), mean = bart_predict_logE2_mean[uY], sd = sigma_E)
    
    A2[uY] <- as.numeric((L22[uY] > 0.2) | A1[uY])
    
    ## generate Y
    invisible(capture.output(bart_predict_Y3_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A2, L1 = L12, L2 = L22, L3 = L32, logE = logE2, t00 = 0, t01 = 0, t02 = 1,
                   t03 = 0, t04 = 0))$prob.test.mean))
    ## generate Y3
    Y3[uY] <- rbinom(sum(uY), 1, bart_predict_Y3_prob[uY])
    ## remaining individuals
    uY <- uY & !Y3
    
    # time 3
    ## generate L1
    invisible(capture.output(bart_predict_t3_L1_prob <-
      predict(bart_covmodel1, newdata =
              data.frame(lag_A = A2, lag_L1 = L12, t01 = 0, t02 = 0, t03 = 1, t04 = 0))$prob.test.mean))
    L13[uY] <- rbinom(sum(uY), 1, bart_predict_t3_L1_prob[uY])
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t3_L2_mean <-
      predict(bart_covmodel2, newdata =
              data.frame(lag_A = A2, lag_L1 = L12, lag_L2 = L22, L1 = L13,  
                         t01 = 0, t02 = 0, t03 = 1, t04 = 0))[1,]))
    L23[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t3_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t3_L3_mean <-
      predict(bart_covmodel3, newdata =
              data.frame(lag_A = A2, lag_L1 = L12, lag_L2 = L22, lag_L3 = L32, L1 = L13, L2 = L23,
                         t01 = 0, t02 = 0, t03 = 1, t04 = 0))[1,]))
    L33[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t3_L3_mean[uY], sd = sigma_L3)
    
    ## generate logE
    invisible(capture.output(bart_predict_logE3_mean <-
      predict(bart_emodel, newdata =
        data.frame(lag_A = A2, lag_L1 = L12, lag_L2 = L22, lag_L3 = L32, lag_logE = logE2, 
                   L1 = L13, L2 = L23, L3 = L33, t01 = 0, t02 = 0, t03 = 1, t04 = 0))[1,]))
    logE3[uY] <- rnorm(sum(uY), mean = bart_predict_logE3_mean[uY], sd = sigma_E)
    
    A3[uY] <- as.numeric((L23[uY] > 0.2) | A2[uY])
    
    invisible(capture.output(bart_predict_Y4_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A3, L1 = L13, L2 = L23, L3 = L33, logE = logE3, t00 = 0, t01 = 0, t02 = 0,
                   t03 = 1, t04 = 0))$prob.test.mean))
    
    ## generate Y4
    Y4[uY] <- rbinom(sum(uY), 1, bart_predict_Y4_prob[uY])
    ## remaining individuals
    uY <- uY & !Y4    

    # time 4
    ## generate L1
    invisible(capture.output(bart_predict_t4_L1_prob <-
      predict(bart_covmodel1, newdata =
              data.frame(lag_A = A3, lag_L1 = L13, t01 = 0, t02 = 0, t03 = 0, t04 = 1))$prob.test.mean))
    L14[uY] <- rbinom(sum(uY), 1, bart_predict_t4_L1_prob[uY])
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t4_L2_mean <-
      predict(bart_covmodel2, newdata =
              data.frame(lag_A = A3, lag_L1 = L13, lag_L2 = L23, L1 = L14,  
                         t01 = 0, t02 = 0, t03 = 0, t04 = 1))[1,]))
    L24[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t4_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t4_L3_mean <-
      predict(bart_covmodel3, newdata =
              data.frame(lag_A = A3, lag_L1 = L13, lag_L2 = L23, lag_L3 = L33, L1 = L14, L2 = L24,
                         t01 = 0, t02 = 0, t03 = 0, t04 = 1))[1,]))
    L34[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t4_L3_mean[uY], sd = sigma_L3)
    
    ## generate logE
    invisible(capture.output(bart_predict_logE4_mean <-
      predict(bart_emodel, newdata =
        data.frame(lag_A = A3, lag_L1 = L13, lag_L2 = L23, lag_L3 = L33, lag_logE = logE3, 
                   L1 = L14, L2 = L24, L3 = L34, t01 = 0, t02 = 0, t03 = 0, t04 = 1))[1,]))
    logE4[uY] <- rnorm(sum(uY), mean = bart_predict_logE4_mean[uY], sd = sigma_E)
    
    A4[uY] <- as.numeric((L24[uY] > 0.2) | A3[uY])
    
    invisible(capture.output(bart_predict_Y5_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A4, L1 = L14, L2 = L24, L3 = L34, logE = logE4,
                   t00 = 0, t01 = 0, t02 = 0, t03 = 0, t04 = 1))$prob.test.mean))
    
    ## generate Y5
    Y5[uY] <- rbinom(sum(uY), 1, bart_predict_Y5_prob[uY])
    
    # compute estimates
    tmpdata <- data.frame(Y1, Y2, Y3, Y4, Y5)
    
    tmpdata$Y2 <- ifelse(tmpdata$Y1==1, 1, tmpdata$Y2)
    tmpdata$Y3 <- ifelse(!is.na(tmpdata$Y2) & tmpdata$Y2==1, 1, tmpdata$Y3)
    tmpdata$Y4 <- ifelse(!is.na(tmpdata$Y3) & tmpdata$Y3==1, 1, tmpdata$Y4)
    tmpdata$Y5 <- ifelse(!is.na(tmpdata$Y4) & tmpdata$Y4==1, 1, tmpdata$Y5)
    
    result_bart[i,] <- colMeans(tmpdata)
    
    if (i %% 100 == 0) {
      cat("\r", paste("Repetition", m, "MC sample", i, sep=" "))
      flush.console()
    }
    
  }
  output <- colMeans(result_bart)
  
  return(output)
}

bart_cov_bs_gform_func(1)
  
