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

source("datagen_dd_sim.R")

bart_cov_gform_func <- function(m){
  
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
  theta0 <- -1; theta1 <- -2; theta2 <- -0.25; theta3 <- -0.5; theta4 <- 0.4
  
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
  Y_covmodel2 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2),]$L2
  
  ## this is to fit the generative model of L2
  invisible(capture.output(bart_covmodel2 <-
    wbart(X_covmodel2, Y_covmodel2, nskip = 1000, ndpost = 10000, ntree = nt,
         nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # the generative model of L3
  X_covmodel3 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, lag_L2, lag_L3, L1, L2, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  Y_covmodel3 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3),]$L3
  
  invisible(capture.output(bart_covmodel3 <-
    wbart(X_covmodel3, Y_covmodel3, nskip = 1000, ndpost = 10000, ntree = nt,
         nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  # the generative model of L4
  X_covmodel4 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2) & 
                          !is.na(dffull$lag_L3) & !is.na(dffull$lag_L4),] %>%
    mutate(t0 = as.factor(t0)) %>%
    select(lag_A, lag_L1, lag_L2, lag_L3, lag_L4, L1, L2, L3, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  ## this is to find the Y's for the generative model of L4
  Y_covmodel4 <- dffull[!is.na(dffull$lag_L1) & !is.na(dffull$lag_L2) & 
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
    select(A, L1, L2, L3, L4, t0) %>%
    as.data.frame() %>%
    model.matrix( ~ .-1, data = .)
  
  invisible(capture.output(bart_ymodel <-
    lbart(X_ymodel, Y_ymodel, nskip = 1000, ndpost = 10000, ntree = nt,
          nkeeptest = 1, nkeeptrain = 1, nkeeptreedraws = 1)))
  
  sigma_L2 <- mean(bart_covmodel2$sigma)
  sigma_L3 <- mean(bart_covmodel3$sigma)
  sigma_L4 <- mean(bart_covmodel4$sigma)
  
  ## Monte Carlo simulation
  N <- 3000
  ids <- as.list(1:N)
  result_bart <- matrix(NA, nrow = N, ncol = K) # storing results in all five periods
  
  for (i in 1:N){
    L10 <- L20 <- L30 <- L40 <- A0 <- Y1 <-
      L11 <- L21 <- L31 <- L41 <- A1 <- Y2 <- 
      L12 <- L22 <- L32 <- L42 <- A2 <- Y3 <- as.numeric(rep(NA, N))
    
    # baseline
    ## bootstrapping the distribution of covariates
    dfbaseline <- dffull[dffull$t0 == 0,]
    idboot <- sample(1:nrow(dfbaseline), N, replace=T)
    uY <- rep(TRUE, N)
    L10 <- dfbaseline$L1[idboot]
    L20 <- dfbaseline$L2[idboot]
    L30 <- dfbaseline$L3[idboot]
    L40 <- dfbaseline$L4[idboot]
    A0 <- as.numeric((L20*L30 + cos(L40)) > 1)
    
    ## samples from BART posteriors of the generative model of Y1, probabilities from lbart
    invisible(capture.output(bart_predict_Y1_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A0, L1 = L10, L2 = L20, L3 = L30, L4 = L40, 
                   t00 = 1, t01 = 0, t02 = 0))$prob.test.mean))
    ## generate Y1
    Y1 <- rbinom(N, 1, bart_predict_Y1_prob)
    ## remaining individuals
    uY <- uY & !Y1 
    
    # time 1
    ## generate L1
    invisible(capture.output(bart_predict_t1_L1_prob <-
      predict(bart_covmodel1, newdata = 
                data.frame(lag_A = A0, lag_L1 = L10, t00 = 1, t01 = 0))$prob.test.mean))
    L11[uY] <- rbinom(sum(uY), 1, bart_predict_t1_L1_prob[uY])
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t1_L2_mean <-
      predict(bart_covmodel2, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, L1 = L11, 
                   t00 = 1, t01 = 0))[1,]))
    L21[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t1_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t1_L3_mean <-
      predict(bart_covmodel3, newdata =
        data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, lag_L3 = L30, 
                   L1 = L11, L2 = L21, t00 = 1, t01 = 0))[1,]))
    L31[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t1_L3_mean[uY], sd = sigma_L3)
    
    ## generate L4
    invisible(capture.output(bart_covmodel4_predict_t1_L4_mean <-
      predict(bart_covmodel4, newdata =
                data.frame(lag_A = A0, lag_L1 = L10, lag_L2 = L20, lag_L3 = L30, 
                           lag_L4 = L40, L1 = L11, L2 = L21, L3 = L31, t00 = 1, t01 = 0))[1,]))
    L41[uY] <- rnorm(sum(uY), mean = bart_covmodel4_predict_t1_L4_mean[uY], sd = sigma_L4)
    
    A1[uY] <- as.numeric(((L21[uY]*L31[uY] + cos(L41[uY])) > 1) | A0[uY])
    
    invisible(capture.output(bart_predict_Y2_prob <-
      predict(bart_ymodel, newdata =
        data.frame(A = A1, L1 = L11, L2 = L21, L3 = L31, L4 = L41, 
                   t00 = 0, t01 = 1, t02 = 0))$prob.test.mean))
    ## generate Y2
    Y2[uY] <- rbinom(sum(uY), 1, bart_predict_Y2_prob[uY])
    ## remaining individuals
    uY <- uY & !Y2 
    
    # time 2
    ## generate L1
    invisible(capture.output(bart_predict_t2_L1_prob <-
      predict(bart_covmodel1, newdata = 
                data.frame(lag_A = A1, lag_L1 = L11, t00 = 0, t01 = 1))$prob.test.mean))
    L12[uY] <- rbinom(sum(uY), 1, bart_predict_t2_L1_prob[uY])
    
    ## generate L2
    invisible(capture.output(bart_covmodel2_predict_t2_L2_mean <-
      predict(bart_covmodel2, newdata =
                data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, L1 = L12, t00 = 0, t01 = 1))[1,]))
    L22[uY] <- rnorm(sum(uY), mean = bart_covmodel2_predict_t2_L2_mean[uY], sd = sigma_L2)
    
    ## generate L3
    invisible(capture.output(bart_covmodel3_predict_t2_L3_mean <-
      predict(bart_covmodel3, newdata =
                data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, lag_L3 = L31, 
                           L1 = L12, L2 = L22, t00 = 0, t01 = 1))[1,]))
    L32[uY] <- rnorm(sum(uY), mean = bart_covmodel3_predict_t2_L3_mean[uY], sd = sigma_L3)
    
    ## generate L4
    invisible(capture.output(bart_covmodel4_predict_t2_L4_mean <-
      predict(bart_covmodel4, newdata =
                data.frame(lag_A = A1, lag_L1 = L11, lag_L2 = L21, lag_L3 = L31, 
                           lag_L4 = L41, L1 = L12, L2 = L22, L3 = L32, t00 = 0, t01 = 1))[1,]))
    L42[uY] <- rnorm(sum(uY), mean = bart_covmodel4_predict_t2_L4_mean[uY], sd = sigma_L4)
    
    A2[uY] <- as.numeric(((L22[uY]*L32[uY] + cos(L42[uY])) > 1) | A1[uY])
    
    invisible(capture.output(bart_predict_Y3_prob <-
      predict(bart_ymodel, newdata =
                data.frame(A = A2, L1 = L12, L2 = L22, L3 = L32, L4 = L42, 
                           t00 = 0, t01 = 0, t02 = 1))$prob.test.mean))
    ## generate Y3
    Y3[uY] <- rbinom(sum(uY), 1, bart_predict_Y3_prob[uY])
    
    # compute estimates
    tmpdata <- data.frame(Y1, Y2, Y3)
    
    tmpdata$Y2 <- ifelse(tmpdata$Y1 == 1, 1, tmpdata$Y2)
    tmpdata$Y3 <- ifelse(!is.na(tmpdata$Y2) & tmpdata$Y2==1, 1, tmpdata$Y3)
    
    result_bart[i,] <- colMeans(tmpdata)
    
    if (i %% 100 == 0) {
      cat("\r", paste("Repetition", m, "MC sample", i, sep=" "))
      flush.console()
    }
    
  }
  output <- colMeans(result_bart)
  
  return(output)
}

bart_cov_gform_func(1)
