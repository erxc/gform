
datagen_dd_nonlinear <- function(i, K, beta0, beta1, beta2, beta3, beta4,
                                 theta0, theta1, theta2, theta3, theta4,
                                 eta0, eta1, eta2, eta3, eta4, sigma) {
  
  id <- as.numeric(i)
  # Data at time 0
  # four covariates; three continuous and one binary
  L1 <- rbinom(1, 1, 0.5) # binary
  L2 <- rnorm(1, mean = beta0, sd = sigma) # continuous
  L3 <- rnorm(1, mean = beta0, sd = sigma) # continuous
  L4 <- rnorm(1, mean = beta0, sd = sigma) # continuous
  
  # linear functional structure for everything
  A <- as.numeric((L2*L3 + cos(L4)) > 1) # treat if L2*L3+cos(L4) is greater than 0.95
  
  ## censoring, depending on the treatment assignment and covariates
  C <- rbinom(1, 1, plogis(eta0 + eta1*A + eta2*L1*cos(eta3*L2) + eta3*L2*L3 + eta4*sin(L4)))
  
  # probability of not being censored
  omega <- 1 - plogis(eta0 + eta1*A + eta2*L1*cos(eta3*L2) + eta3*L2*L3 + eta4*sin(L4)) 
  
  # probability of being treated
  zeta <- as.numeric((L2*L3 + cos(L4)) > 1)
  
  # propensity score
  ps <- omega*zeta
  
  # linear function for the outcome
  Y <- rbinom(1, 1, plogis(theta0 + theta1*A + (theta2*log(abs(L1)) + theta3*L2*L3)^(-1) - 
                             (theta3)*(L1*L2^2) + theta4*abs(L4)))  
  temp <- data.frame(id = id, t0 = 0, L1, L2, L3, L4, A, C, Y, omega, zeta, ps)
  
  if (temp$C == 1) {
    # if censored, then Y is set to NA
    temp$Y <- NA
  } else if (temp$Y != 1) {
    # if not censored and the event has not occurred, proceed
    for (j in 2:K) {
      # time effect since t=0
      t0 <- j-1
      
      # generate new covariates
      L1star <- as.numeric(rbinom(1, 1, 
                   plogis(beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1])) | temp$L1[j-1])
      L2star <- rnorm(1, mean = beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1] + 
                        beta3*temp$L2[j-1]*temp$L3[j-1] + beta3*sin(temp$L2[j-1]), sd = sigma)
      L3star <- rnorm(1, mean = beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1]*temp$L3[j-1] + 
                        beta3*temp$L2[j-1] + beta3*sin(temp$L3[j-1]), sd = sigma)
      L4star <- rnorm(1, mean = beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1]*temp$L3[j-1] + 
                        beta3*temp$L2[j-1] + beta3*sin(temp$L3[j-1]) + beta4*cos(temp$L4[j-1]), 
                      sd = sigma)
      temp_L1 <- c(temp$L1, L1star)
      temp_L2 <- c(temp$L2, L2star)
      temp_L3 <- c(temp$L3, L3star)
      temp_L4 <- c(temp$L4, L4star)
      
      # generate new assignments
      Astar <- as.numeric((L2star*L3star + cos(L4star)) > 1 | temp$A[j-1])
      temp_A <- c(temp$A, Astar)
      
      # generate new censoring status
      Cstar <- 
        rbinom(1, 1, plogis(eta0 + eta1*Astar + eta2*L1star*cos(eta3*L2star) + 
                              eta3*L2star*L3star + eta4*sin(L4star)))
      
      # generate new probability of not being censored
      omegastar <- 1 - plogis(eta0 + eta1*Astar + eta2*L1star*cos(eta3*L2star) +
                                eta3*L2star*L3star + eta4*sin(L4star))
      
      # generate new probability of being treated
      zetastar <- as.numeric((L2star*L3star + cos(L4star)) > 1 | temp$A[j-1])
      
      # generate new propensity score
      psstar <- omegastar*zetastar
      
      if (Cstar == 1) {
        # if censored, Y is set to NA
        Ystar <- NA
        temp <- rbind(temp, c(id, t0, L1star, L2star, L3star, L4star, Astar, Cstar, Ystar, 
                              omegastar, zetastar, psstar))
        break
      } else {
        # if not censored, generate outcome
        Ystar <- rbinom(1, 1, 
                        plogis(theta0 + theta1*Astar + (theta2*log(abs(L1star)) 
                                                        + theta3*L2star*L3star)^(-1) - 
                                 (theta3)*(L1star*L2star^2) + theta4*abs(L4star)))
      }
      
      temp <- rbind(temp, c(id, t0, L1star, L2star, L3star, L4star, Astar, Cstar, Ystar, 
                            omegastar, zetastar, psstar))
      if (Ystar == 1) { break }
    }
  }
  return (temp)
}
