
# Function to create survival data in long form

# Regimen: deterministic & dynamic (monotone treatment)
# Functions: nonlinear

# Two levels of separation can be selected by altering model coefficients
# It simulates ONE INDIVIDUAL's data at t = {0, 1, 2, ..., K-1} (outcomes realized at t = {1, 2, ..., K})

# Y is binary outcome
# Right censoring is present (C)
# A is binary treatment

# No competing risk
# 2 covariates (L1, L2)

# Inputs are i: id, df: data frame specified above, K: max number of time points
# Alphas, betas, thetas, etas, and sigma are user-provided parameters for data generating functions

# treatment assignment strategy is assumed as monotone

# have assumed a pooled model over time, so each set of coefficient applies to the conditional 
# outcome (exposure) model at each time slice

# we have commented out the more complex model where the outcome depends on whole treatment history
# consider specifying nonlinear relationship


datagen_dd_nonlinear <- function(i, K, beta0, beta1, beta2, beta3,
                                 theta0, theta1, theta2, theta3, theta4,
                                 eta0, eta1, eta2, eta3, eta4, sigma) {
  
  id <- as.numeric(i)
  # Data at time 0
  # two covariates; one continuous and one binary
  L1 <- rbinom(1, 1, 0.5) # binary
  L2 <- rnorm(1, mean = beta0, sd = sigma) # continuous
  L3 <- rnorm(1, mean = beta0, sd = sigma) # continuous
  
  # linear functional structure for everything
  A <- as.numeric(L2 > 0.2) # treat if L2 is greater than 0.1
  
  ## censoring, depending on the treatment assignment and covariates
  C <- rbinom(1, 1, plogis(eta0 + eta1*A + eta2*L1*cos(eta3*L2) + eta3*L2*L3))
  
  # probability of not being censored
  omega <- 1 - plogis(eta0 + eta1*A + eta2*L1*cos(eta3*L2) + eta3*L2*L3) 
  
  # probability of being treated
  zeta <- as.numeric(L2 > 0.2)
  
  # propensity score
  ps <- omega*zeta
  
  # linear function for the outcome
  Y <- rbinom(1, 1, plogis(theta0 + theta1*A + theta2*L1 + theta3*L2*L3 - (theta3)*(L1*L2^2)))  
  temp <- data.frame(id = id, t0 = 0, L1, L2, L3, A, C, Y, omega, zeta, ps)
  
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
        plogis(beta0 + beta1*temp$A[j-1]) + beta2*temp$L1[j-1]) | temp$L1[j-1])
      L2star <- rnorm(1, mean = beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1] + 
                        beta3*temp$L2[j-1]*temp$L3[j-1] + beta3*sin(temp$L2[j-1]), sd = sigma)
      L3star <- rnorm(1, mean = beta0 + beta1*temp$A[j-1] + beta2*temp$L1[j-1]*temp$L3[j-1] + 
                        beta3*temp$L2[j-1] + beta3*sin(temp$L3[j-1]), sd = sigma)
      temp_L1 <- c(temp$L1, L1star)
      temp_L2 <- c(temp$L2, L2star)
      temp_L3 <- c(temp$L3, L3star)
      
      # generate new assignments
      Astar <- as.numeric(L2star > 0.2 | temp$A[j-1])
      temp_A <- c(temp$A, Astar)
      
      # generate new censoring status
      Cstar <- 
        rbinom(1, 1, plogis(eta0 + eta1*Astar + eta2*L1star*cos(eta3*L2star) + 
                              eta3*L2star*L3star + eta4*t0))
      
      # generate new probability of not being censored
      omegastar <- 1 - plogis(eta0 + eta1*Astar + eta2*L1star*cos(eta3*L2star) +
                                eta3*L2star*L3star + eta4*t0)
      
      # generate new probability of being treated
      zetastar <- as.numeric(L2star > 0.2 | temp$A[j-1])
      
      # generate new propensity score
      psstar <- omegastar*zetastar
      
      if (Cstar == 1) {
        # if censored, Y is set to NA
        Ystar <- NA
        temp <- rbind(temp, c(id, t0, L1star, L2star, L3star, Astar, Cstar, Ystar, 
                              omegastar, zetastar, psstar))
        break
      } else {
        # if not censored, generate outcome
        Ystar <- rbinom(1, 1, plogis(theta0 + theta1*Astar + theta2*L1star + 
                                       theta3*L2star*L3star - (theta3)*(L1star*L2star^2) + 
                                       theta4*t0))
      }
      
      temp <- rbind(temp, c(id, t0, L1star, L2star, L3star, Astar, Cstar, Ystar, 
                            omegastar, zetastar, psstar))
      if (Ystar == 1) { break }
    }
  }
  return (temp)
}

# # testing the above function
# # set up and parameters
# n <- 1000 # number of participants
# K <- 5 # number of observation periods          
# 
# # exposure model
# eta0 <- -1.5; eta1 <- -1; eta2 <- 0.75; eta3 <- -0.5; eta4 <- 0
# 
# # covariate model
# beta0 <- 0; beta1 <- -2; beta2 <- 0.2; beta3 <- 1
# # continuous covariate sd
# sigma <- 0.1 

# # outcome model
# theta0 <- -2; theta1 <- -3; theta2 <- 1; theta3 <- -6; theta4 <- 0
# 
# # generate data for all 1000 patients
# df <- lapply(
#   as.list(1:n),
#   FUN = function(ind) {
#     datagen_dd_nonlinear(
#       ind,
#       K = K,
#       beta0 = beta0,
#       beta1 = beta1,
#       beta2 = beta2,
#       beta3 = beta3,
#       theta0 = theta0,
#       theta1 = theta1,
#       theta2 = theta2,
#       theta3 = theta3,
#       theta4 = theta4,
#       eta0 = eta0,
#       eta1 = eta1,
#       eta2 = eta2,
#       eta3 = eta3,
#       eta4 = eta4,
#       sigma = sigma
#     )
#   }
# )
# 
# sum(unlist(lapply(df, function(x) 1*(sum(x$C) > 0))))


