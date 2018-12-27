# Example code for implementing the spike-and-slab models
rm(list=ls())

# Set working directory
setwd("Insert_own_directory")

##################################################################################################
##################################################################################################
# This script holds example code for estimating the three HDS models discussed in the paper
# "High-dimensional Econometrics - A Bayesian Perspective", namely: 
# 1) The Standard Linear Regression Model
# 2) Treatment Effects Model
# 3) Finite Mixture Model of Linear Regressions
##################################################################################################
##################################################################################################

##########################################################################################
################################### Preliminaries ######################################## 
##########################################################################################

# Load packages
Packages <- c('mvtnorm', 'coda', 'MCMCpack', 'glmnet')
invisible(lapply(Packages, require, character.only = TRUE))

# Functions to generate models
gen_std_linmod_corr <- function(n = 1000, mu = 1, beta1 = c(1,1), sigma2 = 1, seed = 1) {
  set.seed(seed)
  
  k <- length(beta1)
  # Generate correlation matrix
  corr_x <- matrix(data=0, nrow = k, ncol = k)
  for(i in 1:k){
    for(j in 1:k){
      corr_x[i,j] <- 0.5**abs(i-j)
    }
  }
  
  X <-rmvnorm(n, mean = rep(0,k), sigma = corr_x)
  X <- scale(X)
  e <- rnorm(n, sd = sigma2^0.5)
  ones <- rep(1, n)
  
  Y <- ones*mu + X %*% beta1 + e
  
  return(list(Y=Y, X = X))
}

# Generate linear model for treatment effects
gen_std_linmod_treatment <- function(mu = 1, n = 100, k = 100, alpha_0 = 1, sigma_zeta = 1, sigma_v = 1, seed = 1){
  require('mvtnorm')
  set.seed(seed)
  beta_0 <- c(1/(1:5), rep(0,5), 1/(1:5), rep(0,k-15))
  eta_0 <- c(1/(1:10),rep(0,k-10))
  
  # Generate correlation matrix
  corr_x <- matrix(data=0, nrow = k, ncol = k)
  for(i in 1:k){
    for(j in 1:k){
      corr_x[i,j] <- 0.5**abs(i-j)
    }
  }
  
  x <-rmvnorm(n, mean = rep(0,k), sigma = corr_x)
  zeta <- rnorm(n, mean = 0, sd = sigma_zeta)
  v <- rnorm(n, mean = 0, sd = sigma_v)
  
  d <- x%*%eta_0 + v
  x <- scale(x, scale = FALSE) # Center variables
  d <- scale(d, scale = FALSE) # Center variables
  
  y <- mu + alpha_0*d + x%*%beta_0 + zeta
  
  df <- list('d' = d,'X' = data.frame('X' = x), 'Y' = y)
  return(df)
}

# Generate 2-component mixture model

gen_mixture_model <- function(nobs = 100, pi_corr = 0.5, rho = 0.5, sigma2 = 1, beta1 = c(1,0,0,3,0), 
                                 beta2 = c(-1,2,0,0,3), seed = 1){
  set.seed(seed)
  require(mvtnorm)
  k <- length(beta1)
  
  corr_x <- matrix(data=0, nrow = k, ncol = k)
  for(i in 1:k){
    for(j in 1:k){
      corr_x[i,j] <- pi_corr**abs(i-j)
    }
  }
  
  groups <- rbinom(nobs, 1, rho)
  X <- rmvnorm(nobs, mean = rep(0,k), sigma = corr_x)
  X <- scale(X, scale = FALSE) # Center X
  XB1 <- X %*% beta1
  XB2 <- X %*% beta2
  y <- groups * rnorm(nobs, mean = XB1, sd = sqrt(sigma2)) + (1-groups) * rnorm(nobs, mean = XB2, sd = sqrt(sigma2))
  
  return(list(y=y, X = X))
  
}

##########################################################################################
#### 1) Example with a Dirac spike, independence slab and common beta prior on omega ##### 
##########################################################################################

nobs <- 40 # Number of observations
beta1 <- c(1/1:5, rep(0,nobs-5)) # Sparse, high-dimensional coefficient vector
kvar <- length(beta1)
c <- sqrt(nobs) # Tuning parameter for the independence slab: High c -> High sparsity
seed <- 2018
sigma2 <- 0.1 # Strong signal
iter_MCMC <- 3000 # Number of MCMC iterations
burn_in <- 1000 # Burn-in

# Generate data
data <- gen_std_linmod_corr(mu = 1, beta1 = beta1, n = nobs, sigma2 = sigma2, seed = seed )
X <- data$X; Y <- data$Y

source('ss_common.R')
draws_common <- ss_common(X,Y, iter = iter_MCMC, c = c, print_iter = TRUE)

gamma <- draws_common$deltas # Inclusion indicator variables
post_gamma <- colMeans(gamma) # Posterior means of gamms

# Correctly classified coefficients using the Median Probability Criterion
cat('Correctly classified', sum(post_gamma[0:5] > 0.5), 'out of 5 non-zero coefficients', "\n", 
    'Correctly classified', sum(post_gamma[5:40] < 0.5), 'out of 35 zero coefficients')

# Posterior distributions
betas <- draws_common$betas[burn_in:iter_MCMC,]
betas <- as.mcmc(betas)
plot(betas) # Posterior distributions

##########################################################################################
########################### 2) Example with a Treatment Prior ############################
##########################################################################################

sigma_zeta <- 1 # Variance level corresponding to a signal-to-noise ratio of 7
nobs <- 100
kvar <- 100

source('ss_treatment.R')
df <- gen_std_linmod_treatment(n = nobs, k = kvar, alpha_0 = 1, sigma_zeta = sigma_zeta, sigma_v = 1, seed = seed)
d <- df$d; x <- df$X; Y <- df$Y # d is the treatment variable, x the confounders, y the outcome
x <- as.matrix(x)
X <- as.matrix(cbind(d,x))
c <- sqrt(nobs)

draws_treat <- ss_treatment(d,x,Y, iter = iter_MCMC, c = c,  print_iter = TRUE)

# Posterior distribution of treatment parameter
alpha_zero_hat <- draws_treat$betas[burn_in:iter_MCMC,]
alpha_zero_hat <- as.mcmc(alpha_zero_hat[,1])
plot(alpha_zero_hat)

# Posterior mean estimator
post_mean <- as.numeric(summary(alpha_zero_hat)$'statistics'[1])
cat('The posterior mean estimator is', post_mean)

###########################################################################################
### 3) Example with a finite mixutre model of linear regressions with common beta prior ###
###########################################################################################

rho <- 0.5 # Mixing proportions
pi_corr <- 0.5 # Correlation
beta1 <- c(0,0,1,1/2,1/3,1,rep(0,9))
beta1_true <- beta1 != 0; beta1_false <- (1-beta1_true) == TRUE
beta2 <- c(0,0,0,0,1,-4,1/2,1/3,rep(0,7))
beta2_true <- beta2 != 0; beta2_false <- (1-beta2_true) == TRUE
groups <- 2

df <- gen_mixture_model(nobs = nobs, rho = rho, beta1 = beta1, beta2 = beta2, pi_corr = pi_corr, seed = seed)
X <- df$X; y <- df$y  

source('ss_finite_mixture_reg_individual.R')
draws_mixture <- ss_sparse_mixture(y,X, groups = 2, lambda = 0, iter = iter_MCMC, print_iter = TRUE)

betas <- draws_mixture$betas # Posterior distributions of regression coefficients
betas1 <- betas[,,2][burn_in:iter_MCMC,]
betas2 <- betas[,,1][burn_in:iter_MCMC,] # Clearly, this is component two by looking at the 6th column

deltas <- draws_mixture$deltas # Inclusion indicators
delta_mean1 <- colMeans(deltas[,,2][burn_in:iter_MCMC,])
delta_mean2 <- colMeans(deltas[,,1][burn_in:iter_MCMC,])

wrong_zeros1 <- sum(delta_mean1[beta1_true] < 0.5) 
wrong_zeros2 <- sum(delta_mean2[beta2_true] < 0.5)

right_zeros1 <- sum(delta_mean1[beta1_false] < 0.5)
right_zeros2 <- sum(delta_mean2[beta2_false] < 0.5)

# Number of correctly classified zeros (out of 11) and incorrectly classified zeros (out of 4)
cat('Correct, Component 1:', sum(right_zeros1), 'Incorrect, Component 1:', sum(wrong_zeros1),  "\n",
    'Correct, Component 2:', sum(right_zeros2), 'Incorrect, Component 2', sum(wrong_zeros2))

# Posterior mean estimates
## Component 1
cat('Posterior mean estimates of non-zero coefficients in Component 1 is', "\n", 
    round(colMeans(betas1[,3:6]), digits = 2))
cat('True values of Component 1 is', "\n", beta1[3:6])

## Component 2
cat('Posterior mean estimates of non-zero coefficients in Component 2 is', "\n", colMeans(betas2[,5:8]))
cat('True values of non-zero coefficients in Component 2 is', "\n", beta2[5:8])
