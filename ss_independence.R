# Independence prior slab updated

ss_independence <- function(X, Y, fix = rep(0,ncol(X)), c = 1, s0 = 0, S0 = 0, c0 = 1, d0 = 1, iter = 1000, print_iter = TRUE){
  require(mvtnorm)
  
  #######################################################################
  #######################################################################
  ####################### Auxiliary functions ###########################
  #######################################################################
  #######################################################################
  
  
  ml_y <- function(a0, A0, y, X, delta){
    N <- length(X[,1])
    yc <- y-mean(y)
    sN <- (N-1)/2
    
    if(sum(delta)>0){
      X_d <- X[,delta]
      A0_d <- as.matrix(A0[delta,delta]) # Make sure it's matrix format, if it's e.g. a scalar
      a0_d <- a0[delta]
      
      inv_A0 <- solve(A0_d)
      AN_d <- solve(t(X_d) %*% X_d + inv_A0)
      aN_d <- AN_d %*% (t(X_d) %*% yc + inv_A0 %*% a0_d)
      SN_d <- 0.5 * (t(yc) %*% yc + t(a0_d) %*% inv_A0 %*% a0_d - t(aN_d) %*% solve(AN_d) %*% aN_d)
      
      log_ml_yd <- -0.5*log(N)-0.5*(N-1)*(log(2*pi))+0.5*(log(det(AN_d)))-0.5*(log(det(A0_d)))+ lgamma(sN)-sN*(log(SN_d))
    }
    else{
      SN_d <- 0.5*(t(yc)%*%yc)
      log_ml_yd <- -0.5*log(N)-0.5*(N-1)*(log(2*pi)) + lgamma(sN) - sN * (log(SN_d))
    }
    return(log_ml_yd)
  } # Calculates the marginal likelihood
  
  get_aN_AN_sn_SN <- function(a0,A0,s0,S0,y,X,delta){
    N=length(X[,1])
    yc <- y-mean(y)
    
    if(sum(delta)>0){
      X_d <- X[,delta]
      A0_d <- as.matrix(A0[delta,delta]) # Make sure it's matrix format
      a0_d <- a0[delta]
      inv_A0 <- solve(A0_d)
      
      AN_d <- solve(t(X_d) %*% X_d + inv_A0)
      aN_d <- AN_d %*% (t(X_d) %*% yc + inv_A0 %*% a0_d)
      sN <- s0 + (N-1)/2
      SN_d <- S0 + 0.5 * (t(yc) %*% yc + t(a0_d) %*% inv_A0 %*% a0_d - t(aN_d) %*% solve(AN_d) %*% aN_d)
    }
    else{
      AN_d <- matrix(1)
      aN_d <- c(1)
      sN <- s0+(N-1)/2
      SN_d <-S0+0.5*(t(yc)%*%yc)
    }
    return(list('aN' = aN_d, 'AN' = AN_d, 'sN' = sN, 'SN' = SN_d))
  }
  
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  
  # Remove NA's
  X <- X[,colSums(is.na(X)) != nrow(X)]
  
  k <- length(X[1,])
  N <- length(X[,1])
  
  # Priors for independence slab
  a0 <- rep(0,k)
  A0 <- diag(k)*c
  
  mean_y <- mean(Y)
  
  # Matrices to store chains
  delta_draws <- matrix( data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  mu_draws <- rep(NA, iter)
  alpha_draws <- matrix(data = NA, nrow = iter, ncol = k)
  omega_draws <- rep(NA,iter)
  
  # Initialize delta, omega and alpha
  delta <- rep(FALSE,k)
  alpha <- rep(0,k)
  omega <- rbeta(1,c0,d0) # Initial guess based on the prior
  
  for(j in 1:iter){
    # Step Ia: Deltas
    if(print_iter == TRUE) if (j %%100 == 0) cat("iter =", j, "\r")
    perm <- sample(1:k,k) # Draw a random permutation
    log_yd_prev <- ml_y(a0,A0,Y,X,delta)
    for(i in 1:k){
      log_yd <- rep(NA,2)
      log_delta <- rep(NA,2)
      log_post_delta <- rep(NA,2)
      
      log_yd[delta[perm[i]] + 1] <- log_yd_prev
      delta[perm[i]] <- (1 - delta[perm[i]]) == 1 # Compute the new value
      log_yd[delta[perm[i]] + 1] <- ml_y(a0,A0,Y,X,delta)

      R_j <- exp(log_yd[1] - max(log_yd))/exp(log_yd[2] - max(log_yd)) # Subtract max(log_yd) to avoid very low numbers
      post_delta <- c(1-1/(1+R_j*(1-omega)/omega),1/(1+R_j*(1-omega)/omega))
      delta[perm[i]] <- if(fix[perm[i]] == 0) sample(x=c(0,1),1, replace = TRUE, prob = post_delta) else 1
      log_yd_prev <- log_yd[delta[perm[i]] + 1] # Save value that was selected
      delta <- delta == TRUE # Convert to TRUE/FALSE
    }
    delta_draws[j,] <- delta

    # Auxiliary step: get aN, AN, sN, SN
    c <- get_aN_AN_sn_SN(a0,A0,s0,S0,Y,X,delta)
    aN <- c$aN; AN <- c$AN; sN <- c$sN; SN <- c$SN
    
    # Step Ib: sigmas
    
    sigma2 <- 1/rgamma(1,sN, SN)
    sigma_draws[j] <- sigma2
    
    # Step Ic: mu
    mu <- rnorm(1, mean = mean_y, sd = sqrt(sigma2/N))
    mu_draws[j] <- mu
    
    # Step II: omega
    omega <- rbeta(1,sum(delta) + c0, k - sum(delta) + d0)
    omega_draws[j] <- omega 
    
    # STEP III: Non-zero alphas
    alpha[delta] <- rmvnorm(1, mean = aN, sigma = sigma2 * AN )
    alpha[delta == FALSE] <- 0
    alpha_draws[j,] <- alpha
  }
  
  draws <- list('alphas' = data.frame('alphas' = alpha_draws), 'deltas' = data.frame('deltas' = delta_draws), 'mu' = mu_draws, 'sigma2' = sigma_draws, 
                'omega' = omega_draws)
  
  colnames(draws$deltas) <- colnames(draws$alphas) <- colnames(X)
  
  return(draws)
}



