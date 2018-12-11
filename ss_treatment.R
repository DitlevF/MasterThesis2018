
ss_treatment <- function(d,X, Y, kappa = 10, phi = 0.3, fix_regr = rep(0, ncol(X)-1), c = 1, iter = 1000, print_iter = TRUE){
  require(mvtnorm)
  require(MCMCpack)
  require(glmnet)
  
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
  
  get_aN_AN_sn_SN <- function(a0,A0,y,X,delta){
    N=length(X[,1])
    yc <- y-mean(y)
    sN <- (N-1)/2
    
    if(sum(delta)>0){
      X_d <- X[,delta]
      A0_d <- as.matrix(A0[delta,delta]) # Make sure it's matrix format
      a0_d <- a0[delta]
      inv_A0 <- solve(A0_d)
      
      AN_d <- solve(t(X_d) %*% X_d + inv_A0)
      aN_d <- AN_d %*% (t(X_d) %*% yc + inv_A0 %*% a0_d)
      SN_d <- 0.5 * (t(yc) %*% yc + t(a0_d) %*% inv_A0 %*% a0_d - t(aN_d) %*% solve(AN_d) %*% aN_d)
    }
    else{
      AN_d <- matrix(1)
      aN_d <- c(1)
      SN_d <-0.5*(t(yc)%*%yc)
    }
    return(list('aN' = aN_d, 'AN' = AN_d, 'sN' = sN, 'SN' = SN_d))
  } # Calculates posterior parameters for regressors
  
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  
  X <- as.matrix(cbind(d,X))
  fix <- c(1, fix_regr) # Fix d as well as any wished-for regressors
  k <- length(X[1,])
  N <- length(X[,1])
  
  # Priors for independence slab
  a0 <- rep(0,k)
  A0 <- diag(k)*c
  
  # Matrices to store chains
  delta_draws <- matrix( data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  mu_draws <- rep(NA, iter)
  alpha_draws <- matrix(data = NA, nrow = iter, ncol = k)
  omega_draws <- matrix(data = NA, nrow = iter, ncol = k)
  z_draws <- matrix(data = NA, nrow = iter, ncol = k) # Groups draws
  gamma_draws <- matrix(data = NA, nrow = iter, ncol = k)
  
  xi <- rep(1,groups) 
  p <- rdirichlet(1,xi) # Prior inclusion probability initialization
  
  # Initialize delta, omega and alpha (and gamma and p for sampling omega)
  delta <- rep(TRUE,k) == fix # Initialize to include deltas we want in the model, rest out.
  alpha <- rep(0,k)
  
  # MCMC 1: Obtain p(delta given X), first stage
  ## Approximate this by the lasso for faster computing
  cv.lasso <- cv.glmnet(X[,-1],d)
  coef_lasso <- coef(cv.lasso, s = "lambda.1se")
  indices <- which(coef_lasso !=0)[-1] 
  
  omega <- rep(0,k)
  omega[indices] <- kappa/(1+kappa)
  omega[-indices] <- phi
  
  for(j in 1:iter){ 
    if(print_iter == TRUE) if (j %%100 == 0) cat("iter =", j, "\r")
    
    # Step Ia: Deltas
    perm <- sample(1:k,k) # Draw a random permutation
    log_yd_prev <-  ml_y(a0,A0,Y,X,delta)
    for(i in 1:k){
      if((fix[perm[i]] == 1)) next # Skip the fixed regressors
      log_yd <- rep(NA,2)
      
      log_yd[delta[perm[i]] + 1] <- log_yd_prev
      delta[perm[i]] <- (1 - delta[perm[i]]) == 1 # Compute the new value
      log_yd[delta[perm[i]] + 1] <-  ml_y(a0,A0,Y,X,delta)
      
      R_j <- exp(log_yd[1] - max(log_yd))/exp(log_yd[2] - max(log_yd)) # Subtract max(log_yd) to avoid very low numbers
      post_delta <- c(1-1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])),1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])))
      delta[perm[i]] <- sample(x=c(0,1),1, replace = TRUE, prob = post_delta)
      log_yd_prev <- log_yd[delta[perm[i]] + 1] # Save value that was selected
      delta <- delta == TRUE # Convert to TRUE/FALSE
    }
    delta_draws[j,] <- delta
    
    # Auxiliary step: get aN, AN, sN, SN
    c <- get_aN_AN_sn_SN(a0,A0,Y,X,delta)
    aN <- c$aN; AN <- c$AN; sN <- c$sN; SN <- c$SN
    
    # Step Ib: sigmas
    sigma2 <- 1/rgamma(1,sN, SN)
    sigma_draws[j] <- sigma2
    
    # Step Ic: mu
    mu <- rnorm(1, mean = mean(Y), sd = sqrt(sigma2/N))
    mu_draws[j] <- mu
    
    # STEP III: Non-zero alphas
    alpha[delta] <- rmvnorm(1, mean = aN, sigma = sigma2 * AN )
    alpha[delta == FALSE] <- 0
    alpha_draws[j,] <- alpha
  }
  
  draws <- list('alphas' = data.frame('alphas' = alpha_draws), 'deltas' = data.frame('deltas' = delta_draws), 
                'mu' = mu_draws, 'sigma2' = sigma_draws, 
                'omegas' = data.frame('omegas' = omega_draws), 'z' = data.frame('z' = z_draws))
  
  colnames(draws$deltas) <- colnames(draws$alphas) <- colnames(X)
  
  return(draws)
  
}
