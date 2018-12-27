
ss_mixture <- function(X, Y, groups = 2, fix = rep(0, ncol(X)), c = 1, iter = 1000, print_iter = TRUE){
  require(mvtnorm)
  require(MCMCpack)
  
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
  
  sample_omega <- function(delta, omega, p, s, m, groups = groups){
    options( warn = -1 )
    # Hyperpriors
    delta_prior <- 1; rho <- 0.1; kappa <- 0.00001; sigma2 <-  0.1
    
    # Posterior distribution of s and m function
    log_post_sm2 <- function(s_j, m_j, s_new, m_new, n_j, delta_prior, kappa, rho, omega_j){
      j <- log(exp(m_j) / ((1+exp(m_j))^2))
      j_new <- log(exp(m_new) / ((1+exp(m_new))^2))
      
      m_j <- exp(m_j) / (1+exp(m_j))
      m_new <- exp(m_new) / (1+exp(m_new))
      output <- log((1- exp(-delta_prior*((s_new-2)^2+(m_new-0.5)^2)))) - log((1- exp(-delta_prior*((s_j-2)^2+(m_j-0.5)^2))))
      output <- output + log(exp(-rho/(s_new^2*m_new*(1-m_new)) - kappa * s_new^2/2)) - log(exp(-rho/(s_j^2*m_j*(1-m_j)) - kappa * s_j^2/2))
      output <- output + n_j * log((gamma(s_new)/(gamma(s_new*m_new) * gamma(s_new*(1-m_new))))) - (n_j * log((gamma(s_j)/(gamma(s_j*m_j) * gamma(s_j*(1-m_j))))))
      output <- output + (s_new * m_new) * sum(log(omega_j)) + (s_new*(1-m_new)) * sum(log(1-omega_j)) - ((s_j * m_j) * sum(log(omega_j)) + (s_j*(1-m_j)) * sum(log(1-omega_j)))
      output <- output + j - j_new
      return(output)
    }
    
    dlnormal <- function(x, sigma2, mu){
      output <- (1/(x * sqrt(2*sigma2*pi))) * exp(-((log(x) - mu)^2)/(2*sigma2))
      return(output)
    }
    
    k <- length(omega)
    
    # Storage matrices
    prob <- matrix(NA, nrow = groups, ncol = k) 
    z <- matrix(NA, nrow = groups, ncol = k) 
    
    # Prior on groups
    gamma <- rep(1, groups)
    
    # Sample z
    for(i in 1:groups){
      #prob[i,] <- p[i] * new_beta(omega, s[i], m[i]) 
      prob[i,] <- p[i] * dbeta(omega, shape1 = s[i] * m[i] , shape2 = s[i] * (1-m[i]))
    }
    
    for(j in 1:k){
      prob[,j] <- if(sum(prob[,j]) > 0 & sum(is.na(prob[,j])) == 0 & sum(prob[,j]) != Inf) prob[,j] / sum(prob[,j]) else rep(1/groups,groups)
      z[,j] <- as.integer(rmultinom(1,1,prob[,j]))
    }
    
    n <- rep(NA,groups)
    for(i in 1:groups) n[i] <- sum(z[i,])
    
    Z <- rep(NA,k) # Create categorical group category, e.g. Z[j] = 2 when omega_j belongs to group 2
    for(j in 1:k) Z[j] <- which(z[,j] == 1)
    
    # Sample p
    p <- rdirichlet(1, gamma + n)
    
    # Sampling s and m
    m <- m/(1-m) # Transformation
    
    ## Steps d: Sampling the new hyperpriors
    s_new <- rep(NA,groups)
    m_new <- rep(NA, groups)
    
    # Sample from proposal distribution
    for(i in 1:groups){
      s_new[i] <- exp(log(s[i]) + rnorm(1,sd = sqrt(sigma2)))
      m_new[i] <- exp(log(m[i]) + rnorm(1,sd = sqrt(sigma2)))
      
      # Metropolis-Hastings Step
      log_r <- log_post_sm2(s[i], m[i], s_new[i], m_new[i], n[i], delta_prior, kappa, rho, omega[Z == i])
      r <- exp(log_r)
      r <- r * dlnormal(s[i], sigma2, s_new[i]) / dlnormal(s_new[i], sigma2, s[i])
      u <- runif(1)
      
      
      if(is.na(r) == TRUE) r <- 1
      if(r < u){
        s[i] <- s_new[i]
        m[i] <- m_new[i]
      }
    } 
    m <- m/(1+m) # Return to original formulation
    m_exp <- m
    
    i_delta <- as.integer(delta!=0)
    for(j in 1:k){
      omega[j] <- rbeta(1, shape1 = s[Z[j]] * m[Z[j]]  + i_delta[j], shape2 = s[Z[j]] * (1-m[Z[j]]) + 1 - i_delta[j])
    }
    return(list("omega" = omega, "p" = p, "z" = Z, "s" = s, "m" = m))
  } # Samples omegas via a beta mixture with MH Step
  
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  
  k <- length(X[1,])
  N <- length(X[,1])
  
  # Priors for independence slab
  a0 <- rep(0,k)
  A0 <- diag(k)*c
  
  # Initialization for mixture components
  s <- rep(4, groups)
  m <- c(0.875,0.125,rep(0.5,groups-2))
  p <- rep(1/groups,groups)# Prior inclusion probability initialization
  
  # Matrices to store chains
  delta_draws <- matrix( data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  mu_draws <- rep(NA, iter)
  beta_draws <- matrix(data = NA, nrow = iter, ncol = k)
  omega_draws <- matrix(data = NA, nrow = iter, ncol = k)
  z_draws <- matrix(data = NA, nrow = iter, ncol = k) # Groups draws
  
  s_draws <- matrix(NA, nrow = iter, ncol = groups)
  m_draws <- matrix(NA, nrow = iter, ncol = groups)
  
  # Initialize delta, omega and beta (and gamma and p for sampling omega)
  delta <- rep(FALSE,k) == fix # Initialize to include deltas we want in the model, rest out.
  beta <- rep(0,k)
  
  gamma <- rep(1,groups) 
  omega <- runif(k); omega[fix == 1] <- 1 # Initialize omega
  
  for(j in 1:iter){ 
    if(print_iter == TRUE) if (j %%100 == 0) cat("iter =", j, "\r")
    
    # Step Ia: Deltas
    perm <- sample(1:k,k) # Draw a random permutation
    log_yd_prev <-  ml_y(a0,A0,Y,X,delta)
    for(i in 1:k){
      log_yd <- rep(NA,2)
      
      log_yd[delta[perm[i]] + 1] <- log_yd_prev
      delta[perm[i]] <- (1 - delta[perm[i]]) == 1 # Compute the new value
      log_yd[delta[perm[i]] + 1] <-  ml_y(a0,A0,Y,X,delta)
      
      R_j <- exp(log_yd[1] - max(log_yd))/exp(log_yd[2] - max(log_yd)) # Subtract max(log_yd) to avoid very low numbers
      post_delta <- c(1-1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])),1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])))
      delta[perm[i]] <- if(fix[perm[i]] == 0) sample(x=c(0,1),1, replace = TRUE, prob = post_delta) else 1
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
    
    # Step II: omega
    df <- sample_omega(delta,omega,p,s,m,groups = groups)
    omega <- df$omega; p <- df$p; s <- df$s; m <- df$m; z <- df$z
    omega_draws[j,] <- omega 
    
    z_draws[j,] <- z
    s_draws[j,] <- s
    m_draws[j,] <- m
    
    # STEP III: Non-zero betas
    beta[delta] <- rmvnorm(1, mean = aN, sigma = sigma2 * AN )
    beta[delta == FALSE] <- 0
    beta_draws[j,] <- beta
  }
  draws <- list('betas' = data.frame('betas' = beta_draws), 'deltas' = data.frame('deltas' = delta_draws), 
                'mu' = mu_draws, 'sigma2' = sigma_draws, 
                'omegas' = data.frame('omegas' = omega_draws), 'z' = data.frame('z' = z_draws),
                's' = data.frame('s' = s_draws), 'm' = data.frame('m' = m_draws))
  
  
  colnames(draws$deltas) <- colnames(draws$betas) <- colnames(X)
  
  return(draws)
  
}
