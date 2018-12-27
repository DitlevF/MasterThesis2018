ss_sparse_mixture <- function(y,X, groups = 2, lambda = 0, iter = 1000, print_iter = TRUE){
  require(mvtnorm)
  require(MCMCpack)
  
  # AUXILIARY
  sample_Z_rho <- function(y, X, beta, sigma2, rho, alpha){
    
    N <- length(y)
    XB <- X %*% beta
    groups <- dim(beta)[2]
    
    # Sample z
    prob <- matrix(NA, nrow = groups, ncol = N) 
    z <- matrix(NA, nrow = groups, ncol = N) 
    
    for(i in 1:groups){
      prob[i,] <- rho[i] * dnorm(y, mean = XB[,i] , sd = sqrt(sigma2[i])) 
    }
    for(j in 1:N){
      prob[,j] <- prob[,j] / sum(prob[,j])
      z[,j] <- as.integer(rmultinom(1,1,prob[,j]))
    }
    n <- rep(NA,groups)
    for(i in 1:groups) n[i] <- sum(z[i,])
    
    Z <- rep(NA,N) # Create categorical group category, e.g. Z[j] = 2 when omega_j belongs to group 2
    for(j in 1:N) Z[j] <- which(z[,j] == 1)
    
    # Sample rho
    rho <- rdirichlet(1, alpha + n)
    
    return(list( "z" = Z, "rho" = rho, "nm" = n))
  }
  
  ml_yd_ind_prior <- function(y,X,lambda,delta, gm){
    if(sum(delta) > 0){
      delta <- delta == TRUE
      X <-  if(is.null(dim(X))) t(as.matrix(X[delta][!is.na(X[delta])])) else as.matrix(X[,delta]) # Need to convert into proper format if only one obs
      cols <- dim(X)[2];  N <- dim(X)[1]
      a0 <- rep(0,sum(delta))
      A0 <- sqrt(gm) * diag(sum(delta))
      A0_inv <- solve(A0)
      AN <- solve((t(X) %*% X + A0_inv))
      aN <- AN %*% (t(X) %*% y + A0_inv %*% a0)
      SN <- 0.5 * (t(y) %*% y + t(a0) %*% A0_inv %*% a0 - t(aN) %*% solve(AN) %*% aN)
      sN <- 0.5 * N
      
      log_ml_yd <-  -0.5*N*(log(2*pi)) + 0.5 * (log(det(AN))) - 0.5 * (log(det(A0))) + lgamma(sN) - sN * log(SN)
    } else{
      SN <-0.5*(t(y)%*%y)
      sN <- 0.5 * N
      
      log_ml_yd <-  -0.5*N*(log(2*pi)) + lgamma(sN) - sN * log(SN)
    }
    
    return(log_ml_yd)
  }
  
  get_aN_AN_sn_SN <- function(y,X,delta, gm){
    N <- gm
    yc <- y-mean(y)
    
    if(sum(delta)>0){
      delta <- delta == TRUE
      X <-  if(is.null(dim(X))) t(as.matrix(X[delta][!is.na(X[delta])])) else as.matrix(X[,delta]) # Need to convert into proper format if only one obs
      cols <- dim(X)[2];  N <- dim(X)[1]
      a0 <- rep(0,sum(delta))
      A0 <- sqrt(gm) * diag(sum(delta))
      A0_inv <- solve(A0)
      AN <- solve((t(X) %*% X + A0_inv))
      aN <- AN %*% (t(X) %*% y + A0_inv %*% a0)
      SN <- 0.5 * (t(y) %*% y + t(a0) %*% A0_inv %*% a0 - t(aN) %*% solve(AN) %*% aN)
      sN <- 0.5 * N
    }
    else{
      AN_d <- matrix(1)
      aN_d <- c(1)
      SN <-0.5*(t(y)%*%y)
      sN <- 0.5 * N
    }
    return(list('aN' = aN, 'AN' = AN, 'sN' = sN, 'SN' = SN))
  }
  
  ml_yd_g_prior <- function(y,X,lambda,delta, gm){
    if(sum(delta) > 0){
      delta <- delta == TRUE
      X <-  if(is.null(dim(X))) t(as.matrix(X[delta][!is.na(X[delta])])) else as.matrix(X[,delta]) # Need to convert into proper format if only one obs
      cols <- dim(X)[2];  N <- dim(X)[1]
      a0 <- rep(0,sum(delta))
      A0 <- sqrt(gm) * (solve(t(X) %*% X + diag(sum(delta))))
      A0_inv <- solve(A0)
      AN <- solve((t(X) %*% X + A0_inv))
      aN <- AN %*% (t(X) %*% y + A0_inv %*% a0)
      SN <- 0.5 * (t(y) %*% y + t(a0) %*% A0_inv %*% a0 - t(aN) %*% solve(AN) %*% aN)
      sN <- 0.5 * N
      
      log_ml_yd <-  -0.5*N*(log(2*pi)) + 0.5 * (log(det(AN))) - 0.5 * (log(det(A0))) + lgamma(sN) - sN * log(SN)
    } else{
      SN <-0.5*(t(y)%*%y)
      sN <- 0.5 * N
      
      log_ml_yd <-  -0.5*N*(log(2*pi)) + lgamma(sN) - sN * log(SN)
    }
    
    return(log_ml_yd)
  }
  
  get_aN_AN_sn_SN_gprior <- function(y,X,delta, gm){
    N <- gm
    yc <- y-mean(y)
    
    if(sum(delta)>0){
      delta <- delta == TRUE
      X <-  if(is.null(dim(X))) t(as.matrix(X[delta][!is.na(X[delta])])) else as.matrix(X[,delta]) # Need to convert into proper format if only one obs
      cols <- dim(X)[2];  N <- dim(X)[1]
      a0 <- rep(0,sum(delta))
      A0 <- sqrt(gm) * (solve(t(X) %*% X + diag(sum(delta))))
      A0_inv <- solve(A0)
      AN <- solve((t(X) %*% X + A0_inv))
      aN <- AN %*% (t(X) %*% y + A0_inv %*% a0)
      SN <- 0.5 * (t(y) %*% y + t(a0) %*% A0_inv %*% a0 - t(aN) %*% solve(AN) %*% aN)
      sN <- 0.5 * N
    }
    else{
      AN_d <- matrix(1)
      aN_d <- c(1)
      SN <-0.5*(t(y)%*%y)
      sN <- 0.5 * N
    }
    return(list('aN' = aN, 'AN' = AN, 'sN' = sN, 'SN' = SN))
  }
  
  k <- length(X[1,])
  N <- length(X[,1])
  max_tries <- 10
  
  # Group-specific arrays m = 1,...,M
  Xm <- array(X, dim = c(dim(X)[1],dim(X)[2],groups))
  y_m <- matrix(y, nrow = N, ncol = groups)
  
  # Hyperpriors
  alphas <- rep(2,groups)
  a0 <- b0 <- 0 # Uninformative Jeffrey Prior
  
  q <- 1
  while(q < max_tries){
    # Initialization
    rho <- rep(1/groups,groups)
    z <- c(1,2,1,2)
    beta <- matrix(data = NA, nrow = k, ncol = groups)
    for(i in 1:groups) beta[,i] <- runif(k,-q,q)
    beta_hat <- matrix(data = NA, nrow = k, ncol = groups)
    delta <- matrix(data = 1, nrow = k, ncol = groups); delta <- delta == 1
    sigma2 <- 1/rgamma(groups,10,10)
    omega <- matrix(runif(k * groups), nrow = k, ncol = groups)
    
    # Storage matrices
    z_draws <- matrix(data = NA, nrow = iter, ncol = N)
    rho_draws <- matrix(data = NA, nrow = iter, ncol = groups)
    delta_draws <- array(data = NA, dim = c(iter, k, groups))
    beta_draws <- array(data = NA, dim = c(iter, k, groups))
    test_beta_1 <- matrix(data = NA, nrow = iter, ncol = k)
    sigma_draws <- matrix(data = NA, nrow = iter, ncol = groups)
    omega_draws <- array(data = NA, dim = c(iter, k, groups))
    z_omega_draws <- array(data = NA, dim = c(iter, k, groups))
    
    # Testing
    Rj_draws <- array(data = NA, dim = c(iter, k, groups))
    for(j in 1:iter){
      if(print_iter == TRUE) if (j %%1000 == 0) cat("iter =", j, "\r")
      # Step 1 and 2: Sample z and rho
      if(length(unique(z)) < groups) break # Stop algorithm if trapped at local minimum with zero observations
      
      z_rho <- sample_Z_rho(y,X,beta,sigma2,rho,alphas)
      z <- z_rho$z; rho <- z_rho$rho; qm <- colSums(delta); nm <- z_rho$nm; gm <- nm; am <- qm + nm + a0
      if(sum(gm == 0) > 0) break
      z_draws[j,] <- z
      rho_draws[j,] <- rho
      
      # Step 3,4,5: Sample deltas, betas, sigmas
      Xmd <- list(); ym <- list(); beta_hat <- list(); beta_d <- list(); bm <- rep(NA,groups)
      Sigma <- list(); mu <- list()
      for(i in 1:groups){
        ym[[i]] <- y_m[,i][z == i]
        # Step 3: Sampling delta
        perm <- sample(1:k,k)
        for(t in 1:k){
          log_yd <- rep(NA,2)
          log_delta <- rep(NA,2)
          log_post_delta <- rep(NA,2)
          for(d in 0:1){
            delta[,i][perm[t]] <- d == 1
            if(sum(z == i) > 0){
              log_yd[d+1] <- ml_yd_ind_prior(ym[[i]], Xm[,,i][z == i,], lambda, delta[,i], gm[i])
            } else{
              log_yd[d+1] <- 0
            }
          }
          R_j <- exp(log_yd[1] - max(log_yd))/exp(log_yd[2] - max(log_yd))
          post_delta <- c(1-1/(1+R_j*((1-omega[perm[t],i])/omega[perm[t],i])),1/(1+R_j*((1-omega[perm[t],i])/omega[perm[t],i])))
          delta[,i][perm[t]] <- sample(x=c(0,1),1, replace = TRUE, prob = post_delta)
          delta <- delta == TRUE # Convert to TRUE/FALSE
          Rj_draws[j,perm[t],i] <- R_j # Testing
        }
        delta_draws[j,,i] <- delta[,i]
        
        if(sum(delta[,i]) > 1) {
          Xmd[[i]] <- as.matrix(Xm[,,i][,delta[,i]][z == i,])
        } else if(sum(delta[,i] == 1)){
          Xmd[[i]] <- as.matrix(Xm[,,i][,delta[,i]][z == i])
        } else{
          if(sum(z==i) > 1){
            Xmd[[i]] <- as.matrix(Xm[,,i][z==i,][,delta[,i]])
          } else{
            Xmd[[i]] <- as.matrix(Xm[,,i][z==i,][delta[i]])
          }
        }
        if(nm[i] == 1) Xmd[[i]] <- t(Xmd[[i]])
        
        # Sample omegas
        omega[,i] <- rbeta(length(delta[,i]),delta[,i] + 0.5, 1 + 0.5 - delta[,i])
        omega_draws[j,,i] <- omega[,i]
        
        # Auxiliary step: get aN, AN, sN, SN
        c <- get_aN_AN_sn_SN(ym[[i]], Xm[,,i][z == i,], delta[,i], gm[i])
        aN <- c$aN; AN <- c$AN; sN <- c$sN; SN <- c$SN
        
        # Sample sigma2
        sigma2[i] <- if(sum(z ==i) > 0) 1/rgamma(1,sN, SN) else NA
        sigma_draws[j,i] <- sigma2[i]
        
        # Step 4: Sampling beta
        if(sum(delta[,i]) > 0){
          Sigma[[i]] <-  sigma2[i] * AN
          mu[[i]] <- aN
        } else{
          Sigma[[i]] <-  NA
          mu[[i]] <- NA
        }
        
        
        if(sum(delta[,i]) > 0){ 
          beta_d[[i]] <- t(as.matrix(rmvnorm(1,mean = mu[[i]], sigma = Sigma[[i]])))
          beta[,i][delta[,i]] <- beta_d[[i]]
          beta[,i][delta[,i] == FALSE] <- 0
        } else{
          beta[,i][delta[,i] == FALSE] <- 0
        }
        
        beta_draws[j,,i] <- beta[,i]
        
      }
    }
    t <- q - 1
    q <- if(j == iter) max_tries else q + 1
  }
  if(t > 0) cat('Note: The sampler got stuck', t, 'time(s) during the process with zero 
                observations allocated to one of the groups, and was restarted to provide convergence to the desired number of groups')
  draws <- list('betas' = beta_draws,'Rj' = Rj_draws, 'deltas' = delta_draws, 'z' = data.frame('z_draws' = z_draws),
                'sigmas' = sigma_draws)
  return(draws)
  
}