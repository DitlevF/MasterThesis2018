
ss_NMIG <- function(X,Y, r = 1/10000, iter = 1000, print_iter = TRUE){
  
  # Define location scale transformation of t-distribution function
  dt_ls <- function(x,df = 1000, mean = 0, sd = 1) 1/sd * dt((x - mean)/sd, df)
  
  k <- length(X[1,])
  N <- length(X[,1])
  
  # For now, let hyperpriors be fixed within the function
  a <- 1; b <- 1 # Omega hyperpriors
  v <- 5; Q <- 4 # Psi hyperpriors (from paper)
  yc <- Y - mean(Y)
  
  # Storage matrices
  delta_draws <- matrix( data = NA, nrow = iter, ncol = k)
  psi_draws <- matrix( data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  mu_draws <- rep(NA, iter)
  alpha_draws <- matrix(data = NA, nrow = iter, ncol = k)
  omega_draws <- matrix(data = NA, nrow = iter, ncol = k)

  # Initialize
  sigma2 <- 1/rgamma(1,1,1)
  delta <- rep(FALSE,k)
  psi <- 1/rgamma(k, shape = v, rate = Q)
  rdelta <- rep(NA,k)
  alpha <- runif(k,-0.5,0.5)
  omega <- rbeta(1,a,b)
  
  for(i in 1:iter){
    if(print_iter == TRUE) if (i %%100 == 0) cat("iter =", i, "\r")
  
    # Step 1: Sample mu
    mu <- rnorm(1,yc, sigma2/N)
    mu_draws[i] <- mu
    
    # Step 2: deltas and psis
    for(j in 1:k){
      
      # 2a: Deltas
      Lj <- dt_ls(alpha[j], df = 2*v, mean = 0, sd = sqrt(r*Q/v)) / dt_ls(alpha[j], df = 2*v, mean = 0, sd = sqrt(Q/v))
      post_delta <- c(1 - 1/(1+((1-omega)/omega)*Lj), 1/(1+((1-omega)/omega)*Lj))
      delta[j] <- sample(x=c(0,1),1, replace = TRUE, prob = post_delta)
      rdelta[j] <- if(delta[j] == 1) 1 else r
      
      # 2b: Psis
      psi[j] <- 1/rgamma(1, shape = v + 0.5, rate = Q + alpha[j]^2/(2*rdelta[j]))
    }
    
    delta_draws[i,] <- delta
    psi_draws[i,] <- psi
    
    # Step 3: Sample omega
    omega <- rbeta(1, a + sum(delta), b + k - sum(delta))
    omega_draws[i,] <- omega 
    
    # Step 4: Sample alphas
    D <- diag(rdelta * psi)
    AN <- solve(t(X)%*% X * (1/sigma2) + solve(D))
    aN <- AN %*% t(X) %*% yc * (1/sigma2)
    alpha <- rmvnorm(1, mean = aN, sigma = AN)
    alpha_draws[i,] <- alpha

    # Step 5: Sample sigmas
    sN <- (N-1)/2
    SN <- 0.5 * t((yc - X %*% t(alpha))) %*% (yc - X %*% t(alpha))
  
    sigma2 <- 1/rgamma(1, sN, SN)
    sigma_draws[i] <- sigma2
  }
  
  draws <- list('alphas' = data.frame('alphas' = alpha_draws), 'deltas' = data.frame('deltas' = delta_draws), 
                'mu' = mu_draws, 'sigma2' = sigma_draws, 'omegas' = data.frame('omegas' = omega_draws))
                 
  colnames(draws$deltas) <- colnames(draws$alphas) <- colnames(X)
  return(draws)
  
}