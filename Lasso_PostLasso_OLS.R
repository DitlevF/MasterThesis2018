rm(list=ls())
# Storage
n = 50 #sample size
p = 50 # number of variables
max_s <- 25
iter <- 100
MSE <- array(NA, dim = c(max_s, iter, 3))

require(hdm)

for(s in 1:max_s){
  for(j in 1:iter){
    X = matrix(rnorm(n*p), ncol=p)
    beta = c(rep(5,s), rep(0,p-s))
    Y = X%*%beta + rnorm(n)
    set.seed(s*j*2018)
    
    # Lasso
    ## Using RLasso package
    lasso.reg <- cv.glmnet(X,Y)
    betas_hat <- as.numeric(coef(lasso.reg)[-1])
    
    #lasso.reg = rlasso(Y~X,post=FALSE)
    #betas_hat <- lasso.reg$beta
    MSE[s,j,1] <- sum((betas_hat - beta)^2)
    
    # Post Lasso

    #post.lasso.reg = rlasso(Y~X,post=TRUE) 
    #betas_hat <- post.lasso.reg$beta
    
    indices1 <- which(betas_hat != 0)
    fit.lasso2 <- lm(Y ~X[,indices1])
    betas_hat <- rep(0,length(beta))
    betas_hat[indices1] <- as.numeric(coef(fit.lasso2)[-1])
    
    MSE[s,j,2] <- sum((betas_hat - beta)^2)
    
    # Oracle OLS
    X_oracle <- as.matrix(X[,0:s])
    fit.OLS <- lm.fit(X_oracle,Y)
    betas_OLS_hat <- as.numeric(coef(fit.OLS))
    betas_OLS_hat <- c(betas_OLS_hat,rep(0,p-s))
    
    MSE[s,j,3] <- sum((betas_OLS_hat - beta)^2)
    
  }
  cat("iter =", s,  "\r")
  
}


MSE1 <- as.data.frame(rowMeans(MSE[,,1]))
MSE1 <- cbind(MSE1, 1:max_s)
colnames(MSE1) <- c('MSE1', 's')

MSE2 <- as.data.frame(rowMeans(MSE[,,2]))
colnames(MSE2) <- c('MSE2')

MSE3 <- as.data.frame(rowMeans(MSE[,,3]))
colnames(MSE3) <- c('MSE3')

MSE_total <- cbind(MSE1,MSE2,MSE3)


MSE_PLOT <- ggplot(MSE_total, aes(s)) + 
  geom_line(aes(y = MSE1, linetype = "Lasso"), size = 0.5 ) + 
  geom_line(aes(y = MSE2, linetype = "Post Lasso"), size = 0.5) + 
  geom_line(aes(y = MSE3, linetype = "Oracle OLS"), size = 0.5) + 
  scale_x_continuous(name = "Number of active regressors", limits = c(0,16)) + scale_y_continuous(name = "Bias", limits = c(0,5)) +
  theme(axis.title = element_text(size = 5)) +
  theme(axis.text.x= element_text(size = 3)) + theme(axis.text.y= element_text(size = 3)) +
  theme(legend.text=element_text(size=3)) +
  theme(legend.title=element_blank())

MSE_PLOT

ggsave(file="Github_Final/Updated Marginal Likelihood/post_lasso.png", plot = MSE_PLOT, width=4, height=1.3, dpi=1000)




