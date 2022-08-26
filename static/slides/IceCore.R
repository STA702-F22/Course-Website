
###########################################################################
###########################################################################
################################# Ice core ################################
###########################################################################
###########################################################################

###### Clear environment and load libraries
rm(list = ls())
library(coda)
library(mvtnorm)
#library(psych)


###### Data
#Load data
Data <- read.table("data/icecore.txt",header=TRUE)


###### Now to the sampler
#Data summaries
n <- nrow(Data)
Y <- (Data$tmp-mean(Data$tmp))/sqrt(var(Data$tmp))
X <- cbind(rep(1,n),(Data$co2-mean(Data$co2))/sqrt(var(Data$co2)))
DY <- abs(outer((1:n),(1:n) ,"-"))


#Hyperparameters for the prior
nu_0 <- 1
sigma_0_sq <- 1
mu_0 <- matrix(c(0,0),ncol=1)
Lambda_0 <- diag(1000,nrow=2)
Lambda_0_inv <- solve(Lambda_0)


#Set parameters for proposal density
delta <- 0.1 #use it to tune acceptance ratio


#Initial values for sampler
lmfit <- lm(Y~-1+X) #intercept already in X
beta <- lmfit$coef
sigma_sq <- summary(lmfit)$sigma^2
#fit.gls <- gls(Y ~ X[,2], correlation=corARMA(p=1),method="ML")
rho <- acf(lmfit$res,plot=FALSE)$acf[2]


#First set number of iterations and burn-in, thinning, then set seed
n_iter <- 1000
burn_in <- 0.1*n_iter
thin <- 1
#set.seed(1234)
set.seed(121372)


#Set counter for acceptances
accept_counter <- 0


#Set null matrices to save samples
THETA <- matrix(0,nrow=(n_iter/thin),ncol=4)


#Now, to the Gibbs sampler
for(s in 1:(n_iter+burn_in)){
  
  #sample beta
  C_p <- rho^DY
  C_p_inv <- solve(C_p)
  Lambda_n <- solve( t(X)%*%C_p_inv%*%X/sigma_sq + Lambda_0_inv)
  mu_n <- Lambda_n%*%( t(X)%*%C_p_inv%*%Y/sigma_sq + Lambda_0_inv%*%mu_0)
  beta <- t(rmvnorm(1,mu_n,Lambda_n))
  
  
  #sample sigma_sq
  nu_n <- nu_0 + n
  SSR_beta_rho <- t(Y-X%*%beta)%*%C_p_inv%*%(Y-X%*%beta)
  nu_n_sigma_n_sq <- nu_0*sigma_0_sq + SSR_beta_rho
  sigma_sq <- 1/rgamma(1,(nu_n/2),(nu_n_sigma_n_sq/2))
  
  
  #sample rho
  rho_star <- abs(runif(1,rho-.1,rho+.1)) #if rho_star < 0, reassign as |rho_star|
  rho_star <- min(rho_star, 2-rho_star) #if rho_star > 1, reassign as (2 - rho_star)
  
  XBeta <- X%*%beta
  Sigma <- sigma_sq*C_p
  Sigma_star <- sigma_sq*(rho_star^DY)
  log_r <- dmvnorm(c(Y),XBeta,Sigma_star,log=T) - dmvnorm(c(Y),XBeta,Sigma,log=T)
  
  #Or equivalently...: need #library(psych) to compute trace
  #log_r <- -0.5*(determinant(rho_star^DY,log=TRUE)$mod -
  #                determinant(rho^DY,log=TRUE)$mod +
  #                tr( (Y-X%*%beta)%*%t(Y-X%*%beta)%*%
  #                      (solve(rho_star^DY) -solve(rho^DY)) )/sigma_sq )
  
  if(log(runif(1)) < log_r){
    accept_counter <- accept_counter + 1
    rho <- rho_star
  }
  
  
  if(s > burn_in && s%%thin==0){
    THETA[(s-burn_in)/thin,] <- c(beta,sigma_sq,rho)
    cat(paste("Iteration: ", s,"\n", sep = ""))
    cat(paste("Current acceptance rate: ", round(accept_counter/s,2), "\n\n", sep = ''))
  }
}


colnames(THETA) <- c("beta_0","beta_1","sigma_sq","rho")
#Check acceptance rate
accept_counter/(n_iter+burn_in)


plot(mcmc(THETA))
autocorr.plot(mcmc(THETA))
#Autocorrelation is especially high for sigma_sq and rho + trace plots are horrid.
#Let's try keeping only every 20th scan and also add a burn-in.
#Also, increase n_iter to 50000
apply(THETA,2,effectiveSize)


round(apply(THETA,2,mean),2)
round(apply(THETA,2,function(x) quantile(x,probs=c(0.025,0.975))),2)
mean(THETA[,"beta_1"]>0)
#Posterior of beta indicates that CO2 measures are positively correlated with temperature.
#The posterior mean is 0.28 with a 95% credible interval of (0.09; 0.49).
#Each 1 increase in CO2 is associated with a 0.28 increase in temperature on average.
#The posterior probability that beta_1 > 0 is 0.999.







