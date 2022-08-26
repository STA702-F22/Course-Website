
################# Gibbs Sampler Function #################
rm(list = ls())
fit_mixture <- function(x,k,n_iter,burn_in){
  ###Initialize
  p = ncol(x)
  n = nrow(x)
  mu_0 = matrix(0,p,1)
  m = 1
  MU_all_k = (matrix(seq(-3,3,length.out =k),k,p)*matrix(rep(sqrt(diag(cov(x))),each=k),k,p))+
    matrix(colMeans(x),k,p,byrow = TRUE)
  #MU_all_k = matrix(0,k,p)
  Phi = diag(p)
  nu_0 = p+2
  SIGMA_all_k = matrix(1000*diag(p),p^2,k)
  z = matrix(sample(1,k,replace=TRUE))
  post_prob_z = matrix(0,n,k)
  a_alpha = b_alpha = 0.25
  alpha = 1
  U = matrix(rbeta(k,1,alpha),nrow=k)
  U[k]=1
  one_min_U = 1L-U
  one_min_U = c(1,cumprod(one_min_U[1:(k-1)]))
  lambda = U*one_min_U
  
  ###Lastly,
  LAMBDA = matrix(0,nrow=(n_iter-burn_in),ncol=k)
  MU = matrix(0,nrow=(n_iter-burn_in),ncol=k*p)
  SIGMA = matrix(0,nrow=(n_iter-burn_in),ncol=k*p*p)
  Z = matrix(0,nrow=(n_iter-burn_in),ncol=n)
  K = matrix(0,nrow=(n_iter-burn_in),ncol=1)
  
  ###Start Gibbs
  for(mc in 1:n_iter){
    #Sample z_i, the mixture i.d. for each crab
    for(j in 1:k){
      post_prob_z[,j] = lambda[j]*dmvnorm(x,mean=MU_all_k[j,],
                                         sigma=matrix(SIGMA_all_k[,j],p,p))
    }
    post_prob_z = post_prob_z/matrix(rowSums(post_prob_z),nrow=n,ncol=k)
    Ran_unif_z = runif(nrow(post_prob_z))
    cumul_z = post_prob_z%*%upper.tri(diag(ncol(post_prob_z)),diag=TRUE)
    z = rowSums(Ran_unif_z>cumul_z) + 1L
    
    #Sample mu and sigma for each mixture component
    for(jj in 1:k){
      if(c(table(factor(z,levels=c(1:k))))[jj]!=0){
        n_j = c(table(factor(z,levels=c(1:k))))[jj]
        x_bar_j = colMeans(matrix(x[which(z==jj),],ncol=p))
        m_n = m + n_j
        nu_n = nu_0 + n_j
        S = ((t(matrix(x[which(z==jj),],ncol=p)) -
                x_bar_j)%*%t(t(matrix(x[which(z==jj),],ncol=p)) - x_bar_j))/n_j
        x_bar_min_mu_0 = (t(t(x_bar_j))-mu_0)%*%t(t(t(x_bar_j))-mu_0)
        Phi_n = Phi + (n_j*S) + (n_j*m*x_bar_min_mu_0/m_n)
        SIGMA_all_k[,jj] = solve(matrix(rWishart(1,nu_n,solve(Phi_n)),p,p))
        mu_n = ((n_j*x_bar_j) + (m*mu_0))/m_n
        MU_all_k[jj,] = rmvnorm(1,mu_n,(matrix(SIGMA_all_k[,jj],p,p)/m_n))
      } 
    }
    #Constraint: ad-hoc ordering to avoid label switching
    mu_ordering = order(MU_all_k[,1],decreasing=TRUE)
    MU_all_k =  matrix(MU_all_k[mu_ordering,],nrow=k) 
    SIGMA_all_k = matrix(SIGMA_all_k[,mu_ordering],ncol=k)
    
    #Sample U and lambda
    n_f = matrix(summary(factor(z,levels=c(1:k))),ncol=1)
    U[k]=1
    U[1:(k-1),1] = rbeta((k-1),(1L+n_f[1:(k-1)]),(alpha+(n-cumsum(n_f[-k]))))
    if(length(which(U[-k]==1))>0){
      U[which(U[-k]==1)] = 0.99999
    }
    one_min_U = 1L-U
    one_min_U_prod = c(1,cumprod(one_min_U[1:(k-1)]))
    lambda = U*one_min_U_prod
    
    #Sample alpha
    alpha = rgamma(1,shape=(a_alpha+k-1),rate=(b_alpha-log(lambda[k])))
    
    #Print iteration info
    if(sum(mc == seq(1,n_iter,by=500))==1){
      cat("Iteration =", mc,"\t"," ",
          "k =",formatC(length(unique(z)), width=2, flag=" "),"\t",
          "alpha =",alpha,"\n",sep = " ")
    } else if(mc==n_iter){
      cat("Iteration =", mc,"\t"," ",
          "k =",formatC(length(unique(z)), width=2, flag=" "),"\t",
          "alpha =",alpha,"\n",sep = " ")
    }
    
    #Sace results
    if(mc > burn_in){
      LAMBDA[(mc-burn_in),] = lambda 
      MU[(mc-burn_in),] = c(t(MU_all_k))
      SIGMA[(mc-burn_in),] = c(SIGMA_all_k) 
      Z[(mc-burn_in),] = z
      K[(mc-burn_in)] = length(unique(z)) }
    
  }
  
  #Return results
  list(LAMBDA=LAMBDA,MU=MU,SIGMA=SIGMA,Z=Z,K=K)
}
################# End of Gibbs Sampler Function #################




#################################### START ####################################
library(MASS)
library(dplyr)
library(DirichletReg)
library(coda)
library(mvtnorm)

################# Load Data
data(crabs)

################# Fit model for Question 1
#blue_crabs = crabs %>% filter(sp == 'B')
#Data = as.matrix(blue_crabs[,c("RW","CL")])
Data = as.matrix(crabs[,c("RW","CL","FL","BD","CW")])
k = 10
n_iter = 10000
burn_in = 0.3*n_iter
set.seed(4321)
ModelFit <- fit_mixture(x=Data,k=k,n_iter=n_iter,burn_in=burn_in)
p <- ncol(Data)
n <- nrow(Data)

###Posterior density for number of occupied clusters
plot(density(ModelFit$K))

###Posterior means for occupied clusters only
LAMBDA_occ = round(colMeans(ModelFit$LAMBDA),2)
occ_index = which(LAMBDA_occ>0)
LAMBDA_occ = matrix(LAMBDA_occ[occ_index],ncol=1)
MU_occ = matrix(round(colMeans(ModelFit$MU),2),ncol=p,byrow = T)
MU_occ = matrix(MU_occ[occ_index,],ncol=p)
SIGMA_occ = matrix(round(colMeans(ModelFit$SIGMA),2),ncol=p^2,byrow = T)
SIGMA_occ = matrix(SIGMA_occ[occ_index,],ncol=p^2)

###For each individual, select posterior mode for cluster id
z_i = apply(ModelFit$Z, 2, function(x) which(as.data.frame(table(factor(x,levels=c(1:k))))$
                                      Freq==max(as.data.frame(table(x))$Freq)))

###For each individual, calculate probability of being in each cluster
z_post_pr = t(apply(ModelFit$Z, 2, function(x) (as.data.frame(table(factor(x,levels=c(1:k))))$Freq/
                                         (n_iter-burn_in))))
z_post_pr = round(z_post_pr,3)
z_post_pr = z_post_pr[,occ_index]

###Plot the probability of being in each occupied cluster
for(h in 1:ncol(z_post_pr)){
  plot(z_post_pr[,h],xlab="Crab Index",ylab="Probability",ylim=c(0,1),main=
         paste0("Post. Inclusion Prob. of being in Cluster ",h),
       col = c(rep("cyan2",50),rep("blue4",50),rep("yellow",50),rep("red",50)),
       pch=17,cex=0.8,lwd=1.5)
  legend(170,0.9,c("B//M","B//F","O//M","O//F"),pch=c(17,17),
         col=c("cyan2","blue4","yellow","red"),bty="n")
}

###Plot the posterior mode for cluster id for all observations
y_lim = c(min(sort(unique(z_i))),max(sort(unique(z_i))))
plot(z_i,xlab="Crab Index",ylab="Cluster",main=
       "Posterior Classification",ylim=y_lim,yaxt="n",
     col = c(rep("cyan2",50),rep("blue4",50),rep("yellow",50),rep("red",50)),
     pch=17,cex=0.8,lwd=1.5)
axis(2, at = sort(unique(z_i)), labels=c(1:4))
legend(170,max(sort(unique(z_i))),c("B//M","B//F","O//M","O//F"),pch=17,
       col=c("cyan2","blue4","yellow","red"),bty="n")





