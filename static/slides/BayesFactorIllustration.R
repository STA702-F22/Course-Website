
theta <- 0.25
n<-1000
y<-rbinom(n,1,prob=theta)
logbayes<-c()
for (i in 1:n){
  #_sum <- i*theta
  #y_sum <- sum(rbinom(i,1,prob=theta))
  temp_logbayes <- (lbeta(1+y_sum, 1+i-y_sum) - i*log(1/2))
  logbayes <- c(logbayes,temp_logbayes)
}
plot(x=c(1:n), y=exp(logbayes))
