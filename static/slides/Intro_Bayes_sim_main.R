
######## Simulation Exercise
rm(list = ls())
PlotPriorPlusPosterior <- function(a,b,ones,zeros){
  curve(dbeta(x,a+ones,b+zeros),col="green4",xlab=expression(theta),ylab="density",lwd=1,ylim=c(0,15))
  curve(dbeta(x,a,b),col="red3",add=TRUE,lwd=1)
  legend("topright", legend=c(expression(paste(pi,"(", phi, "|",x,")")),expression(paste(pi,"(", phi,")"))),
         col=c("green4", "red3"), lwd=2, cex=1)
  
  post_mean <- (a+ones)/(a+ones+b+zeros)
  p <- c(0.75, 7)
  points(t(p), pch=16)
  text(t(p), eval(expression(paste("Post. Mean = ",round(post_mean,3)))), adj=-0.05)
}
PlotPriorPlusPosterior(a=2,b=25,ones=100,zeros=1000)

set.seed(360602)
n <- 300
theta_true <- 0.5
y <- rbinom(n,1,theta_true)
#y <- rep(0,n)
mean(y)

ones <- 0; zeros <- 0
for(i in 1:n){
  if(y[i]==1){
    ones <- ones + 1
  } else {zeros <- zeros + 1}
  PlotPriorPlusPosterior(a=2,b=25,ones,zeros)
  points(t(c(0.75, 8)), pch=16)
  text(t(c(0.75, 8)), paste("Observation",i," = ",y[i]), adj=-0.05)
  if(i < 8){
    Sys.sleep(5)
  }
  if(i > 8 & i < 25){
    Sys.sleep(2)
  }
  if(i > 25){
    Sys.sleep(0.05)
  }
}