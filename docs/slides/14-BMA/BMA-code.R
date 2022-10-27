library(HH)
data(usair)
#pairs(usair,diag.panel=panel.hist )
pairs(usair)
poll.lm = lm(log(SO2) ~ temp + log(mfgfirms) + log(popn) + wind + precip+ raindays, data=usair)
par(mfrow=c(2,2))
plot(poll.lm, ask=F)
library(car)
help(av.plots)
par(mfrow=c(2,3))
av.plots(poll.lm, ask=F, one.page=T, iden=F)

#library(devtools)
#install_github("merliseclyde/BAS")
# install from CRAN

library(BAS)
poll.bma = bas.lm(log(SO2) ~ temp + log(mfgfirms) + log(popn) +
                    wind + precip+ raindays, data=usair, prior="g-prior", alpha=41, n.models=2^6, modelprior=uniform())

pdf("poll-bma-sum.pdf")
par(mfrow=c(2,2))
plot(poll.bma, ask=F)
dev.off()

par(mfrow=c(1,1))
pdf("poll-image.pdf")
image(poll.bma)
dev.off()

beta = coef(poll.bma)
pdf("poll-beta.pdf")
par(mfrow=c(2,3))
plot(beta, subset=2:7,ask=F)
dev.off()


poll.bb.bma = bas.lm(log(SO2) ~ temp + log(mfgfirms) + log(popn) + wind + precip+ raindays, data=usair, prior="g-prior", alpha=41, n.models=2^6, modelprior=beta.binomial(1,6), update=50, initprobs="Uniform")
par(mfrow=c(1,1))
pdf("poll-bb-image.pdf")
image(poll.bb.bma)
dev.off()

# pollution & mortality data

#ex1217 = read.table("http://stat.duke.edu/courses/Fall10/sta290/datasets/sleuthdata/ex1217.csv", sep=",", header=T)
library(Sleuth3)
data(ex1217)
attach(ex1217)
ex1217$logHC = log(HC)
ex1217$logNOX = log(NOX)
ex1217$logSO2 = log(SO2)
detach(ex1217)
mortality = data.frame(ex1217[,c(2:14, 18:20)])

mort.bma = bas.lm(Mortality ~ ., data=mortality, prior="ZS-null", alpha=60, n.models=2^15, update=100, initprobs="eplogp")

pdf("mort-sum.pdf")
par(mfrow=c(2,2))
plot(mort.bma, ask=F)
dev.off()

par(mfrow=c(1,1))
pdf("mort-image.pdf")
image(mort.bma)
dev.off()

mort.beta = coef(mort.bma)
pdf("mort-beta1.pdf")
par(mfrow=c(2,4))
plot(mort.beta, subset=1:8, ask=F)
dev.off()
pdf("mort-beta2.pdf")
par(mfrow=c(2,4))
plot(mort.beta, subset=9:16, ask=F)
dev.off()

which.mat = list2matrix.which(mort.bma,1:(2^15))
poll.in = (which.mat[, 14:16] %*% rep(1, 3)) > 0
sum(poll.in * mort.bma$postprob)




SO2.bmaZS = bas.lm(logSO2 ~ PRECIP +HUMIDITY+JANTEMP+JULYTEMP+OVER65+HOUSE+EDUC+SOUND+DENSITY+NONWHITE+WHITECOL+POOR, data=mortality, prior="ZS-null", alpha=60, n.models=2^15, update=100, initprobs="eplogp")
NOX.bmaZS = bas.lm(logNOX ~ PRECIP +HUMIDITY+JANTEMP+JULYTEMP+OVER65+HOUSE+EDUC+SOUND+DENSITY+NONWHITE+WHITECOL+POOR, data=mortality, prior="ZS-null", alpha=60, n.models=2^15, update=100, initprobs="eplogp")
HC.bmaZS = bas.lm(logHC ~ PRECIP +HUMIDITY+JANTEMP+JULYTEMP+OVER65+HOUSE+EDUC+SOUND+DENSITY+NONWHITE+WHITECOL+POOR, data=mortality, prior="ZS-null", alpha=60, n.models=2^15, update=100, initprobs="eplogp")


adjmort.bmaZS = bas.lm(MORTALITY ~ ., data=mortality, prior="ZS-null", alpha=60, n.models=2^15, update=100, initprobs="eplogp", modelprior=Bernoulli(probs=c(1, 1, 1, 1, 1, .5, 1, .5, 1, 1, .5, .5, .5, .5, .5)))


library(BEAU)
ttt = BAC(mortality[,c(1,15, 2:13)], num_its=20000, burn=1000)
plotBAC(ttt)
summaryBAC(ttt)
