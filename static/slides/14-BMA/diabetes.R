set.seed(8675309)
source("yX.diabetes.train.txt")
diabetes.train = as.data.frame(diabetes.train)

library(BAS)
diabetes.bas = bas.lm(y ~ ., data=diabetes.train, method="MCMC", n.models = 2500000, thin = 10, initprobs="eplogp")

diagnostics(diabetes.bas, type="pip")

summary(diabetes.bas$sampleprobs)

image(diabetes.bas)

source("yX.diabetes.test.txt")
diabetes.test = as.data.frame(diabetes.test)
colnames(diabetes.test)[1] = "y"

out = predict(diabetes.bas, newdata=diabetes.test, estimator="BMA")
cv.summary.bas(out$fit, diabetes.test$y)

library(lars)
diabetes.lasso = lars(as.matrix(diabetes.train[, -1]), diabetes.train[,1], type="lasso")

plot(diabetes.lasso)
best = which.min(diabetes.lasso$Cp)
out.lasso = predict(diabetes.lasso, as.matrix(diabetes.test[,-1]))
cv.summary.bas(out.lasso$fit[, best], diabetes.test$y)

