rm(list=ls())
library(fusedanova)
data(aves)

fa <- fusedanova(aves$weight, aves$family, gamma = .1)

print(fa$path)
par(mfrow=c(1,3))
plot(fa, main="fused-ANOVA new")
plot(fa, "AIC", main="fused-ANOVA new - AIC", type="l")
plot(fa, "BIC", main="fused-ANOVA new - BIC", type="l")
par()
