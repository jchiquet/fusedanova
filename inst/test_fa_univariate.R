rm(list=ls())
library(fusedanova)
data(aves)

fa <- fusedanova(aves$weight, aves$family, gamma = 0)

print(fa$path)
par(mfrow = c(1,2))
plot(fa, main = "fused-ANOVA tree")
plot(fa, "AIC", main = "AIC/BIC", type = "l", col= "blue")
lines(BIC(fa),  col= "red")
par(mfrow = c(1,1))
