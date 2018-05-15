rm(list=ls())
library(fusedanova)
library(tidyverse)
data(iris)

fa <- fusedanova(iris %>% select(-Species))

print(fa$path)
par(mfrow=c(1,3))
plot(fa, main="fused-ANOVA new")
plot(AIC(fa), main="fused-ANOVA new - AIC", type="l")
plot(BIC(fa), main="fused-ANOVA new - BIC", type="l")
par()
