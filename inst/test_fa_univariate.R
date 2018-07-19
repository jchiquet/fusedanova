rm(list=ls())
library(fusedanova)
data(aves)

fa <- fusedanova(aves$weight, aves$family, gamma = 0)

plot(fa, main = "fused-ANOVA tree")
