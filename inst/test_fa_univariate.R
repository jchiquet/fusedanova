rm(list=ls())
library(fusedanova)
data(aves)

fa <- fusedanova(aves$weight, aves$family, gamma = 0)

print(fa$path)

plot(as.hclust(fa), main = "fused-ANOVA tree")
