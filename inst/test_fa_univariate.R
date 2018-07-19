rm(list=ls())
library(fusedanova)
data(aves)

fa <- fusedanova(aves$weight, aves$family, gamma = 0)

print(fa$path)
print(as.hclust(fa))

plot(as.hclust(fa), main = "fused-ANOVA tree")
