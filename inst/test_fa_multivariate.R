rm(list=ls())
library(fusedanova)
library(tidyverse)
data(iris)

fa_mult <- fusedanova(select(iris, Sepal.Length, Petal.Length))
fa_uni1 <- fusedanova(pull(iris, Sepal.Length))
fa_uni2 <- fusedanova(pull(iris, Petal.Length))

par(mfrow = c(1,3))
plot(as.hclust(fa_uni1), main = "Sepal")
plot(as.hclust(fa_uni2), main = "Petal")
plot(fa_mult, main = "both")
par(mfrow = c(1,1))
