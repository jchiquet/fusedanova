rm(list=ls())
library(fusedanova)
library(tidyverse)
data(iris)

fa_mult <- fusedanova(select(iris, -Species))
fa_univ <- fusedanova(pull(iris, Sepal.Length))

par(mfrow = c(1,2))
plot(fa_univ)
plot(fa_mult)

