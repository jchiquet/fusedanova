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

# 
# hc_list <- list(as.hclust(fa_uni1), as.hclust(fa_uni2))
# order_hc_list(hc_list)

# 
# rm(list=ls())
# library(fusedanova)
# library(Rmergetrees)
# library(tidyverse)
# data(iris)
# 
# ward1 <- hclust(dist(pull(iris, Sepal.Length)), method = "ward.D2")
# ward2 <- hclust(dist(pull(iris, Petal.Length)), method = "ward.D2")
# 
# ward_merged <- mergeTrees(list(ward1, ward2))
# 
# par(mfrow = c(1,3))
# plot(ward1, main= "Ward")
# plot(ward2, main= "complete")
# plot(ward_merged, main= "merged")
# 
# fa1 <- as.hclust(fusedanova(pull(iris, Sepal.Length)))
# fa2 <- as.hclust(fusedanova(pull(iris, Petal.Length)))
# 
# fa_merged <- mergeTrees(list(fa1, fa2))
# par(mfrow = c(1,3))
# plot(fa1, main= "Ward")
# plot(fa2, main= "complete")
# plot(fa_merged, main= "merged")
