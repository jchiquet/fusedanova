rm(list=ls())
library(fusedanova)
data(aves)
source("inst/functions_Audrey.R")

set.seed(111)
# x <- aves$weight
# group <- aves$family
# n <- length(unique(group))
n <- 10
x <- rnorm(n)
group <- 1:n

gamma <- 1
weights <- "laplace"
standardize <- TRUE

fa1 <- fusedanova_old(x, group, weights = weights, gamma = gamma, standardize = standardize)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda>0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]
hc <- hclust.one.dim(n, fa1)

fa2 <- fusedanova(x, group, weighting = weights, gamma = gamma, standardize = standardize)

print(fa2$path)
print(fa1_path)
par(mfrow=c(1,2))
plot(hc, main ="fused-ANOVA old")
plot(fa2, main="fused-ANOVA new")
par()

