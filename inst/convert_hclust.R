library(fusedanova)
source("inst/functions_Audrey.R")

set.seed(111)
n <- 10
x <- rnorm(n)
w <- 1:n
gamma <- 3
weights <- "laplace"
standardize <- FALSE

fa1 <- fusedanova(x, w, weights = weights, gamma = gamma, standardize = standardize)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda>0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]
hc <- hclust.one.dim(n, fa1)

fa2 <- fusedanova2(x, w, weighting = weights, gamma = gamma, standardize = standardize)

print(fa2$path)
print(fa1_path)
par(mfrow=c(1,2))
plot(hc, main ="fa1")
plot(fa2$hc, main="fa2")
par()

