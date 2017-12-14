library(fusedanova)
source("inst/functions_Audrey.R")
data(aves)

fa1 <- fusedanova(aves$weight, aves$family)
rule_fa1 <- RecuperationRegles(40, fa1)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda != 0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]

fa2 <- fusedanova2(aves$weight, aves$family)
fa2_path <- fa2$path

print(sum((fa1_path- fa2_path[, 1:5])^2))

print(tail(fa2_path, 8))
print(tail(fa1_path, 8))

n <- 1e5
x <- sort(rnorm(n))
w <- 1:n
 
fa1 <- fusedanova(x, w)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda != 0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]

fa2 <- fusedanova2(x, w)
fa2_path <- fa2$path
print(sum((fa1_path- fa2_path[, 1:5])^2))

## New version is faster for sufficiently large data
out <- microbenchmark::microbenchmark(old = fusedanova(x, w), new = fusedanova2(x, w), times = 20)
ggplot2::autoplot(out)
