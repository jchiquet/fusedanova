rm(list=ls())
library(fusedanova)

n <- 1e5
x <- sort(rnorm(n))
w <- 1:n
 
fa1 <- fusedanova_old(x, w, gamma=1)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda != 0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]

fa2 <- fusedanova(x, w)
fa2_path <- fa2$path
print(sum((fa1_path[, c(1,2,4,5)]- fa2_path[, 1:4])^2))

## New version is faster for sufficiently large data
out <- microbenchmark::microbenchmark(old = fusedanova_old(x, w), new = fusedanova(x, w), times = 20)
ggplot2::autoplot(out)
