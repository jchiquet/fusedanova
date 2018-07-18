rm(list=ls())
library(fusedanova)

set.seed(1111)
n1 <- 100; n2 <- 200; n3 <- 150; n <- n1 + n2 + n3
mu1 <- -2; mu2 <- 0; mu3 <- 2
sigma <- .25
trait <- setNames(
  c(rnorm(n1, mu1, sigma), rnorm(n2, mu2, sigma), rnorm(n3, mu3, sigma)),
  nm = paste0("ind",1:n)
  )
group <- rep(1:3, c(n1,n2,n3))

# sampling ther data order
o <- sample(1:(n1 + n2 + n3))
group <- group[o]
trait <- trait[o]

fa <- fusedanova(trait, gamma = 1)
plot(fa)

clusters <- cutree(fa$hc, 3)
clusters <- lapply(split(names(clusters), clusters), sort)

laricode::ARI(cutree(fa$hc, 3), group)
