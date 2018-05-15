rm(list=ls())
library(fusedanova)
library(tidyverse)
library(microbenchmark)

n <- 1e5
x <- sort(rnorm(n))
w <- 1:n

timings <- microbenchmark(fusedanova(x, w), times = 20)

autoplot(timings)
summary(timings, unit = "s")
