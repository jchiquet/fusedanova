rm(list=ls())
library(fusedanova)
library(multifusedanova)
library(dplyr)

y <- c(2,1,5,7,9,16,13)
n <- length(y); K <- n; nk <- 1

fa0 <- fusedanova(y, standardize=FALSE) 
cl_val <- get.fa.vec.classifC(fa0)

## new code
fa  <- fusedanova2(y, standardize=FALSE) 
cl_new  <- fa$cut_tree()


