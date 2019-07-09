rm(list=ls())
library(fusedanova)

## ward.D2 link function 
inertia <- function(gr1){
  sum((gr1 - mean(gr1))^2)
}
ward.D2.link <- function(gr1, gr2){
  inertia(c(gr1, gr2)) - inertia(gr1) - inertia(gr2)
}

## data: random sorted
x <- sort(rnorm(10))
or <- hclust(dist(x), method = "ward.D2")
plot(or)
or$height

w1 <- fusedanova:::ward1d.numeric(as.double(x))
plot(w1)
w1$height

lw <- ward.D2.link(x[7:8], x[9:10]); sqrt(lw*2)

print(sum((or$height - w1$height)^2))
