library(fusedanova)
source("inst/functions_Audrey.R")

set.seed(111)
n <- 10
x <- sort(rnorm(n))
w <- 1:n
 
fa1 <- fusedanova(x, w)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda>0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]
hc <- hclust.one.dim(n, fa1)
plot(as.dendrogram(hc), type = "triangle", center=TRUE, ylab = "penalty")

fa2 <- fusedanova2(x, w)
fa2_path <- fa2$path

as.merge <- function(path) {
  n <- nrow(path) + 1
  last_time_seen <- rep(NA,n)
  merge <- -path[, c("down", "high")]
  for (i in 1:(n - 1)) {

    # print(last_time_seen)    
    if (!is.na(last_time_seen[abs(merge[i, 1])]))
      merge[i, 1] <- last_time_seen[abs(merge[i, 1])]
    
    if (!is.na(last_time_seen[abs(merge[i, 2])]))
      merge[i, 2] <- last_time_seen[abs(merge[i, 2])]

    last_time_seen[as.numeric(abs(merge[i, ]))] <- i
    
  }
  t(apply(merge, 1, sort))
}