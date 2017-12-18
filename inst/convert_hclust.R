library(fusedanova)
source("inst/functions_Audrey.R")

set.seed(111)
n <- 10
x <- rnorm(n)
w <- 1:n
 
fa1 <- fusedanova(x, w)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda>0, ]
fa1_path <- fa1_path[order(fa1_path$lambda), ]
hc <- hclust.one.dim(n, fa1)
dendro_fa1 <- reorder(as.dendrogram(hc), fa1@result[[1]]$order)
#dendro_fa1 <- as.dendrogram(hc)

fa2 <- fusedanova2(x, w)

par(mfrow=c(1,2))
plot(dendro_fa1, main ="fa1")
plot(fa2$dendrogram, main="fa2")
# plot(reorder(fa2$dendrogram, fa2$order), main="fa2")
par()

