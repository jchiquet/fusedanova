stopif <- function(expr, message) {
  if (expr) stop(message)  
}

get_norm <- function(x, group, n, k, nk) {
  if (k != n) {
    # if any initial grouping, normalize withing each group
    s <- (rowsum(x^2,group) - (1/nk) * (rowsum(x,group))^2) / (nk - 1)
    s[nk == 1] <- 0
    s <- sqrt(sum(s*(nk - 1))/(n - k))
  } else {
    # if no grouping (one guy per class0) normalize at the vector scale
    s <- sd(x)
  }
  s
}
