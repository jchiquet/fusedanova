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

loglik_ANOVA <- function(group, x) {
  n  <- length(x)
  nk <- tabulate(group)
  k  <- length(nk)
  betak <- rowsum(x, group)/nk
  RSS <- sum((x - betak[group])^2)
  sigma2 <- RSS/(n - k)
  ### Check the loglikelihood
  loglik <- -.5 * (n*log(2*pi) + n*sum(log(sigma2)) + (n - k))
  loglik
}
