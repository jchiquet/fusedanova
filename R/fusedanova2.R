##' @export 
fusedanova2 <- function(x, group = factor(rep(1,length(x))),
                       weighting = c("laplace", "gaussian", "adaptive"),
                       gamma = 0, standardize = FALSE, W = NULL) {
  
  ## overwrite default parametrs with user's
  weighting <- match.arg(weighting)
  if (!is.null(W)) weighting <- "personal" else W <- matrix(0,0,0)
  # conversion of group ot a factor
  if (!is.factor(group)) group <- as.factor(group)
    
  ## problem dimension
  n  <- length(x)
  k  <- length(unique(group))
  nk <- tabulate(group)
  
  ## check to partilly avoid crashes of th C++ code
  stopif(!is.numeric(x)     , "x must be a numeric vector.")
  stopif(gamma < 0          , "gamma must be non-negative.")
  stopif(any(is.na(x))      , "NA value in x not allowed.")
  stopif(n != length(group) ,  "x and group length do not match")
  stopif(k == 1             , "x has only one level, and there is no point in fusing one group, isn't it?")
  if (weighting == "personal") stopif(nrow(W) != n | nrcol(W) != n, "W must be a square matrix.")
  
  # data standardization
  if (standardize) {
    if (k != n) {
      # if any initial grouping, normalize withing each group
      s <- (rowsum(x^2,group) - (1/nk) * (rowsum(x,group))^2) / (nk - 1)
      s[nk == 1] <- 0
      s <- sqrt(sum(s*(nk - 1))/(n - k))
    } else {
      # if no grouping (one guy per class0) normalize at the vector scale
      s <- sd(x)
    }
    x <- (x - mean(x))/s
  }
  
  # data compression and ordering
  mean_k <- rowsum(x, group)/nk
  order  <- order(mean_k) 

  slopes <- get_slopes(mean_k[order], nk[order], weighting, gamma, W)
  out    <- fuse(mean_k[order], slopes, nk[order]) 

  hc <- list(merge = out$merge, height = out$path$lambda, labels = levels(group)[order], order = out$order)
  class(hc) <- "hclust"

  res <- list(path = out$path, hc = hc)
  res
}

