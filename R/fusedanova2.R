##' Fit a Fused ANOVA model
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total Variation) penality, 
##' ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the difference of each coefficient. 
##' 
##' @param x vector of observation for n individuals.
##'
##' @param class vector or factor giving the initial class of each individual. If missing, 
##' \code{1:length(x)} is used (clustering mode with one individual per class).
##'
##' @param weights character; which type of weights is supposed to be used.
##' The supported weights are: \code{"default"}, \code{"laplace"}, \code{"gaussian"},
##'  \code{"adaptive"}, \code{"naivettest"}, \code{"ttest"}, \code{"welch"} and \code{"personal"}. 
##' See details below. By default, its value is \code{"default"}.
##' 
##' @param standardize logical; should the vector be standardized before computation?
##' Default is \code{TRUE}.
##' 
##' @param ... list of additional parameters to control the optimization procedure. Include :
##' \itemize{%
##' \item{\code{gamma}: } {non-negative scalar; the \eqn{\gamma}{gamma} parameter needed for
##' \code{"laplace"}, \code{"gaussian"} and \code{"adaptive"} weights. Default is 1.}
##'
##' \item{\code{W}: } {numeric matrix; the matrix of weights needed if the \code{"personal"}
##' weights were selected. By default, a matrix with zero row and zero column.}
##' 
##' \item{\code{checkargs}: }{logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is \code{TRUE}.
##' Automatically set to \code{FALSE} when a call is made
##' from cross-validation}
##'
##' \item{\code{verbose}: } {boolean; should the code print out its progress. By default, FALSE. }
##'
##' \item{\code{mxSplitSize}: } {integer; the maximum size for a group for being checked
## for eventual splitting. By default, 100.
##' the cores available. }
##'
##' }
##'
##' @return an object with class \code{fusedanova}, see the
##' documentation page \code{\linkS4class{fusedanova}} for details.
##'
##' The optimization problem solved by fused-ANOVA is
##' \if{latex}{\deqn{%
##' \hat{\beta}_{\lambda} = \arg \min_{\beta}
##' \left\{\sum_{k=1}^K \sum_{i=1}^{n_k} \left(Y_{ik}-\beta_k \right)^2
##' + \lambda \sum_{k,\ell} w_{kl} \left|\beta_k - \beta_\ell \right|\right\}}}
##' \if{html}{\out{ <center> &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub></sub> =
##' argmin<sub>&beta;</sub> sum<sub>k</sub> sum_i (Y<sub>ik</sub> - &beta<sub>k</sub>)<sup>2</sup>
##' + &lambda; sum<sub>k,l</sub> w<sub>k,l</sub>
##' &#124; &beta;<sub>k</sub> - &beta;<sub>l</sub> &#124;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda) = argmin_beta sum_k sum_i (Y_ik - beta_k)^2
##' + lambda sum_k sum_l w_kl | beta_k - beta_l|,}}
##'
##' where \eqn{Y_{ik}}{Y_ik} is the intensity of a continuous random
##' variable for sample \eqn{i}{i} in condition \eqn{k}{k} and
##' \eqn{\beta_k}{beta_k} is the mean parameter of condition
##' \eqn{k}{k}. We denote by \eqn{K}{K} the total number of conditions
##' and \eqn{n_k}{n_k} the number of sample in each condition.
##'
##' More details related to the weights are coming...
##'
##' @seealso See also \code{\linkS4class{fusedanova}},
##' \code{\link{plot,fusedanova-method}} and \code{\link{cv.fa}}.
##' @name fusedanova
##' @rdname fusedanova
##' @keywords models, regression
##'
##' @examples \dontrun{
##' data(aves)
##' fa.laplace <- fusedanova2(x=aves$weight, group=aves$family, weighting="laplace", gamma=5)
##' plot(fa.laplace, labels=aves$order)
##'
##' fa.gaussian <- fusedanova(x=aves$weight, class=aves$family, weights="gaussian", gamma=5)
##' plot(fa.ttest, labels=aves$order)
##'
##' fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
##' plot(fa.ada, labels=aves$order)
##' }
##'
##' @export 
fusedanova2 <- function(x, group = 1:length(x),
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
  
  ## the vector of means is the 
  mean_k <- rowsum(x, group)/nk
  
  # data compression and ordering
  order  <- order(mean_k) 
  slopes <- get_slopes(mean_k[order], nk[order], weighting, gamma, W)
  out    <- fuse(mean_k[order], slopes, nk[order]) 

  hc <- structure(list(merge  = out$merge,
                  height = out$path$lambda, 
                  labels = levels(group)[order],
                  order = out$order), class = "hclust")

  res <- list(path = out$path, hc = hc, slopes = slopes)
  res
}

slopes <- function(x, group, gamma = 1) {
  
  nk <- tabulate(group)  
  k <- length(nk)
  mean_k <- rowsum(x, group)/nk
  order <- order(mean_k)

  nk <- nk[order]
  mean_k <- mean_k[order]
  
  ## as fast à C++
  ## Laplace weights (nk.nl exp(- gamma | yk - yl|)), computation in O(n)/O(K)
  c1 <- rev(cumsum(c(0,rev(nk * exp(-gamma*mean_k))[-k])))
  c2 <- cumsum(c(0,(nk * exp(gamma*mean_k))[-k]))
  w <- exp(gamma*mean_k) * c1 - exp(-gamma*mean_k) * c2
  w
}

# slopes0 <- function(x, group, gamma = 1) {
#   
#   
#   nk <- tabulate(group)  
#   k <- length(nk)
#   mean_k <- rowsum(x, group)/nk
#   order <- order(mean_k)
# 
#   nk <- nk[order]
#   mean_k <- mean_k[order]
#   
#   ## as fast à C++
#   ## Laplace weights (nk.nl exp(- gamma | yk - yl|)), computation in O(n)/O(K)
#   s <- numeric(k)
#   w <- matrix(0,k,k)
#   for (i in 1:k) {
#     for (j in 1:k) {
#       w[i,j] <- nk[j] * exp(-gamma * abs(mean_k[i] - mean_k[j]) )
#     }
#     s[i] <- - sum( w[i, ] * sign (i - 1:k))
#   }
#   s
# }
