##' Fit a Fused ANOVA model
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total Variation) penality, 
##' ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the difference of each coefficient. 
##' 
##' @param x vector of observation for n individuals.
##'
##' @param group vector or factor giving the initial group of each individual. If missing, 
##' \code{1:length(x)} is used (clustering mode with one individual per group).
##'
##' @param weighting character; which type of weights is supposed to be used.
##' The supported weights are: \code{"laplace"}, \code{"gaussian"} or  \code{"adaptive"}.
##' See details below. By default, its value is \code{"laplace"}. Ignore if \code{W} is not \code{NULL}.
##' 
##' @param standardize logical; should the vector be standardized before computation?
##' Default is \code{FALSE}. If \code{TRUE}, the vector is centered and scaled according to the pooled variance of initial groups.
##' 
##' @param gamma non-negative scalar; the \eqn{\gamma}{gamma} parameter needed for
##' \code{"laplace"}, \code{"gaussian"} and \code{"adaptive"} weights. Default is 0.
##'
##' @param W  a numeric matrix of weights of user defined weights. Default is \code{NULL}. 
##' If not \code{NULL}, should be a k x k matrix (with k the initial number of groups) that
##' will overwrite the \code{weighting} parameter.
##' 
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
##' @name fusedanova
##' @rdname fusedanova
##' @keywords models, regression
##'
##' @examples \dontrun{
##' data(aves)
##' fa.laplace <- fusedanova2(aves$weight, aves$family, gamma=0)
##' plot(fa.laplace)
##'
##' fa.ada <- fusedanova(aves$weight, aves$family, weights="adaptive", gamma=2)
##' plot(fa.ada)
##' }
##'
##' @export 
fusedanova <- function(x, group = 1:length(x),
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
  out    <- fusedanova_cpp(mean_k[order], slopes, nk[order]) 

  hc <- structure(list(merge  = out$merge,
                  height = out$path$lambda, 
                  labels = levels(group)[order],
                  order = out$order), class = "hclust")

  res <- structure(list(x_bar = mean_k, group = group, lambda = out$pathlambda,
                        path = out$path, hc = hc, call = match.call), class = "fusedanova")
  res
}

#' plot a fusedanova object
#' 
#' plot a fusedanova object
#' 
#' @export
#' 
plot.fusedanova <- function(x, type = c("dendrogram", "BIC"), ...) {
  stopifnot(inherits(x, "fusedanova"))
  type <- match.arg(type)
  if (type == 'dendrogram')
    plot(x$hc, ...)
  if (type == 'BIC')
    cat("\nNot yet implemented")
}

#' compute loglikelihood of a fusedanova object
#' 
#' compute loglikelihood of a fusedanova object
#' 
#' @export
#' 
logLik.fusedanova <- function(object, ngroups = 1:(nlevels(object$group) - 1), heights = NULL) {
  groups <- cutree(object$hc, k = ngroups, h = heights)
  loglik <- apply(groups, 2, loglik.ANOVA, object$x)
  loglik
}

#' AIC of a fusedanova object
#' 
#' compute AIC of a fusedanova object
#' 
#' @export
#' 
AIC.fusedanova <- function(object, ngroups = 1:nrow(object$path), heights = NULL, k=2) {
  groups <- cutree(object$hc, k = ngroups, h = heights)
  loglik <- apply(groups, 2, loglik.ANOVA, object$x)
  df <- apply(groups, 2, function(grp) length(unique(grp)))
  AIC <- -2 * loglik + k*df 
  AIC
}

#' BIC of a fusedanova object
#' 
#' compute BIC of a fusedanova object
#' 
#' @export
#' 
BIC.fusedanova <- function(object, ngroups = 1:nrow(object$path), heights = NULL) {
  BIC <- AIC.fusedanova(object, heights = heights, ngroups = ngroups, k = log(length(object$x)))  
  BIC
}

loglik.ANOVA <- function(group, x){
  n  <- length(x)
  nk <- tabulate(group)
  k <- length(nk)
  betak <- rowsum(x, group)/nk
  RSS <- sum((x - betak[group])^2)
  sigma2 <- RSS/(n-k)
  loglik <- -.5 * (n*log(2*pi) + n*sum(log(sigma2)) + (n-k))
  loglik
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

