##' Fit a Fused ANOVA model
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total Variation) penality, 
##' ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the difference of each coefficient. 
##' 
##' @param x a vector, matrix or data.frame of observation for n individuals.
##'
##' @param group vector or factor giving the initial group of each individual. If missing, 
##' each individual are set in a single group is used (clustering mode).
##'
##' @param weighting character; which type of weights is supposed to be used.
##' The supported weights are: \code{"laplace"}, \code{"gaussian"} or  \code{"adaptive"}.
##' See details below. By default, its value is \code{"laplace"}. Ignore if \code{W} is not \code{NULL}
##' 
##' @param standardize logical; should the vector be standardized before computation?
##' Default is \code{FALSE}. If \code{TRUE}, the vector is centered and scaled according to the pooled variance of initial groups.
##' 
##' @param gamma non-negative scalar (or vector if x is a matrix or a data frame);
##'  the \eqn{\gamma}{gamma} parameter needed for
##' \code{"laplace"}, \code{"gaussian"} and \code{"adaptive"} weights. Default is 0.
##'
##' @param W  a numeric matrix of weights of user defined weights. Default is \code{NULL}. 
##' If not \code{NULL}, should be a k x k matrix (with k the initial number of groups) that
##' will overwrite the \code{weighting} parameter.
##' 
##' @return an S3 object with class \code{fusedanova}.
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
##' @name fusedanova
##' @rdname fusedanova
##' @keywords models, regression
##'
##' @examples
##' data(aves)
##' fa.laplace <- fusedanova(aves$weight, aves$family, gamma = 0)
##' plot(fa.laplace)
##'
##' fa.ada <- fusedanova(aves$weight, aves$family, "adaptive", gamma = 2)
##' plot(fa.ada)
##'
##' @include fusedanova_univariate.R fusedanova_multivariate.R
##' @export 
fusedanova <- function(x, ...) UseMethod("fusedanova", x)

##' @rdname fusedanova
##' @export 
fusedanova.matrix <- 
  function(x, group = 1:nrow(x),
           weighting = c("laplace", "gaussian", "adaptive"),
           gamma = rep(0,ncol(x)), standardize = TRUE, W = NULL) {
  res <- fusedanova.data.frame(
    as.data.frame(x), 
    group = group, 
    weighting = weighting, 
    gamma = gamma, 
    standardize = standardize, 
    W = W
  )
  res
}

#' plot a fusedanova object
#' 
#' plot a fusedanova object
#' 
#' @export
#' 
plot.fusedanova <- function(x, ...) {
  stopifnot(inherits(x, "fusedanova"))
  plot(as.hclust.fusedanova(x), ...)
}

#' export to hclust format
#'
#' export a fusedanova pobject to an hclust object
#'
#' @export
#'
as.hclust.fusedanova <- function(object, ...) {
  merge <- export_merge(object$path$parent1, object$path$parent2)
  order <- export_order(merge, object$path$sizes)
  hc <- structure(
    list(
      merge  = merge,
      height = object$path$lambda, 
      labels = object$labels,
      order  = order
    ), class = "hclust")
  hc
}

#' #' compute loglikelihood of a fusedanova object
#' #' 
#' #' compute loglikelihood of a fusedanova object
#' #' 
#' #' @export
#' #' 
#' logLik.fusedanova <- function(object, ngroups=1:nrow(object$path)) {
#'   groups <- cutree(object$hc, k = ngroups)
#'   loglik <- apply(groups, 2, loglik_ANOVA, object$means)
#'   loglik
#' }
#' 
#' #' AIC of a fusedanova object
#' #' 
#' #' compute AIC of a fusedanova object
#' #' 
#' #' @export
#' #' 
#' AIC.fusedanova <- function(object, ngroups = 1:nrow(object$path), k = 2) {
#'   loglik <- logLik.fusedanova(object, ngroups)
#'   AIC <- -2 * loglik + k * ngroups
#'   AIC
#' }
#' 
#' #' BIC of a fusedanova object
#' #' 
#' #' compute BIC of a fusedanova object
#' #' 
#' #' @export
#' #' 
#' BIC.fusedanova <- function(object, ngroups = 1:nrow(object$path)) {
#'   BIC <- AIC.fusedanova(object, ngroups = ngroups, k = log(length(object$means)))  
#'   BIC
#' }

# slopes <- function(x, group, gamma = 1) {
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
#   c1 <- rev(cumsum(c(0,rev(nk * exp(-gamma*mean_k))[-k])))
#   c2 <- cumsum(c(0,(nk * exp(gamma*mean_k))[-k]))
#   w <- exp(gamma*mean_k) * c1 - exp(-gamma*mean_k) * c2
#   w
# }

