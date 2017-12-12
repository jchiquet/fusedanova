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
##' @parma weights character; which type of weights is supposed to be used.
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
##' fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
##' plot(fa.laplace, labels=aves$order)
##'
##' fa.ttest <- fusedanova(x=aves$weight, class=aves$family, weights="naivettest")
##' plot(fa.ttest, labels=aves$order)
##'
##' fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
##' plot(fa.ada, labels=aves$order)
##' }
##'
##' @export
fusedanova <- function(x, class = 1:length(x),
                       weights = c("default", "laplace", "gaussian", "adaptive", "naivettest", "ttest", "welch", "personal"),
                       standardize = TRUE, ...) {
  
  ## overwrite default parametrs with user's
  weights <- match.arg(weights)
  args <- fusedANOVA_args(weights, standardize, list(...))

  ## check to partilly avoid crashes of th C++ code
  stopif(!is.numeric(x)                 , "x must be a numeric vector")
  stopif(any(is.na(x))                  , "NA value in x not allowed.")
  stopif(length(x) != length(class)     , "data and class dimensions do not match")
  stopif(length(unique(class)) == 1     , "y has only one level.")
  
  # conversion of class ot a factor
  if (!is.factor(class)) class <- as.factor(class)
  
  # normalization 
  if (standardize) x <- normalize(x, class)

  res <- calculatepath(x, class, args)

  res <- list(table = res$table, order = res$order)
  
  ## small warning on last beta
  if (standardize == TRUE && abs(res$table[1,1]) > 10^(-8))
    warning("There may be some approximation errors (Beta(lambda_max) far from 0). You may want to lower the gamma if you are using one.")
  
  return(new("fusedanova",
             result = list(res),
             classes = class,
             weights = weights,
             algorithm = "No Split"))
  
}

#########################################
# Calculate path with or without splits #
#########################################
calculatepath <- function(x, group, args) {

  ngroup <- tabulate(group)# vector of number by group
  xm <- rowsum(x,group)/ngroup
  xv <- rep(0,length(xm))
  if (args$weights %in% c("welch", "naivettest", "ttest")){
    ## var needed if weights are of welch or ttest type
    xv <- ngroup/(ngroup-1)*(rowsum(x^2,group)/ngroup - xm^2)
  }

  o <- order(xm)
  xm <- xm[o] # sort from the smallest beta to the highest
  ngroup <- ngroup[o]
  xv <- xv[o]

  slopes <- get_slopes(xm, ngroup, xv, args$weights, args$gamma, args$W)  
  res  <- fuse_old(xm, slopes, ngroup)

  return(list(table = res, order = o))
  
}

#############################
# normalization of a vector
#############################
normalize <- function(x,group){
  ## if Pooled 
  if (nlevels(group) == length(x)) {
    s <- sd(x)
  } else {
    n <- length(x)
    ngroup <- tabulate(group)
    s <- 1/(ngroup-1)*(rowsum(x^2,group) - (1/ngroup)*(rowsum(x,group))^2)
    s[which(ngroup == 1)] <- 0
    if (sum(s)==0){
      s <- sd(x)
    } else{
      s <- sqrt(sum(s*(ngroup-1))/(n-length(ngroup)))
    }
  }
  res <- (x - mean(x))/s
  res
}

##' @export 
fusedanova2 <- function(x, class = rep(1:length(x)),
                       weights = c("default", "laplace", "gaussian", "adaptive", "naivettest", "ttest", "welch", "personal"),
                       standardize = TRUE, ...) {
  
  ## overwrite default parametrs with user's
  weights <- match.arg(weights)
  args <- fusedANOVA_args(weights, standardize, list(...))
  
  ## check to partilly avoid crashes of th C++ code
  stopif(!is.numeric(x)                 , "x must be a numeric vector")
  stopif(any(is.na(x))                  , "NA value in x not allowed.")
  stopif(length(x) != length(class)     , "data and class dimensions do not match")
  stopif(length(unique(class)) == 1     , "y has only one level.")

  # conversion of class ot a factor
  if (!is.factor(class)) class <- as.factor(class)
  
  myFA <- fusedANOVA$new(data = x, class0 = class, weighting = weights, standardize = standardize)
  myFA$get_path(args)
  myFA
}

