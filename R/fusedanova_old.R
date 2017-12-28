##' Fit a Fused ANOVA model with old version of the package
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total Variation) penality, 
##' ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the difference of each coefficient. 
##' 
##' @export
fusedanova_old <- function(x, class = 1:length(x),
                           weights = c("default", "laplace", "gaussian", "adaptive", "personal"),
                           standardize = FALSE, ...) {
  
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
  
  return(list(result = list(res),
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

  o <- order(xm)
  xm <- xm[o] # sort from the smallest beta to the highest
  ngroup <- ngroup[o]

  slopes <- get_slopes(xm, ngroup, args$weights, args$gamma, args$W)  
  res  <- fusedanova_cpp_old(xm, slopes, ngroup)

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

