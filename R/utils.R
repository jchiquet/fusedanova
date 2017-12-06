fusedANOVA_args <- function(weights, standardize, user_args) {
  args <- modifyList(
    list(
      W           = matrix(nrow = 0, ncol = 0),
      weights     = weights,
      gamma       = 1,
      lambdalist  = numeric(0),
      checkargs   = TRUE,
      verbose     = FALSE,
      mxSplitSize = 100,
      epsilon     = 1e-10
    ), user_args)
  args
}

default.args.cv <- function() {
  return(list(
    weights = "default",
    W=matrix(nrow=0, ncol=0),
    gamma = 1,
    standardize = TRUE,
    splits = 0,
    epsilon =10^-10,
    checkargs = TRUE,
    nlambda =100,
    log.scale = TRUE,
    min.ratio = 1e-8,
    mc.cores = detectCores(),
    verbose = FALSE,
    mxSplitSize = 100
  ))
}

stopif <- function(expr, message) {
  if (expr) stop(message)  
}
