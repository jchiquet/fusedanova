fusedANOVA_args <- function(weights, standardize, user_args) {
  args <- modifyList(
    list(
      W           = matrix(nrow = 0, ncol = 0),
      weights     = weights,
      gamma       = 1,
      verbose     = FALSE
    ), user_args)
  args
}

stopif <- function(expr, message) {
  if (expr) stop(message)  
}
