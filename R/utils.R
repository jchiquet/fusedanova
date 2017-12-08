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

# default.args.cv <- function() {
#   return(list(
#     weights = "default",
#     W=matrix(nrow=0, ncol=0),
#     gamma = 1,
#     standardize = TRUE,
#     nlambda =100,
#     log.scale = TRUE,
#     min.ratio = 1e-8,
#     verbose = FALSE
#   ))
# }

stopif <- function(expr, message) {
  if (expr) stop(message)  
}
