#' Fit a fast version of the 1 dimensional Ward clustering algorithm
#' 
#' The Ward clustering algorithm is the more statistically grounded hierarchical clustering approach.
#' This version is dedicated to clustering of univariate data (i.e., vector).
#' 
#' @param x a numeric vector of observation for n individuals.
#'
#' @param hclust boolean: should the result be outputed as an hclust object? Default is \code{TRUE}.
#' 
#' @return an S3 object with class \code{hclust} or a data frame of the succesive fusions.
#'
#' @examples
#' data(aves)
#' ward1d <- ward_1d(aves$weight)
#' plot(ward1d)
#'
#' @export 
ward_1d <- function(x, hclust = TRUE) {

  ## check to partially avoid crashes of the C++ code
  stopif(!is.numeric(x)     , "x must be a numeric vector.")
  stopif(anyNA(x)           , "NA value in x not allowed.")
  
  n <- length(x)
  if (is.null(names(x))) names(x) <- 1:n
  
  # data ordering
  ordering <- order(x)
  sums     <- x[ordering]
  sum2s    <- (x^2)[ordering]

  ## call to C++
  fusion <- ward1d_cpp(sums, sum2s, rep(1, n))

  ## creation of the fused-anova object
  res <- structure(
    list(fusionTree = fusion,
         weighting  = "none",
         method     = "Ward 1D",
         labels     = names(x),
         ordering   = ordering,
         call       = match.call()), class = "univarclust")
  if (hclust) res <- as_hclust(res)
  res
}
