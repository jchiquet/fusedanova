#' Fit a fast version of the 1 dimensional Ward clustering algorithm
#' 
#' The Ward clustering algorithm is the more statistically grounded hierarchical clustering approach.
#' This version is dedicated to clustering of univariate data (i.e., vector).
#' 
#' @param x a numeric vector of observation for n individuals.
#'
#' @param group an optional vector or factor giving the initial grouping. If missing, 
#' each individual are set in a single group.
#'
#' @param hclust boolean: should the result be outputed as an hclust object? Default is \code{TRUE}.
#' 
#' @return an S3 object with class \code{hclust} or a data frame of the succesive fusions.
#'
#' @examples
#' data(aves)
#' ward1d <- ward_1d(aves$weight, aves$family)
#' plot(ward1d)
#'
#' @export 
ward_1d <- function(x, group, hclust = TRUE) {

  ## group vector: default or/and conversion to a factor
  if (missing(group)) {
    if (is.null(names(x))) {
      group <- factor(paste0(1:length(x)))
    } else {
      group <- factor(names(x))        
    }
  } else if (!is.factor(group)) group <- as.factor(group)

  ## problem dimensions
  n  <- length(x)
  k  <- length(unique(group))
  nk <- tabulate(group)
  
  ## check to partially avoid crashes of the C++ code
  stopif(!is.numeric(x)     , "x must be a numeric vector.")
  stopif(anyNA(x)           , "NA value in x not allowed.")
  stopif(n != length(group) , "x and group length do not match")
  stopif(k == 1             , "x only has one level: there's no point in fusing one group, you know...")

  # data compression and ordering
  group_sums  <- rowsum(x  , group)
  group_sum2s <- rowsum(x^2, group)
  group_means <- group_sums / nk
  ordering    <- order(group_means)
  group_names <- levels(group)[ordering]
  group_sums  <- group_sums[ordering]
  group_sum2s <- group_sum2s[ordering]
  group_means <- group_means[ordering]
  group_sizes <- nk[ordering]

  ## call to C++
  fusion <- ward1d_cpp(group_sums, group_sum2s, group_sizes)

  ## creation of the fused-anova object
  res <- structure(
    list(fusionTree = fusion,
         weighting  = "none",
         labels     = group_names, 
         call       = match.call()), class = "univarclust")
  if (hclust) res <- as_hclust(res, "Ward 1D")
  res
}
