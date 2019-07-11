#' Fit a fast version of the weighted l1-clusterpath algorithm
#' 
#' The l1-clusterpath algorithm is a convex clustering algorithm with fused-LASSO 
#' (or Total Variation) penality, ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the 
#' difference of each coefficient. 
#' 
#' @param x a numeric vector of observation for n individuals.
#'
#' @param group an optional vector or factor giving the initial grouping. If missing, 
#' each individual are set in a single group.
#'
#' @param weighting character; which type of weights is supposed to be used.
#' The supported weights are: \code{"laplace"}, \code{"gaussian"} or  \code{"adaptive"}.
#' See details below. \code{"laplace"} by default.
#' 
#' @param gamma non-negative scalar ; the \eqn{\gamma}{gamma} parameter is needed for
#' \code{"laplace"}, \code{"gaussian"} and \code{"adaptive"} weights. Default is 0.
#'
#' @param hclust boolean: should the result be outputed as an hclust object? Default is \code{TRUE}.
#' 
#' @return an S3 object with class \code{hclust} or a data frame of the succesive fusions.
#'
#' The optimization problem solved is
#' \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda} = \arg \min_{\beta}
#' \left\{\sum_{k=1}^K \sum_{i=1}^{n_k} \left(Y_{ik}-\beta_k \right)^2
#' + \lambda \sum_{k,\ell} w_{kl} \left|\beta_k - \beta_\ell \right|\right\}}}
#' \if{html}{\out{ <center> &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub></sub> =
#' argmin<sub>&beta;</sub> sum<sub>k</sub> sum_i (Y<sub>ik</sub> - &beta<sub>k</sub>)<sup>2</sup>
#' + &lambda; sum<sub>k,l</sub> w<sub>k,l</sub>
#' &#124; &beta;<sub>k</sub> - &beta;<sub>l</sub> &#124;, </center> }}
#' \if{text}{\deqn{beta.hat(lambda) = argmin_beta sum_k sum_i (Y_ik - beta_k)^2
#' + lambda sum_k sum_l w_kl | beta_k - beta_l|,}}
#'
#' where \eqn{Y_{ik}}{Y_ik} is the intensity of a continuous random
#' variable for sample \eqn{i}{i} in condition \eqn{k}{k} and
#' \eqn{\beta_k}{beta_k} is the mean parameter of condition
#' \eqn{k}{k}. We denote by \eqn{K}{K} the total number of conditions
#' and \eqn{n_k}{n_k} the number of sample in each condition.
#'
#' @name clusterpath_l1
#' 
#' @references
#' Chiquet J, Gutierrez P, Rigaill G: Fast tree inference with weighted fusion penalties,
#'  Journal of Computational and Graphical Statistics 205–216, 2017.
#'  
#' T. Hocking, J.-P. Vert, F. Bach, and A. Joulin. Clusterpath: an
#' Algorithm for Clustering using Convex Fusion Penalties, ICML,
#' 2011.
#'
#' @examples
#' data(aves)
#' clpath <- clusterpath_l1(aves$weight, aves$family, gamma = 0)
#' plot(clpath)
#'
#' @export
clusterpath_l1 <- function(x, group,
                           weighting = c("laplace", "gaussian", "adaptive"),
                           gamma = 0, hclust = TRUE) {

  ## group vector: default or/and conversion to a factor
  if (missing(group)) {
    if (is.null(names(x))) {
      group <- factor(paste0(1:length(x)))
    } else {
      group <- factor(names(x))        
    }
  } else if (!is.factor(group)) group <- as.factor(group)

  ## overwrite default parameters with user's
  weighting <- match.arg(weighting)

  ## problem dimensions
  n  <- length(x)
  k  <- length(unique(group))
  nk <- tabulate(group)
  
  ## check to partially avoid crashes of the C++ code
  stopif(!is.numeric(x)     , "x must be a numeric vector.")
  stopif(gamma < 0          , "gamma must be non-negative.")
  stopif(anyNA(x)           , "NA value in x not allowed.")
  stopif(n != length(group) , "x and group length do not match")
  stopif(k == 1             , "x only has one level: there's no point in fusing one group, you know...")

  # data compression and ordering
  group_means <- rowsum(x, group) / nk
  ordering    <- order(group_means)
  group_names <- levels(group)[ordering]
  group_means <- group_means[ordering]
  group_sizes <- nk[ordering]

  ## call to fused-ANOVA cpp
  slopes <- get_slopes(group_means, group_sizes, gamma, weighting)
  fusion <- fusedanova_cpp(group_means, slopes, group_sizes)

  ## creation of the univarclust object
  res <- structure(
    list(fusionTree = fusion,
         weighting  = weighting,
         labels     = group_names, 
         call       = match.call()), class = "univarclust")
  if (hclust) res <- as_hclust(res, "l1-clusterpath")
  res
}
