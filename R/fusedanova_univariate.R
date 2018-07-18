##' @rdname fusedanova
##' @export 
fusedanova.numeric <- function(x, 
                        group,
                        weighting = c("laplace", "gaussian", "adaptive"),
                        gamma = 0, standardize = FALSE, W = NULL) {

  ## group vector: default or/and conversion to a factor
  if (missing(group)) {
    if (is.null(names(x))) {
      group <- factor(paste0("ind",1:length(x)))
    } else {
      group <- factor(names(x))        
    }
  } else if (!is.factor(group)) group <- as.factor(group)

  ## overwrite default parameters with user's
  weighting <- match.arg(weighting)
  if (!is.null(W)) weighting <- "personal" else W <- matrix(0, 0, 0)
  
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
  if (weighting == "personal") stopif(nrow(W) != n | ncol(W) != n, "W must be a square matrix.")
  
  # data standardization
  if (standardize) {
    s <- get_norm(x, group, n, k, nk)
    m <- mean(x)
    x <- (x - m) / s
  }
  
  # data compression and ordering
  group_means <- rowsum(x, group) / nk
  ordering    <- order(group_means)
  group_names <- levels(group)[ordering]
  group_means <- group_means[ordering]
  group_sizes <- nk[ordering]
  names(group_means) <- group_names

  ## call to fused-ANOVA
  slopes <- get_slopes(group_means, group_sizes, gamma, weighting, W)
  fa_out <- fusedanova_cpp(group_means, slopes, group_sizes) 

  ## creating the hc object
  hc <- structure(
    list(
      merge  = fa_out$merge,
      height = fa_out$path$lambda, 
      labels = group_names,
      order  = fa_out$order # order for plotting the dendrogram
    ), class = "hclust")

  ## creating the fused-anova object
  fa_object <- structure(
    list(
      means = group_means,
      path  = fa_out$path, 
      hc    = hc,
      call  = match.call
    ),
    class = "fusedanova")
  
  fa_object
}
