##' @rdname fusedanova
##' @export 
fusedanova.numeric <- function(x, 
                        group,
                        weighting = c("laplace", "gaussian", "adaptive"),
                        gamma = 0, standardize = TRUE, W = NULL) {

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
    x <- (x - mean(x)) / s
  }
  
  # data compression and ordering
  mean_k   <- rowsum(x, group) / nk
  ordering <- order(mean_k)
  
  ## call to fused-ANOVA
  slopes <- get_slopes(mean_k[ordering], nk[ordering], gamma, weighting, W)
  fa_out <- fusedanova_cpp(mean_k[ordering], slopes, nk[ordering]) 

  ## creating the hc object
  hc <- structure(
    list(
      merge  = fa_out$merge,
      height = fa_out$path$lambda, 
      labels = levels(group)[ordering],
      order  = fa_out$order # order for plotting the dendrogram
    ), class = "hclust")

  ## creating the fused-anova object
  fa_object <- structure(
    list(
      x_bar   = mean_k,
      order   = ordering, 
      path    = fa_out$path, 
      hc      = hc,
      call    = match.call
    ),
    class = "fusedanova")
  
  fa_object
}
