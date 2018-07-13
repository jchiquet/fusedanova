##' @rdname fusedanova
##' @export 
fusedanova.data.frame <- 
  function(x, group = 1:nrow(x),
           weighting = c("laplace", "gaussian", "adaptive"),
           gamma = rep(0,ncol(x)), standardize = TRUE, W = NULL) {

    ## overwrite default parameters with user's
    weighting <- match.arg(weighting)
    if (!is.null(W)) weighting <- "personal" else W <- matrix(0,0,0)
    # conversion of group to a factor
    if (!is.factor(group)) group <- as.factor(group)
    
    ## problem dimensions
    n  <- nrow(x)
    p  <- ncol(x)
    k  <- length(unique(group))
    nk <- tabulate(group)
    
    ## check to partially avoid crashes of the C++ code
    stopif(!is.data.frame(x)  , "x must be a numeric vector.")
    stopif(length(gamma) != p , "gamma and ncol(x) do not match.")
    stopif(any(gamma < 0)     , "gamma must be non-negative.")
    stopif(anyNA(x)           , "NA value in x not allowed.")
    stopif(n != length(group) , "x and group length do not match")
    stopif(k == 1             , "x only has one level: there's no point in fusing one group, you know...")
    if (weighting == "personal") stopif(nrow(W) != n | ncol(W) != n, "W must be a square matrix.")
    
    # data standardization
    if (standardize) {
      m <- colMeans(x)
      s <- sapply(x, get_norm, group, n, k, nk)
      x <- sweep(x, 2, m, "-")
      x <- sweep(x, 2, s, "/")
    }
    
    ## data compression and ordering
    mean_k    <- sweep(rowsum(x, group), 1, nk, "/")
    orderings <- lapply(mean_k, order)
    mean_k    <- lapply(1:p, function(j) mean_k[orderings[[j]], j])
    nk        <- lapply(orderings, function(ordering) nk[ordering])
    
    ## call to fused-ANOVA
    slopes <- mapply(get_slopes, mean_k, nk, gamma, 
                     MoreArgs = list(weights = weighting, W = W), 
                     SIMPLIFY = FALSE)
    fa_out <- mapply(fusedanova_cpp, mean_k, slopes, nk, SIMPLIFY = FALSE)
    
    ## extract lists of rules and lambdas
    Rules  <- lapply(fa_out, function(fa) 
      list(rules = as.matrix(subset(fa$path, select = c(down, split, up)))[(k - 1):1, ],
           order = fa$order)
    )
    Lambdas  <- lapply(fa_out, function(fa) rev(fa$path$lambda))
    o_lambda <- order(unlist(Lambdas), decreasing = TRUE)
    orderRules <- as.matrix(cbind(rep(1:(k - 1), p),rep(1:p, each = k - 1))[o_lambda,])
    
    ## build the aggregating rules 
    aggregation <- pruneSplits(Rules, orderRules, k, p)

    ## height in the tree are obatined from the lambdas
    heights <- sapply(rev(aggregation$rule), FUN = function(rule){
      Lambdas[[orderRules[rule,2]]][orderRules[rule,1]]
    })
    
    ## creating hclust object
    hc <- structure(list(
      merge  = CreationMatriceMerge(subset(aggregation, select = 1:4)),
      height = heights[-length(heights)], 
      labels = levels(group)[fa_out[[1]]$order],
      order  = OrdreIndividus(aggregation)), class = "hclust")
    
    res <- structure(
        list(
          order  = fa_out[[1]]$order, 
          path  = NULL,
          hc    = hc, 
          call  = match.call),
        class = "fusedanova")
    res
  }
