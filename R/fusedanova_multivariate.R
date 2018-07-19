##' @rdname fusedanova
##' @export 
fusedanova.data.frame <- 
  function(x, group,
           weighting = c("laplace", "gaussian", "adaptive"),
           gamma = rep(0,ncol(x)), standardize = FALSE, W = NULL) {
    
    ## group vector: default or/and conversion to a factor
    if (missing(group)) {
      if (is.null(rownames(x))) {
        group <- factor(paste0("ind",1:nrow(x)))
      } else {
        group <- factor(rownames(x))        
      }
    } else if (!is.factor(group)) group <- as.factor(group)
    
    
    ## overwrite default parameters with user's
    weighting <- match.arg(weighting)

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
      s <- sapply(x, get_norm, group, n, k, nk)
      m <- colMeans(x)
      x <- sweep(x, 2, m, "-")
      x <- sweep(x, 2, s, "/")
    }

    ## MULTIPLE CALLS TO UNIVAIRATE FUSED-ANOVA
    fa_objs <- mapply(
      fusedanova.numeric, 
      as.list(x), gamma, 
      MoreArgs = list(group = group, weighting = weighting, standardize = FALSE, W = W),
      SIMPLIFY = FALSE
    )
    
    ## list of dendrograms (hclust)
    hc_objs <- lapply(fa_objs, as.hclust.fusedanova)

    ### AGGREGATION
    ## extract lists of rules and lambdas
    Rules  <- mapply(function(fa, hc) {
        list(
          rules = as.matrix(subset(fa$path, select = c(down, split, up)))[(k - 1):1, ],
          order = hc$order
        )
    }, fa_objs, hc_objs, SIMPLIFY = FALSE)

    Lambdas  <- lapply(fa_objs, function(fa) rev(fa$path$lambda))
    
    o_lambda <- order(unlist(Lambdas), decreasing = TRUE)
    orderRules <- as.matrix(cbind(rep(1:(k - 1), p),rep(1:p, each = k - 1))[o_lambda,])
    
    ## build the aggregating rules 
    aggregation <- pruneSplits(Rules, orderRules, k, p)

    ## height in the tree are obatined from the lambdas
    heights <- sapply(rev(aggregation$rule), FUN = function(rule){
      Lambdas[[orderRules[rule,2]]][orderRules[rule,1]]
    })

### TODO: create fa object (path)
### Is it possible from aggregation?
    # res <- structure(
    #     list(
    #       labels  = fa_out[[1]]$labels, 
    #       path  = NULL),
    #     class = "fusedanova")
    # res
        
    ## creating hclust object
    hc <- structure(list(
      merge  = CreationMatriceMerge(subset(aggregation, select = 1:4)),
      height = heights[-length(heights)], 
      labels = fa_objs[[1]]$labels,
      order  = OrdreIndividus(aggregation)), class = "hclust")
    hc
  }
