##' An R6 Class for representing and performing optimization of a fused-ANOVA model
##'
##' @description
##'
##' @importFrom R6 R6Class
##' @importFrom dplyr arrange desc filter %>%
##' @export
##' 
fusedANOVA <-
   R6Class(classname = "fusedANOVA",
    public = list(
      data      = NULL,
    	cl0       = NULL,
	    weighting = NULL,
	  	path      = NULL,
	  	order     = NULL,
	  	initialize = function(data, class0, weighting, standardize) {
        self$data      <- data
        private$n      <- length(data)
        self$weighting <- weighting
        self$cl0       <- class0
        private$nk <- tabulate(class0)
        private$k  <- length(private$nk)
        if (standardize) private$standardize()
	  	}
    ),
    private = list(
      n      = NULL, # sample size
      k      = NULL, # number of groups in class0
      nk     = NULL, # group sizes in class0
      fusion = NULL  # path with fusion only
    ), 
    active = list(
      penalties = function() {self$path$lambda},
      fusionTable = function() {private$fusion}
      
    )
  )

fusedANOVA$set("private", "standardize",
  function() {
     if (private$k != private$n) {
       # if any initial grouping, normalize withing each group
       s <- (rowsum(self$data^2,self$cl0) - (1/private$nk) * (rowsum(self$data,self$cl0))^2) / (private$nk - 1)
       s[private$nk == 1] <- 0
       s <- sqrt(sum(s*(private$nk - 1))/(private$n - private$k))
     } else {
       # if no grouping (one guy per class0) normalize at the vector scale 
       s <- sd(self$data)
     }
     self$data <- (self$data - mean(self$data))/s
  }
)

fusedANOVA$set("public", "get_path", 
  function(args) {
    mean_k <- rowsum(self$data, self$cl0)/private$nk
    var_k  <- rep(0,private$k)
    if (self$weighting %in% c("welch", "naivettest", "ttest"))
      var_k <- private$nk/(private$nk - 1)*(rowsum(self$data^2,self$cl0)/private$nk - mean_k^2)
    self$order <- order(mean_k)

    slopes <- get_slopes(mean_k[self$order], private$nk[self$order], var_k[self$order], self$weighting, 1, matrix(0,0,0) )
    
    self$path  <- fuse(mean_k[self$order], slopes, private$nk[self$order])

    private$fusion <- self$path %>% filter(idown != iup) %>% arrange(desc(lambda))
      
    invisible(self)
  }
)

fusedANOVA$set("public", "cut_tree", 
  function(heights = NULL) {
    if (is.null(heights)) {
      heights <- self$penalties
    } else {
      stopifnot(all(heights %in% fusion$lambda))
    }
    heights <- sort(unique(heights), decreasing = TRUE)        

    cl <- get_clustering(heights, private$fusion$lambda, private$fusion$idown, private$fusion$iup, private$k)
    cl <- cl[self$order, ]
    list(cl = cl, heights = heights)
  }
)

fusedANOVA$set("public", "get_slopes", 
  function() {
    mean_k <- rowsum(self$data, self$cl0)/private$nk
    var_k  <- rep(0,private$k)
    if (self$weighting %in% c("welch", "naivettest", "ttest"))
      var_k <- private$nk/(private$nk - 1)*(rowsum(self$data^2,self$cl0)/private$nk - mean_k^2)
    self$order <- order(mean_k)
    
    w <- get_slopes(mean_k[self$order], private$nk[self$order], var_k[self$order], self$weighting, 1, matrix(0,0,0) )  
    w
  }
)

fusedANOVA$set("public", "get_slopes2", 
  function() {
    mean_k <- rowsum(self$data, self$cl0)/private$nk
    var_k  <- rep(0,private$k)
    if (self$weighting %in% c("welch", "naivettest", "ttest"))
      var_k <- private$nk/(private$nk - 1)*(rowsum(self$data^2,self$cl0)/private$nk - mean_k^2)
    self$order <- order(mean_k)

    nk <- private$nk[self$order]    
    mean_k <- mean_k[self$order]    
    var_k <- var_k[self$order]    
    
    ## as fast à C++
    if (self$weighting == "laplace") {
      ## Laplace weights (nk.nl exp(- gamma | yk - yl|)), computation in O(n)/O(K)
      c1 <- rev(cumsum(c(0,rev(nk * exp(-mean_k))[-private$k])))
      c2 <- cumsum(c(0,(nk * exp(mean_k))[-private$k]))
      w <- exp(mean_k) * c1 - exp(-mean_k) * c2
    } else {
      ## default weights (nk.nl), computation in O(n)/O(K)
      w <- sum(nk) + nk - 2 * cumsum(nk)
    }
    w
  }
)

