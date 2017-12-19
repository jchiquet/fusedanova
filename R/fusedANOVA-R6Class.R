##' An R6 Class for representing and performing optimization of a fused-ANOVA model
##'
##' @description
##'
##' @importFrom R6 R6Class
##' @export
##' 
fusedANOVA <-
   R6Class(classname = "fusedANOVA",
    public = list(
      data       = NULL,
    	cl0        = NULL,
	    weighting  = NULL,
      hclust     = NULL,
	  	path       = NULL,
	  	merge      = NULL,
	  	order      = NULL,
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
      nk     = NULL  # group sizes in class0
    ), 
    active = list(
      penalties = function() {self$path$lambda}
    )
  )

# fusedANOVA$set("public", "get_slopes2", 
#   function() {
#     mean_k <- rowsum(self$data, self$cl0)/private$nk
#     var_k  <- rep(0,private$k)
#     if (self$weighting %in% c("welch", "naivettest", "ttest"))
#       var_k <- private$nk/(private$nk - 1)*(rowsum(self$data^2,self$cl0)/private$nk - mean_k^2)
#     self$order <- order(mean_k)
# 
#     nk <- private$nk[self$order]    
#     mean_k <- mean_k[self$order]    
#     var_k <- var_k[self$order]    
#     
#     ## as fast à C++
#     if (self$weighting == "laplace") {
#       ## Laplace weights (nk.nl exp(- gamma | yk - yl|)), computation in O(n)/O(K)
#       c1 <- rev(cumsum(c(0,rev(nk * exp(-mean_k))[-private$k])))
#       c2 <- cumsum(c(0,(nk * exp(mean_k))[-private$k]))
#       w <- exp(mean_k) * c1 - exp(-mean_k) * c2
#     } else {
#       ## default weights (nk.nl), computation in O(n)/O(K)
#       w <- sum(nk) + nk - 2 * cumsum(nk)
#     }
#     w
#   }
# )
