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
	  	},
    get_path = function(args) {
      xm <- rowsum(self$data, self$cl0)/private$nk
      xv <- rep(0,private$k)
      if (self$weighting %in% c("welch", "ttest")) {
        ## var needed if weights are of welch or ttest type
        xv <- private$nk/(private$nk - 1)*(rowsum(self$data^2,self$cl0)/private$nk - xm^2)
      }
      
      self$order <- order(xm)
      xm <- xm[self$order] # sort from the smallest beta to the highest
      nk <- private$nk[self$order]
      xv <- xv[self$order]
      
      if (args$splits) {
        if (args$verbose) cat("\nPath calculated with possible splits")
        res  <- .Call("withSplit", R_x = xm, R_xv = xv, R_ngroup = nk, R_args = args, PACKAGE = "fusedanova")
      } else {
        if (args$verbose) cat("\nPath calculated without split")
        res  <- .Call("noSplit"  , R_x = xm,R_xv = xv,R_ngroup = nk, R_args = args, PACKAGE = "fusedanova")
      }
      
      self$path <- res$res
    }
    ),
    private = list(
      n  = NULL,
      k  = NULL,
      nk = NULL,
      standardize = function() {
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
  )

