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
	  	result    = NULL,
	    penalties = NULL,
		  algorithm = NULL,
  		
	  	initialize = function(data, class0, weighting, standardize) {
        self$data      <- data
        private$n      <- nrow(data)
        private$p      <- ncol(data)
        self$weighting <- weighting
        self$cl0       <- class0
        private$nk <- tabulate(class0)
        private$k  <- length(private$nk)
  		  if (standardize) private$standardize()
	  	}
	  	
    ),
    private = list(
      n  = NULL,
      p  = NULL,
      k  = NULL,
      nk = NULL,
      standardize = function() {
        self$data <- lapply(self$data, function(x) {
          ## Expect in the clustering case ... (one guy per class0)
          if (private$k == private$n) {
            s <- sd(x)
          } else {
            ## Use Pooled variance
            n <- length(x)
            s <- (rowsum(x^2,self$cl0) - (1/private$nk) * (rowsum(x,self$cl0))^2) / (private$nk - 1)
            s[self$cl0_size == 1] <- 0
            s <- sqrt(sum(s*(private$nk - 1))/(private$n - private$k))
          }
          res <- (x - mean(x))/s
          res
          })
        },
      homotopy = function() {
        
      }
    )
  )

