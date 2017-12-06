##' An R6 Class for representing and performing optimization of a fused-ANOVA model
##'
##' @description
##'
##' @importFrom R6 R6Class
##' @import dplyr
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
      penalties = function() {self$path$lambda}
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
    
    if (args$splits) {
      if (args$verbose) cat("\nPath calculated with possible splits\n")
      self$path  <- .Call("withSplit",
                      R_x      = mean_k[self$order],
                      R_xv     = var_k[self$order],
                      R_ngroup = private$nk[self$order],
                      R_args = args, PACKAGE = "fusedanova")$res
    } else {
      if (args$verbose) cat("\nPath calculated without split\n")
      self$path  <- .Call("noSplit",
                      R_x      = mean_k[self$order],
                      R_xv     = var_k[self$order],
                      R_ngroup = private$nk[self$order],
                      R_args = args, PACKAGE = "fusedanova")$res
    }
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
  })

