stopif <- function(expr, message) {
  if (expr) stop(message)  
}

# export to hclust format
#
# export a fusedanova_cpp output to an hclust object - internal used only
as_hclust <- function(object) {
  
  merge <- export_merge(object$fusionTree$child1, object$fusionTree$child2)
  mergeReordered <- merge
  mergeReordered[merge < 0] <- -object$ordering[-merge[merge < 0]]
  dendo_order <- export_order(mergeReordered, object$fusionTree$sizes)
      
  hc <- structure(
    list(
      merge    = mergeReordered,
      height   = object$fusionTree$lambda, 
      labels   = object$labels,
      order    = dendo_order,
      method   = object$method,
      dist.method = paste(object$weighting, "weights"),
      call = object$call
    ), class = "hclust")
  hc
}
