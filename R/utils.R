stopif <- function(expr, message) {
  if (expr) stop(message)  
}

# export to hclust format
#
# export a fusedanova_cpp output to an hclust object - internal used only
as_hclust <- function(object, method) {
  merge <- export_merge(object$fusionTree$child1, object$fusionTree$child2)
  mergeReordered <- merge
  ## ind <- match(1:length(object$labels), order(as.numeric(as.factor(object$labels))))
  ind <- 1:length(object$labels)
  mergeReordered[merge < 0] <- -ind[-merge[merge < 0]]
  order <- export_order(mergeReordered, object$fusionTree$sizes)
  hc <- structure(
    list(
      merge    = mergeReordered,
      height   = object$fusionTree$lambda, 
      labels   = object$labels,
      order    = order,
      method   = method,
      dist.method = paste(object$weighting, "weights"),
      call = object$call
    ), class = "hclust")
  hc
}
