#' @export
print.backbone <- function(x, ...) {
  
  #Call
  cat("Call -\n")
  cat(deparse(x$call, width.cutoff = 500L)))
  cat("\n\n")
  
  #Source
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source: matrix-class bipartite network")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source: Matrix-class bipartite network")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source: igraph-class bipartite network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source: matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source: matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source: igraph-class weighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source: matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source: matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source: igraph-class unweighted network")}
  cat(source,"\n")
  
  #Result
  if (methods::is(x$backbone, "matrix")) {result <- paste0("Source: matrix-class unweighted network")}
  if (methods::is(x$backbone, "Matrix")) {result <- paste0("Source: Matrix-class unweighted network")}
  if (methods::is(x$backbone, "igraph")) {result <- paste0("Source: igraph-class unweighted network")}
  cat(source,"\n")
  
  #Model
  if (x$model == "sdsm") {model <- paste0("Backbone Model: Stochastic Degree Sequence Model")}
  if (x$model == "fdsm") {model <- paste0("Backbone Model: Fixed Degree Sequence Model")}
  if (x$model == "fixedrow") {model <- paste0("Backbone Model: Fixed Row Model")}
  if (x$model == "fixedcol") {model <- paste0("Backbone Model: Fixed Column Model")}
  if (x$model == "fixedfill") {model <- paste0("Backbone Model: Fixed Fill Model")}
  if (x$model == "disparity") {model <- paste0("Backbone Model: Disparity Filter")}
  if (x$model == "lans") {model <- paste0("Backbone Model: Locally Adaptive Network Sparsification")}
  if (x$model == "mlf") {model <- paste0("Backbone Model: Marginal Likelihood Filter")}
  if (x$model == "global") {model <- paste0("Backbone Model: Global Threshold")}
  if (x$model == "skeleton") {model <- paste0("Backbone Model: Skeleton")}
  if (x$model == "lspar") {model <- paste0("Backbone Model: Local Sparsification")}
  if (x$model == "gspar") {model <- paste0("Backbone Model: Global Sparsification")}
  if (x$model == "simmelian") {model <- paste0("Backbone Model: Simmelian Sparsification")}
  if (x$model == "jaccard") {model <- paste0("Backbone Model: Jaccard Sparsification")}
  if (x$model == "meetmin") {model <- paste0("Backbone Model: Meetmin Sparsification")}
  if (x$model == "geometric") {model <- paste0("Backbone Model: Geometric Sparsification")}
  if (x$model == "hyper") {model <- paste0("Backbone Model: Hypergeometric Sparsification")}
  if (x$model == "degree") {model <- paste0("Backbone Model: Local Degree")}
  if (x$model == "quadrilateral") {model <- paste0("Backbone Model: Quadrilateral Simmelian Sparsification")}
  if (x$model == "custom") {model <- paste0("Backbone Model: Custom unweighted sparsification model")}
  cat(model,"\n\n")
  
  invisible(x)
}
