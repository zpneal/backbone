#' @export
print.backbone <- function(x, ...) {

  #Call
  cat("Call -\n")
  cat(deparse(x$call, width.cutoff = 500L))
  cat("\n\n")

  #Source
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source:   matrix-class bipartite network")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source:   Matrix-class bipartite network")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source:   igraph-class bipartite network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source:   matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source:   matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source:   igraph-class weighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "matrix"))) {source <- paste0("Source:   matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "Matrix"))) {source <- paste0("Source:   matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("Source:   igraph-class unweighted network")}
  cat(source,"\n")

  #Result
  if (methods::is(x$backbone, "matrix")) {result <- paste0("Backbone: matrix-class unweighted network")}
  if (methods::is(x$backbone, "Matrix")) {result <- paste0("Backbone: Matrix-class unweighted network")}
  if (methods::is(x$backbone, "igraph")) {result <- paste0("Backbone: igraph-class unweighted network")}
  cat(result,"\n")

  #Model
  if (x$model == "sdsm") {model <- paste0("Model:    Stochastic Degree Sequence Model")}
  if (x$model == "fdsm") {model <- paste0("Model:    Fixed Degree Sequence Model")}
  if (x$model == "fixedrow") {model <- paste0("Model:    Fixed Row Model")}
  if (x$model == "fixedcol") {model <- paste0("Model:    Fixed Column Model")}
  if (x$model == "fixedfill") {model <- paste0("Model:    Fixed Fill Model")}
  if (x$model == "disparity") {model <- paste0("Model:    Disparity Filter")}
  if (x$model == "lans") {model <- paste0("Model:    Locally Adaptive Network Sparsification")}
  if (x$model == "mlf") {model <- paste0("Model:    Marginal Likelihood Filter")}
  if (x$model == "global") {model <- paste0("Model:    Global Threshold")}
  if (x$model == "skeleton") {model <- paste0("Model:    Skeleton")}
  if (x$model == "lspar") {model <- paste0("Model:    Local Sparsification")}
  if (x$model == "gspar") {model <- paste0("Model:    Global Sparsification")}
  if (x$model == "simmelian") {model <- paste0("Model:    Simmelian Sparsification")}
  if (x$model == "jaccard") {model <- paste0("Model:    Jaccard Sparsification")}
  if (x$model == "meetmin") {model <- paste0("Model:    Meetmin Sparsification")}
  if (x$model == "geometric") {model <- paste0("Model:    Geometric Sparsification")}
  if (x$model == "hyper") {model <- paste0("Model:    Hypergeometric Sparsification")}
  if (x$model == "degree") {model <- paste0("Model:    Local Degree")}
  if (x$model == "quadrilateral") {model <- paste0("Model:    Quadrilateral Simmelian Sparsification")}
  if (x$model == "custom") {model <- paste0("Model:    Custom unweighted sparsification model")}
  cat(model,"\n\n")

  invisible(x)
}
