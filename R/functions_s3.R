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
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "matrix"))) {source <- paste0("Source:   matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "Matrix"))) {source <- paste0("Source:   matrix-class weighted network")}
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "igraph"))) {source <- paste0("Source:   igraph-class weighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$unweighted, "matrix"))) {source <- paste0("Source:   matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$unweighted, "Matrix"))) {source <- paste0("Source:   matrix-class unweighted network")}
  if ("unweighted" %in% names(x) & (methods::is(x$unweighted, "igraph"))) {source <- paste0("Source:   igraph-class unweighted network")}
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

#' @export
summary.backbone <- function(object, ...) {

  #Column headers
  cat("\n")
  cat("   NETWORK       TYPE       NODES      EDGES     MIN     MAX\n")
  cat(paste(strrep("\u2500",60),"\n",collapse=""))

  #Source
  if ("bipartite" %in% names(object)) {
    if (methods::is(object$bipartite,"igraph")) {source <- igraph::as_biadjacency_matrix(object$bipartite, sparse = FALSE)} else {source <- as.matrix(object$bipartite)}
    type <- "Bipartite"
    nodes <- paste0(dim(source)[1], "/", dim(source)[2])
    edges1 <- sum(source)
    min <- "---"
    max <- "---"
  }
  if ("weighted" %in% names(object)) {
    if (methods::is(object$weighted,"igraph")) {source <- igraph::as_adjacency_matrix(object$weighted, sparse = FALSE)} else {source <- as.matrix(object$weighted)}
    type <- "Weighted"
    nodes <- as.character(dim(source)[1])
    if (isSymmetric(source)) {edges1 <- sum(source!=0)/2} else {edges1 <- sum(source!=0)}
    min <- as.character(round(min(source[source!=0]),2))
    max <- as.character(round(max(source[source!=0]),2))
  }
  if ("unweighted" %in% names(object)) {
    if (methods::is(object$unweighted,"igraph")) {source <- igraph::as_adjacency_matrix(object$unweighted, sparse = FALSE)} else {source <- as.matrix(object$unweighted)}
    type <- "Unweighted"
    nodes <- as.character(dim(source)[1])
    edges1 <- sum(source)/2
    min <- "---"
    max <- "---"
  }
  cat(sprintf("%10s", "Original"), sprintf("%10s", type), sprintf("%11s", nodes), sprintf("%10i", edges1), sprintf("%7s", min), sprintf("%7s", max), "\n")

  #Projection (if present)
  if ("projection" %in% names(object)) {
    if (methods::is(object$projection,"igraph")) {projection <- igraph::as_adjacency_matrix(object$projection, sparse = FALSE)} else {projection <- as.matrix(object$projection)}
    diag(projection) <- 0
    type <- "Weighted"
    nodes <- dim(projection)[1]
    edges1 <- sum(projection!=0)/2
    min <- round(min(projection[projection!=0]),2)
    max <- round(max(projection[projection!=0]),2)
    cat(sprintf("%10s", "Projection"), sprintf("%10s", type), sprintf("%11s", nodes), sprintf("%10i", edges1), sprintf("%7s", min), sprintf("%7s", max), "\n")
  }

  #Backbone
  if (methods::is(object$backbone,"igraph")) {backbone <- igraph::as_adjacency_matrix(object$backbone, sparse = FALSE)} else {backbone <- as.matrix(object$backbone)}
  type <- "Unweighted"
  nodes <- dim(backbone)[1]
  if (isSymmetric(backbone)) {edges2 <- sum(backbone!=0)/2} else {edges2 <- sum(backbone!=0)}
  min <- "---"
  max <- "---"
  cat(sprintf("%10s", "Backbone"), sprintf("%10s", type), sprintf("%11s", nodes), sprintf("%10i", edges2), sprintf("%7s", min), sprintf("%7s", max), "\n")

  #Footer
  cat(paste(strrep("\u2500",60),"\n",collapse=""))

  #Model
  if (object$model == "sdsm") {model <- paste0("Stochastic Degree Sequence Model (\u03B1 = ", object$alpha, ")")}
  if (object$model == "fdsm") {model <- paste0("Fixed Degree Sequence Model (\u03B1 = ", object$alpha, ")")}
  if (object$model == "fixedrow") {model <- paste0("Fixed Row Model (\u03B1 = ", object$alpha, ")")}
  if (object$model == "fixedcol") {model <- paste0("Fixed Column Model (\u03B1 = ", object$alpha, ")")}
  if (object$model == "fixedfill") {model <- paste0("Fixed Fill Model (\u03B1 = ", object$alpha, ")")}
  if (object$model == "disparity") {model <- paste0("Disparity Filter (\u03B1 = ", object$alpha, ")")}
  if (object$model == "lans") {model <- paste0("Locally Adaptive Network Sparsification (\u03B1 = ", object$alpha, ")")}
  if (object$model == "mlf") {model <- paste0("Marginal Likelihood Filter (\u03B1 = ", object$alpha, ")")}
  if (object$model == "global") {model <- paste0("Global Threshold (parameter = ", object$parameter, ")")}
  if (object$model == "skeleton") {model <- paste0("Skeleton (parameter = ", object$parameter, ")")}
  if (object$model == "lspar") {model <- paste0("Local Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "gspar") {model <- paste0("Global Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "simmelian") {model <- paste0("Simmelian Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "jaccard") {model <- paste0("Jaccard Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "meetmin") {model <- paste0("Meetmin Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "geometric") {model <- paste0("Geometric Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "hyper") {model <- paste0("Hypergeometric Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "degree") {model <- paste0("Local Degree (parameter = ", object$parameter, ")")}
  if (object$model == "quadrilateral") {model <- paste0("Quadrilateral Simmelian Sparsification (parameter = ", object$parameter, ")")}
  if (object$model == "custom") {model <- paste0("Custom unweighted sparsification model (parameter = ", object$parameter, ")")}
  cat(model,"\n")

  #Edge reduction
  change <- round(((edges1 - edges2) / edges1) * 100,1)
  cat(paste0(edges1 - edges2, " (", change, "%) edges removed\n\n"))

  #Narrative
  cat("NARRATIVE SUMMARY -\n")
  cat(object$narrative)
  cat("\n")

  invisible(object)
}

#' @export
#' @importFrom graphics plot
plot.backbone <- function(x, ...) {
  mfrow <- graphics::par("mfrow")  #Save existing plot settings
  graphics::par(mfrow = c(1, 2))  #Plot side-by-side

  #Get original network
  if (!is.null(x$projection)) {original <- x$projection}
  if (!is.null(x$weighted)) {original <- x$weighted}
  if (!is.null(x$unweighted)) {original <- x$unweighted}
  if (methods::is(original, "matrix") | methods::is(original, "Matrix")) {original <- igraph::graph_from_adjacency_matrix(original, weighted = TRUE, diag = FALSE, mode = "undirected")}

  #Get backbone network
  backbone <- x$backbone
  if (methods::is(backbone, "matrix") | methods::is(backbone, "Matrix")) {backbone <- igraph::graph_from_adjacency_matrix(backbone, weighted = FALSE, diag = FALSE, mode = "undirected")}

  #Plot original network
  plot(original, main = "Original", ...)

  #Plot backbone network
  plot(backbone, main = "Backbone", ...)

  graphics::par(mfrow = mfrow)  #Restore plot settings
}
