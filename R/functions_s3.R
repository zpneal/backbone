#' @export
print.backbone <- function(x, ...) {
  #Header
  cat("\n=== SUMMARY OF BACKBONE OBJECT ===\n")
  
  #Source network
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "matrix") | methods::is(x$bipartite, "Matrix"))) {source <- paste0("ORIGINAL:   Bipartite with ", nrow(x$bipartite), " agents, ", ncol(x$bipartite), " artifacts, and ", sum(x$bipartite), " edges")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {source <- paste0("ORIGINAL:   Bipartite with ", as.numeric(table(igraph::V(x$bipartite)$type)[1]), " agents, ", as.numeric(table(igraph::V(x$bipartite)$type)[2]), " artifacts, and ", igraph::gsize(x$bipartite), " edges")}
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "matrix") | methods::is(x$weighted, "Matrix"))) {if (isSymmetric(x$weighted)) {source <- paste0("ORIGINAL:   Weighted with ", nrow(x$weighted), " nodes and ", sum(x$weighted > 0)/2, " edges (weights: ", round(min(x$weighted[which(x$weighted!=0)]),2), " - ", round(max(x$weighted[which(x$weighted!=0)]),2), ")")}}
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "matrix") | methods::is(x$weighted, "Matrix"))) {if (!isSymmetric(x$weighted)) {source <- paste0("ORIGINAL:   Weighted with ", nrow(x$weighted), " nodes and ", sum(x$weighted > 0), " edges (weights: ", round(min(x[which(x$weighted!=0)]),2), " - ", round(max(x$weighted[which(x$weighted!=0)]),2), ")")}}
  if ("weighted" %in% names(x) & (methods::is(x$weighted, "igraph"))) {source <- paste0("ORIGINAL:   Weighted with ", igraph::gorder(x$weighted), " nodes and ", igraph::gsize(x$weighted), " edges (weights: ", round(min(igraph::E(x$weighted)$weight),2), " - ", round(max(igraph::E(x$weighted)$weight),2), ")")}
  if ("unweighted" %in% names(x) & (methods::is(x$unweighted, "matrix") | methods::is(x$unweighted, "Matrix"))) {source <- paste0("ORIGINAL:   Unweighted with ", nrow(x$unweighted), " nodes and ", sum(x$unweighted > 0), " edges")}
  if ("unweighted" %in% names(x) & (methods::is(x$unweighted, "igraph"))) {source <- paste0("ORIGINAL:   Unweighted with ", igraph::gorder(x$unweighted), " nodes and ", igraph::gsize(x$unweighted), " edges")}
  cat(source,"\n\n")
  
  #Projection (if present)
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "matrix") | methods::is(x$bipartite, "Matrix"))) {projection <- paste0("PROJECTION: Weighted with ", nrow(x$projection), " nodes and ", sum(x$projection > 0)/2, " edges (weights: ", round(min(x$projection[which(x$projection!=0)]),2), " - ", round(max(x$projection[which(x$projection!=0)]),2), ")")}
  if ("bipartite" %in% names(x) & (methods::is(x$bipartite, "igraph"))) {projection <- paste0("PROJECTION: Weighted with ", igraph::gorder(x$projection), " nodes and ", igraph::gsize(x$projection), " edges (weights: ", round(min(igraph::E(x$projection)$weight),2), " - ", round(max(igraph::E(x$projection)$weight),2), ")")}
  if ("bipartite" %in% names(x)) {cat(projection,"\n\n")}
  
  #Backbone
  if ("lower" %in% names(x$pvalues)) {type <- "Signed with "} else {type <- "Unweighted with "}
  if (methods::is(x$backbone, "matrix") | methods::is(x$backbone, "Matrix")) {if (isSymmetric(x$backbone)) {backbone <- paste0("BACKBONE:   ", type, nrow(x$backbone), " nodes and ", sum(x$backbone!=0)/2, " edges")}}
  if (methods::is(x$backbone, "matrix") | methods::is(x$backbone, "Matrix")) {if (!isSymmetric(x$backbone)) {backbone <- paste0("BACKBONE:   ", type, nrow(x$backbone), " nodes and ", sum(x$backbone!=0), " edges")}}
  if (methods::is(x$backbone, "igraph")) {backbone <- paste0("BACKBONE:   ", type, igraph::gorder(x$backbone), " nodes and ", igraph::gsize(x$backbone), " edges")}
  cat(backbone,"\n")
  
  #Model
  if (x$model == "sdsm") {model <- paste0("   Model:   Stochastic Degree Sequence with alpha = ", x$alpha)}
  if (x$model == "fdsm") {model <- paste0("   Model:   Fixed Degree Sequence with alpha = ", x$alpha)}
  if (x$model == "fixedrow") {model <- paste0("   Model:   Fixed Row with alpha = ", x$alpha)}
  if (x$model == "fixedcol") {model <- paste0("   Model:   Fixed Column with alpha = ", x$alpha)}
  if (x$model == "fixedfill") {model <- paste0("   Model:   Fixed Fill with alpha = ", x$alpha)}
  if (x$model == "disparity") {model <- paste0("   Model:   Disparity Filter with alpha = ", x$alpha)}
  if (x$model == "lans") {model <- paste0("   Model:   Locally Adaptive Network Sparsification with alpha = ", x$alpha)}
  if (x$model == "mlf") {model <- paste0("   Model:   Marginal Likelihood Filter with alpha = ", x$alpha)}
  if (x$model == "global") {model <- paste0("   Model:   Global Threshold with threshold = ", x$parameter)}
  if (x$model == "skeleton") {model <- paste0("   Model:   Skeleton with parameter = ", x$parameter)}
  if (x$model == "lspar") {model <- paste0("   Model:   Local Sparsification with parameter = ", x$parameter)}
  if (x$model == "gspar") {model <- paste0("   Model:   Global Sparsification with parameter = ", x$parameter)}
  if (x$model == "simmelian") {model <- paste0("   Model:   Simmelian with parameter = ", x$parameter)}
  if (x$model == "jaccard") {model <- paste0("   Model:   Jaccard with parameter = ", x$parameter)}
  if (x$model == "meetmin") {model <- paste0("   Model:   Meetmin with parameter = ", x$parameter)}
  if (x$model == "geometric") {model <- paste0("   Model:   Geometric with parameter = ", x$parameter)}
  if (x$model == "hyper") {model <- paste0("   Model:   Hypergeometric with parameter = ", x$parameter)}
  if (x$model == "degree") {model <- paste0("   Model:   Local Degree with parameter = ", x$parameter)}
  if (x$model == "quadrilateral") {model <- paste0("   Model:   Quadrilateral Simmelian with parameter = ", x$parameter)}
  if (x$model == "custom") {model <- paste0("   Model:   Custom unweighted backbone with parameter = ", x$parameter)}
  cat(model,"\n")
  
  #Call
  call <- paste0("    Call:   ", deparse(x$call, width.cutoff = 500L))
  cat(call,"\n\n")
  
  #Narrative
  cat("DESCRIPTION -\n")
  cat(x$narrative,"\n")
  
  #Footer
  cat("==================================")
  invisible(x)
}
