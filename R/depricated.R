#' Depricated function to extract SDSM backbone
#'
#' \code{sdsm()} was replaced by \code{backbone_from_bipartite()}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \code{\link{igraph}} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
sdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_bipartite(model = \"sdsm\")")

  return(
  backbone_from_bipartite(B,
                          model = "sdsm",
                          alpha = alpha,
                          signed = signed,
                          mtc = mtc,
                          missing_as_zero = missing.as.zero,
                          narrative = TRUE)
  )
}
