#' Extract the backbone from a bipartite projection
#'
#' `backbone_from_bipartite` extracts the unweighted backbone from the weighted projection of a bipartite network.
#'
#' @param B An unweighted bipartite network as a matrix, \code{\link{Matrix}}, or \code{\link{igraph}} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param model string: backbone model. This must be onf of: \code{"sdsm"}, \code{"fdsm"}, \code{"fixedrow"}, \code{"fixedcol"}, or \code{"fixedfill"}
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction; can be either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param missing_as_zero boolean: treat missing edges be treated as edges with zero weight and test them for significance
#' @param only_pvalues boolean: return only a matrix of edgewise p-values; all parameters except `model` are ignored
#' @param narrative boolean: display suggested text & citations
#'
#' @return TBD
#' @export
#'
#' @examples
#' #TBD
backbone_from_bipartite <- function(B,
                                    alpha = 0.05,
                                    model = "sdsm",
                                    signed = FALSE,
                                    mtc = "none",
                                    missing_as_zero = FALSE,
                                    only_pvalues = FALSE,
                                    narrative = FALSE) {

}
