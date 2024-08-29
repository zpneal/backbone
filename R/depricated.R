#' Depricated function to extract Stochastic Degree Sequence Model (SDSM) backbone
#'
#' \code{sdsm()} was replaced by \code{backbone_from_bipartite(model="sdsm")}.
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

#' Depricated function to extract Fixed Row backbone
#'
#' \code{fixedrow()} was replaced by \code{backbone_from_bipartite(model = "fixedrow")}.
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
fixedrow <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_bipartite(model = \"fixedrow\")")

  return(
    backbone_from_bipartite(B,
                            model = "fixedrow",
                            alpha = alpha,
                            signed = signed,
                            mtc = mtc,
                            missing_as_zero = missing.as.zero,
                            narrative = TRUE)
  )
}

#' Depricated function to extract Fixed Column backbone
#'
#' \code{fixedcol()} was replaced by \code{backbone_from_bipartite(model = "fixedcol")}.
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
fixedcol <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_bipartite(model = \"fixedcol\")")

  return(
    backbone_from_bipartite(B,
                            model = "fixedcol",
                            alpha = alpha,
                            signed = signed,
                            mtc = mtc,
                            missing_as_zero = missing.as.zero,
                            narrative = TRUE)
  )
}

#' Depricated function to extract Fixed Fill backbone
#'
#' \code{fixedfill()} was replaced by \code{backbone_from_bipartite(model = "fixedfill")}.
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
fixedfill <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_bipartite(model = \"fixedfill\")")

  return(
    backbone_from_bipartite(B,
                            model = "fixedfill",
                            alpha = alpha,
                            signed = signed,
                            mtc = mtc,
                            missing_as_zero = missing.as.zero,
                            narrative = TRUE)
  )
}

#' Depricated function to extract Fixed Degree Sequence Model (FDSM) backbone
#'
#' \code{fdsm()} was replaced by \code{backbone_from_bipartite(model = "fdsm")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \code{\link{igraph}} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#' @param trials numeric: the number of bipartite graphs generated using fastball to approximate the edge weight distribution
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
fdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE, trials = NULL){
  .Deprecated("backbone_from_bipartite(model = \"fdsm\")")

  return(
    backbone_from_bipartite(B,
                            model = "fdsm",
                            alpha = alpha,
                            signed = signed,
                            mtc = mtc,
                            missing_as_zero = missing.as.zero,
                            narrative = TRUE)
  )
}
