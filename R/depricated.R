#' Depricated function to extract Stochastic Degree Sequence Model (SDSM) backbone
#'
#' \code{sdsm()} was replaced by \code{backbone_from_projection(model="sdsm")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
sdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_projection(model = \"sdsm\")")

  return(
  backbone_from_projection(B,
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
#' \code{fixedrow()} was replaced by \code{backbone_from_projection(model = "fixedrow")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
fixedrow <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_projection(model = \"fixedrow\")")

  return(
    backbone_from_projection(B,
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
#' \code{fixedcol()} was replaced by \code{backbone_from_projection(model = "fixedcol")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
fixedcol <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_projection(model = \"fixedcol\")")

  return(
    backbone_from_projection(B,
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
#' \code{fixedfill()} was replaced by \code{backbone_from_projection(model = "fixedfill")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
fixedfill <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_projection(model = \"fixedfill\")")

  return(
    backbone_from_projection(B,
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
#' \code{fdsm()} was replaced by \code{backbone_from_projection(model = "fdsm")}.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing.as.zero boolean: treat missing edges as edges with zero weight and test them for significance
#' @param narrative boolean: display suggested text & citations
#' @param trials numeric: the number of bipartite graphs generated using fastball to approximate the edge weight distribution
#'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
fdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE, trials = NULL){
  .Deprecated("backbone_from_projection(model = \"fdsm\")")

  return(
    backbone_from_projection(B,
                            model = "fdsm",
                            alpha = alpha,
                            signed = signed,
                            mtc = mtc,
                            missing_as_zero = missing.as.zero,
                            narrative = TRUE,
                            trials = trials)
  )
}

#' Depricated function to extract disparity filter backbone
#'
#' \code{disparity()} was replaced by \code{backbone_from_weighted(model = "disparity")}.
#'
#' @param W A positively-weighted unipartite graph as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by [p.adjust()].
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' #'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
disparity <- function(W, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_weighted(model = \"disparity\")")

  return(
    backbone_from_weighted(W,
                           model = "disparity",
                           alpha = alpha,
                           signed = signed,
                           mtc = mtc,
                           missing_as_zero = missing.as.zero,
                           narrative = TRUE)
  )
}

#' Depricated function to extract LANS backbone
#'
#' \code{lans()} was replaced by \code{backbone_from_weighted(model = "lans")}.
#'
#' @param W A positively-weighted unipartite graph as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by [p.adjust()].
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' #'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
lans <- function(W, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_weighted(model = \"lans\")")

  return(
    backbone_from_weighted(W,
                           model = "lans",
                           alpha = alpha,
                           signed = signed,
                           mtc = mtc,
                           missing_as_zero = missing.as.zero,
                           narrative = TRUE)
  )
}

#' Depricated function to extract MLF backbone
#'
#' \code{mlf()} was replaced by \code{backbone_from_weighted(model = "mlf")}.
#'
#' @param W A positively-weighted unipartite graph as a binary incidence matrix or a binary bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by [p.adjust()].
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' #'
#' @details
#' See backbone v2.1.4 for original documentation
#'
#' @export
mlf <- function(W, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE){
  .Deprecated("backbone_from_weighted(model = \"mlf\")")

  return(
    backbone_from_weighted(W,
                           model = "mlf",
                           alpha = alpha,
                           signed = signed,
                           mtc = mtc,
                           missing_as_zero = missing.as.zero,
                           narrative = TRUE)
  )
}