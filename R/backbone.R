#' The backbone package
#'
#' backbone is an R package for extracting network backbones.
#' 
#' @name backbone-package
#' @keywords internal
#' @aliases backbone-package backbone
#'
#' @description The backbone package implements methods for extracting an unweighted and sparse network
#'    (i.e., a backbone) that contains only the most "important" edges from:
#'    * a weighted network using [backbone_from_weighted()]
#'    * a weighted projection of a bipartite network or hypergraph using [backbone_from_projection()]
#'    * an unweighted network using [backbone_from_unweighted()]
#'
#'    For an overview of the package with examples, please see the \href{../doc/backbone.html}{Introduction to Backbone}
#'    using `vignette("backbone")`. For a detailed empirical example, please see the \href{../doc/senate.html}{U.S. Senate Example}
#'    using `vignette("senate108")`.
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
"_PACKAGE"
NULL

## usethis namespace: start
#' @useDynLib backbone, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
