#' backbone: Extracts the Backbone from Graphs
#'
#' @description Provides methods for extracting an unweighted and sparse subgraph (i.e., a backbone)
#'    that contains only the most "important" edges from:
#'    * a weighted network using [backbone_from_weighted()]
#'    * a weighted projection of a bipartite network or hypergraph using [backbone_from_projection()]
#'    * an unweighted network using [backbone_from_unweighted()]
#'
#'    For a detailed illustration of these methods, please see \href{../doc/backbone.html}{the vignette} using `vignette("backbone")`.
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#'
#' @name backbone
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
