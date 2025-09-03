#' backbone: Extracts the Backbone from Graphs
#'
#' @description Provides methods for extracting an unweighted and sparse subgraph (i.e., a backbone)
#'    that contains only the most "important" edges from:
#'    * a weighted network using [backbone_from_weighted()]
#'    * a weighted bipartite projection using [backbone_from_bipartite()]
#'    * an unweighted network using [backbone_from_unweighted()]
#'
#'    For a detailed illustration of these methods, please see the vignette using [vignette("backbone")].
#'
#' @references {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
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
