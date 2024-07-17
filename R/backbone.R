#' backbone: Extracts the Backbone from Graphs
#'
#' @description Provides methods for extracting from an unweighted and sparse subgraph (i.e., a backbone)
#'    that contains only the most "important" edges in a weighted bipartite projection, a non-projection
#'    weighted network, or an unweighted network.
#'
#' @references {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#'
#' @docType package
#' @aliases backbone-package
#' @name backbone
NULL

## usethis namespace: start
#' @useDynLib backbone, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
