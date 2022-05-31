#' backbone: Extracts the Backbone from Graphs
#'
#' @description Provides methods for extracting from an unweighted and sparse subgraph (i.e., a backbone)
#'    that contains only the most "important" edges in a weighted bipartite projection, a non-projection
#'    weighted network, or an unweighted network.
#'
#'    Available backbone extraction functions include:
#'    \itemize{
#'    \item For weighted bipartite projections of weighted bipartite networks: [osdsm()].
#'    \item For weighted bipartite projections of binary bipartite networks: [fixedfill()], [fixedrow()], [fixedcol()], [sdsm()], and [fdsm()].
#'    \item For non-projection weighted networks: [global()], [disparity()].
#'    \item For unweighted networks: [sparsify()], [sparsify.with.skeleton()], [sparsify.with.gspar()], [sparsify.with.lspar()], [sparsify.with.simmelian()], [sparsify.with.jaccard()], [sparsify.with.meetmin()], [sparsify.with.geometric()], [sparsify.with.hypergeometric()], [sparsify.with.localdegree()], [sparsify.with.quadrilateral()].
#'    \item For all networks: [backbone.suggest()] will examine the data and suggest an appropriate backbone function
#'    }
#'
#'    The package also includes some utility functions:
#'    \itemize{
#'    \item [fastball()] - Fast marginal-preserving randomization of binary matrices
#'    \item [bicm()] - Compute probabilities under the bipartite configuration model
#'    }
#'
#' For additional documentation and background on the package functions, see \href{../doc/backbone.html}{\code{vignette("backbone")}}.
#'     For updates, papers, presentations, and other backbone news, please see \href{https://www.zacharyneal.com/backbone}{www.rbackbone.net}
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#'
#' @docType package
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
