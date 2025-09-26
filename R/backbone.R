#' Extract the backbone from a network
#'
#' \code{backbone()} extracts the backbone from a weighted or unweighted network
#'
#' @param N A network represented as a matrix, \link[Matrix]{Matrix}, or \link[igraph]{igraph} object
#' @param ... Optional arguments
#'
#' @details Given a weighted and/or dense network, the backbone is an sparse and unweighted subgraph
#'    that contains only the most "important" edges.
#'
#'    \code{backbone()} is a wrapper that detects the type of network in `N`, then extracts the backbone
#'    using the appropriate `backbone_from_*()` function:
#'
#'    * If `N` is a weighted network, [backbone_from_weighted()]
#'    * If `N` is a bipartite network or hypergraph, [backbone_from_projection()]
#'    * If `N` is an unweighted network, [backbone_from_unweighted()]
#'
#'    Designed as a user-friendly wrapper, `backbone()` uses the defaults for the underlying functions,
#'    which should work in many cases, and which can be modified by passing optional arguments. Unlike
#'    the underlying functions, by default `backbone()` displays narrative text describing what it did
#'    so that you know.
#'
#'    For details about the backbone models, see the documentation for the underlying functions above. For
#'    an overview of the package with examples, please see the \href{../doc/backbone.html}{Introduction to
#'    Backbone} using `vignette("backbone")`. For a detailed empirical example, please see the
#'    \href{../doc/senate.html}{U.S. Senate Example} using `vignette("senate108")`.
#'
#' @return A backbone in the same class as `N`
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#'
#' @export
#'
#' @examples
#' N <- igraph::sample_gnp(100, .3)  #A random unweighted network
#' backbone(N)
#'
backbone <- function(N, ...) {

  #### Coerce input to matrix ####
  if (methods::is(N,"matrix")) {mat <- N}  #matrix --> matrix
  if (methods::is(N,"Matrix")) {mat <- as.matrix(N)}  #Matrix --> matrix
  if (methods::is(N,"igraph")) {
    if (igraph::is_bipartite(N)) {
      if ("weight" %in% igraph::edge_attr_names(N)) {mat <- igraph::as_biadjacency_matrix(N, names = FALSE, sparse = FALSE, attr = "weight")}  #weighted bipartite igraph --> weighted incidence matrix
      if (!("weight" %in% igraph::edge_attr_names(N))) {mat <- igraph::as_biadjacency_matrix(N, names = FALSE, sparse = FALSE)}  #unweighted bipartite igraph --> binary incidence
    }
    if (!igraph::is_bipartite(N)) {
      if ("weight" %in% igraph::edge_attr_names(N)) {mat <- igraph::as_adjacency_matrix(N, names = FALSE, sparse = FALSE, attr = "weight")}  #weighted unipartite igraph --> weighted adjacency
      if (!("weight" %in% igraph::edge_attr_names(N))) {mat <- igraph::as_adjacency_matrix(N, names = FALSE, sparse = FALSE)}  #unweighted unipartite igraph --> binary adjacency
    }
  }

  #Detect and extract backbone
  if (is.numeric(mat) & dim(mat)[1]==dim(mat)[2] & all(mat %in% c(0,1))) {  #Numeric, square, binary
    return(backbone_from_unweighted(N, narrative = TRUE, ...))
  } else if (is.numeric(mat) & dim(mat)[1]==dim(mat)[2] & !all(mat %in% c(0,1))) {  #Numeric, square, valued
      return(backbone_from_weighted(N, narrative = TRUE, ...))
  } else if (is.numeric(mat) & dim(mat)[1]!=dim(mat)[2] & all(mat %in% c(0,1))) {  #Numeric, non-square, binary
      return(backbone_from_projection(N, narrative = TRUE, ...))
  } else {stop("`N` does not seem to contain a supported type of network.")}

}

## usethis namespace: start
#' @useDynLib backbone, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
