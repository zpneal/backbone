#' Extract the backbone from an unweighted, undirected network
#'
#' \code{backbone_from_unweighted()} extracts the unweighted backbone from an unweighted, undirected network
#'
#' @param U An unweighted, undirected network as an adjacency matrix or an unweighted unipartite \code{igraph} object
#' @param model string: backbone model
#' @param parameter real: filtering parameter
#' @param escore string: Method for scoring edges' importance
#' @param normalize string: Method for normalizing edge scores
#' @param filter string: Type of filter to apply
#' @param umst boolean: TRUE if the backbone should include the union of minimum spanning trees, to ensure connectivity
#' @param narrative boolean: display suggested text & citations
#' @param return string: return either only the \code{"backbone"} or \code{"everything"}
#'
#' @details
#' The \code{backbone_from_unweighted} function extracts the backbone from an unweighted unipartite network. The backbone is an
#' unweighted unipartite network that contains only edges preserved by a backbone model.
#'
#' The following backbone models are available using the \code{model} parameter:
#' * \code{skeleton} - Karger's (1999) Skeleton backbone
#' * \code{gspar} - Satuluri et al's (2011) Global Sparsification backbone
#' * \code{lspar} - Satuluri et al's (2011) Local Sparsification backbone
#' * \code{simmelian} - Nick et al's (2013) Simmelian backbone
#' * \code{jaccard} - Goldberg and Roth's (2003) Jaccard backbone
#' * \code{meetmin} - Goldberg and Roth's (2003) MeetMin backbone
#' * \code{geometric} - Goldberg and Roth's (2003) Geometric backbone
#' * \code{hyper} - Goldberg and Roth's (2003) Hypergeometric backbone
#' * \code{degree} - Hamann et al.'s (2016) Local Degree backbone
#' * \code{quadrilateral} - Nocaj et al.'s (2015) Quadrilateral Simmelian backbone
#' * \code{custom} - A custom backbone model specified by \code{escore}, \code{normalize},\code{filter}, and \code{umst}
#'
#' The \code{escore} parameter determines how an unweighted edge's importance is calculated.
#' Unless noted below, scores are symmetric and larger values represent more important edges.
#' * \code{random}: a random number drawn from a uniform distribution
#' * \code{betweenness}: edge betweenness
#' * \code{triangles}: number of triangles that include the edge
#' * \code{jaccard}: jaccard similarity coefficient of the neighborhoods of an edge's endpoints, or alternatively, triangles normalized by the size of the union of the endpoints neighborhoods
#' * \code{dice}: dice similarity coefficient of the neighborhoods of an edge's endpoints
#' * \code{quadrangles}: number of quadrangles that include the edge
#' * \code{quadrilateral}: geometric mean normalization of quadrangles
#' * \code{degree}: degree of neighbor to which an edge is adjacent (asymmetric)
#' * \code{meetmin}: triangles normalized by the smaller of the endpoints' neighborhoods' sizes
#' * \code{geometric}: triangles normalized by the product of the endpoints' neighborhoods' sizes
#' * \code{hypergeometric}: probability of the edge being included at least as many triangles if edges were random, given the size of the endpoints' neighborhoods (smaller is more important)
#'
#' The \code{normalize} parameter determines whether edge scores are normalized.
#' * \code{none}: no normalization is performed
#' * \code{rank}: scores are normalized by neighborhood rank, such that the strongest edge in a node's neighborhood is ranked 1 (asymmetric)
#' * \code{embeddedness}: scores are normalized using the maximum Jaccard coefficient of the top k-ranked neighbors of each endpoint, for all k
#'
#' The \code{filter} parameter determines how edges are filtered based on their (normalized) edge scores.
#' * \code{threshold}: Edges with scores >= `s` are retained in the backbone
#' * \code{proportion}: Specifies the approximate proportion of edges to retain in the backbone
#' * \code{degree}: Retains each node's d^`s` most important edges, where d is the node's degree (requires that \code{normalize = "rank"})
#' * \code{disparity}: Applies the disparity filter using [backbone_from_weighted()]
#' * \code{lans}: Applies locally adaptive network sparsification using [backbone_from_weighted()]
#' * \code{mlf}: Applies the marginal likelihood filter using [backbone_from_weighted()]
#'
#' @return If \code{return = "backbone"}, a backbone in the same class as \code{B}. If \code{return = "everything"}, then the backbone
#' is returned as an element in a list that also includes the original weighted network, a narrative description, and (for statistical
#' backbone models) the edgewise p-values.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references skeleton: {Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. *Mathematics of Operations Research, 24*, 383-413. \doi{10.1287/moor.24.2.383}}
#' @references gspar and lspar: {Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). \doi{10.1145/1989323.1989399}}
#' @references simmelian: {Nick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). \doi{10.1145/2492517.2492569}}
#' @references jaccard, meetmin, geometric, hyper: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @references degree: {Hamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. *Social Network Analysis and Mining, 6*, 22. \doi{10.1007/s13278-016-0332-2}}
#' @references quadrilateral: {Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. *Journal of Graph Algorithms and Applications, 19*, 595-618. \doi{10.7155/jgaa.00370}}
#' @export
#'
#' @examples
backbone_from_unweighted <- function(U,
                                     model = "degree",
                                     parameter = 0.6,
                                     escore,
                                     normalize,
                                     filter,
                                     umst,
                                     narrative = TRUE,
                                     return = "backbone") {

  #### Check parameters ####
  #All models
  if (!(model %in% c("custom", "skeleton", "gspar", "lspar", "simmelian", "jaccard", "meetmin", "geometric", "hyper", "degree", "quadrilateral"))) {stop("`model` must be one of: \"custom\", \"skeleton\", \"gspar\", \"lspar\", \"simmelian\", \"jaccard\", \"meetmin\", \"geometric\", \"hyper\", \"degree\", \"quadrilateral\"")}
  if (!is.numeric(parameter)) {stop("`parameter` must be a numeric value between 0 and 1")}
  if (parameter < 0 | parameter > 1) {stop("`parameter` must be a numeric value between 0 and 1")}
  if (!is.logical(narrative)) {stop("`narrative` must be either TRUE or FALSE")}
  if (!(return %in% c("backbone", "everything"))) {stop("`return` must be one of: \"backbone\", \"everything\"")}

  #If existing model specification, set model parameters
  if (model == "skeleton") {escore <- "random"; normalize <- "none"; filter <- "proportion"; umst <- FALSE}
  if (model == "gspar") {escore <- "jaccard"; normalize <- "none"; filter <- "proportion"; umst <- FALSE}
  if (model == "lspar") {escore <- "jaccard"; normalize <- "rank"; filter <- "degree"; umst <- FALSE}
  if (model == "simmelian") {escore <- "triangles"; normalize <- "embeddedness"; filter <- "threshold"; umst <- FALSE}
  if (model == "jaccard") {escore <- "jaccard"; normalize <- "none"; filter <- "threshold"; umst <- FALSE}
  if (model == "meetmin") {escore <- "meetmin"; normalize <- "none"; filter <- "threshold"; umst <- FALSE}
  if (model == "geometric") {escore <- "geometric"; normalize <- "none"; filter <- "threshold"; umst <- FALSE}
  if (model == "hyper") {escore <- "hypergeometric"; normalize <- "none"; filter <- "threshold"; umst <- FALSE}
  if (model == "degree") {escore <- "degree"; normalize <- "rank"; filter <- "degree"; umst <- FALSE}
  if (model == "quadrilateral") {escore <- "quadrilateral"; normalize <- "embeddedness"; filter <- "threshold"; umst <- TRUE}

  #If custom model specification, check model parameters
  if (model == "custom") {
    if (!(escore %in% c("random", "betweenness", "triangles", "jaccard", "dice", "quadrangles", "quadrilateral", "degree", "meetmin", "geometric" , "hypergeometric"))) {stop("`escore` must be one of: \"random\", \"betweenness\", \"triangles\", \"jaccard\", \"dice\", \"quadrangles\", \"quadrilateral\", \"degree\", \"meetmin\", \"geometric\" , \"hypergeometric\"")}
    if (!(normalize %in% c("none", "rank", "embeddedness"))) {stop("`normalize` must be one of: \"none\", \"rank\", \"embeddedness\"")}
    if (!(filter %in% c("threshold", "proportion", "degree", "disparity", "lans", "mlf"))) {stop("`filter` must be one of: \"threshold\", \"proportion\", \"degree\", \"lans\", \"mlf\"")}
    if (!is.logical(umst)) {stop("`umst` must be either TRUE or FALSE")}
  }

  #### Check and format input ####
  #Check that input is a weighted adjacency matrix or weighted unipartite igraph
  if (!methods::is(U,"matrix") & !methods::is(U,"igraph")) {stop("`U` must be an adjacency matrix or igraph object")}

  if (methods::is(U,"matrix")) {
    if (dim(U)[1] != dim(U)[2]) {stop("`U` must be a symmetric adjacency matrix")}
    if (!all(U %in% c(0,1))) {stop("The entries of `U` must be either 0 or 1")}
    if (!isSymmetric(U)) {stop("`U` must be a symmetric adjacency matrix")}
  }

  if (methods::is(U,"igraph")) {
    if (igraph::is_bipartite(U)) {stop("`U` must be an undirected unipartite igraph object")}
    if (igraph::is_directed(U)) {stop("`U` must be an undirected unipartite igraph object")}
    if ("weight" %in% igraph::edge_attr_names(U)) {stop("An edge weight attribute is present in `U`, but will be ignored")}
  }

  #Convert input to adjacency matrix
  if (methods::is(U,"matrix")) {A <- U}  #matrix --> matrix
  if (methods::is(U,"igraph")) {A <- igraph::as_adjacency_matrix(U, names = FALSE, sparse = FALSE)}

  #### Compute edge scores ####
  G <- .escore(A, escore = escore)

  #### Apply edge score normalization ####
  G <- .normalize(G, normalize = normalize)

  #### Apply filter ####
  G <- .filter(G, filter = filter, parameter = parameter)

  #### Symmetrize ####  ==> REQUIRES TESTING
  G[lower.tri(G)] <- pmax(G[lower.tri(G)],t(G)[lower.tri(t(G))])
  G[upper.tri(G)] <- t(G)[upper.tri(G)]

  #### Add UMST #### ==> REQUIRES TESTING
  if (umst) {
    tree <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")  #Convert original to igraph
    tree <- igraph::mst(tree)  #Find the UMST
    tree <- igraph::as_adjacency_matrix(tree, sparse = FALSE)  #Convert back to matrix
    G <- (G | tree)*1  #Include an edge if it is in either the sparsified graph or the tree
  }

  #### FOLLOW OTHER CODE AS TEMPLATE

}
