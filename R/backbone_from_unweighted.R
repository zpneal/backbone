#' Extract the backbone from a weighted network
#'
#' \code{backbone_from_unweighted()} extracts the unweighted backbone from an unweighted network
#'
#' @param U An unweighted network as an adjacency matrix or an unweighted unipartite \code{igraph} object
#' @param model string: backbone model
#' @param escore string: Method for scoring edges' importance
#' @param normalize string: Method for normalizing edge scores
#' @param filter string: Type of filter to apply
#' @param symmetrize boolean: TRUE if the result should be symmetrized
#' @param umst boolean: TRUE if the backbone should include the union of minimum spanning trees, to ensure connectivity
#' @param parameter real: parameter used to control structural backbone models (see details)
#' @param narrative boolean: display suggested text & citations
#' @param return string: return either only the \code{"backbone"} or \code{"everything"}
#'
#' @details
#' The \code{backbone_from_unweighted} function extracts the backbone from an unweighted unipartite network. The backbone is an
#' unweighted unipartite network that contains only edges preserved by a backbone model. A model can be chosen from a list of models
#' described in the literature using the \code{model} parameter, or a custom modelcan be specified using the \code{escore},
#' \code{normalize},\code{filter}, \code{symmetrize}, and \code{umst} parameters.
#'
#' The following backbone models are available using the \code{model} parameter:
#' * skeleton - Karger's (1999) Skeleton backbone
#' * gspar - Satuluri et al's (2011) Global Sparsification backbone
#' * lspar - Satuluri et al's (2011) Local Sparsification backbone
#' * simmelian - Nick et al's (2013) Simmelian backbone
#' * jaccard - Goldberg and Roth's (2003) Jaccard backbone
#' * meetmin - Goldberg and Roth's (2003) MeetMin backbone
#' * geometric - Goldberg and Roth's (2003) Geometric backbone
#' * hypergeometric - Goldberg and Roth's (2003) Hypergeometric backbone
#' * localdegree - Hamann et al.'s (2016) Local Degree backbone
#' * quadrilateral - Nocaj et al.'s (2015) Quadrilateral Simmelian backbone
#'
#' The \code{escore} parameter determines how an unweighted edge's importance is calculated.
#' Unless noted below, scores are symmetric and larger values represent more important edges.
#' * \code{random}: a random number drawn from a uniform distribution
#' * \code{betweenness}: edge betweenness
#' * \code{triangles}: number of triangles that include the edge
#' * \code{jaccard}: jaccard similarity coefficient of the neighborhoods of an edge's endpoints, or alternatively, triangles normalized by the size of the union of the endpoints neighborhoods
#' * \code{dice}: dice similarity coefficient of the neighborhoods of an edge's endpoints
#' * \code{quadrangles}: number of quadrangles that include the edge
#' * \code{quadrilateral embeddedness}: geometric mean normalization of quadrangles
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
#' * \code{disparity}: Applies the disparity filter using [disparity()]
#'
#' Using \code{escore == "degree"} or \code{normalize == "rank"} can yield an assymmetric network. When \code{symmetrize == TRUE} (default),
#'   after applying a filter, the network is symmetrized by such that i-j if i->j or i<-j.
#'
#' @return If \code{return = "backbone"}, a backbone in the same class as \code{B}. If \code{return = "everything"}, then the backbone
#' is returned as an element in a list that also includes the original weighted network, a narrative description, and (for statistical
#' backbone models) the edgewise p-values.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references skeleton: {Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. *Mathematics of Operations Research, 24*, 383-413. \doi{10.1287/moor.24.2.383}}
#' @references gspar and lspat: {Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). \doi{10.1145/1989323.1989399}}
#' @references simmelian: {Nick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). \doi{10.1145/2492517.2492569}}
#' @references jaccard, meetmin, geometric, hypergeometric: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @references localdegree: {Hamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. *Social Network Analysis and Mining, 6*, 22. \doi{10.1007/s13278-016-0332-2}}
#' @references quadrilateral: {Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. *Journal of Graph Algorithms and Applications, 19*, 595-618. \doi{10.7155/jgaa.00370}}
#' @export
#'
#' @examples
#' bb <- backbone_from_unweighted()
backbone_from_unweighted <- function(U,
                                     model,
                                     s,
                                     escore,
                                     normalize,
                                     filter,
                                     symmetrize = TRUE,
                                     umst = FALSE,
                                     narrative = TRUE,
                                     return = "backbone") {

  #ADD OUTLINE HERE

}
