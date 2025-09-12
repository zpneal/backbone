#' Extract the backbone from an unweighted, undirected network
#'
#' \code{backbone_from_unweighted()} extracts the unweighted backbone from an unweighted, undirected network
#'
#' @param U An unweighted, undirected network as an adjacency matrix or \link[Matrix]{Matrix}, or an unweighted unipartite \link[igraph]{igraph} object
#' @param model string: backbone model
#' @param parameter real: filtering parameter
#' @param escore string: Method for scoring edges' importance
#' @param normalize string: Method for normalizing edge scores
#' @param filter string: Type of filter to apply
#' @param umst boolean: TRUE if the backbone should include the union of maximum spanning trees, to ensure connectivity
#' @param narrative boolean: display suggested text & citations
#' @param return string: return either only the \code{"backbone"} or \code{"everything"}
#'
#' @details
#' The \code{backbone_from_unweighted} function extracts the backbone from an unweighted unipartite network. The backbone is an
#' unweighted unipartite network that contains only edges preserved by a backbone model.
#'
#' The following backbone models are available using the \code{model} parameter:
#' * \code{skeleton} - Skeleton backbone (Karger, 1999) 
#' * \code{gspar} - Global Sparsification (Satuluri et al., 2011) 
#' * \code{lspar} - Local Sparsification (Satuluri et al., 2011)
#' * \code{simmelian} - Simmelian backbone (Nick et al., 2013) 
#' * \code{jaccard} - Jaccard backbone (Goldberg and Roth, 2003) 
#' * \code{meetmin} - MeetMin backbone (Goldberg and Roth, 2003) 
#' * \code{geometric} - Geometric backbone (Goldberg and Roth, 2003) 
#' * \code{hyper} - Hypergeometric backbone, (Goldberg and Roth, 2003) 
#' * \code{degree} - Local Degree backbone (Hamann et al, 2016) 
#' * \code{quadrilateral} - Quadrilateral Simmelian backbone (Nocaj et al, 2015) 
#' * \code{custom} - A custom backbone model specified by \code{escore}, \code{normalize}, \code{filter}, and \code{umst}
#'
#' The \code{escore} parameter determines how an unweighted edge's importance is calculated.
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
#' * \code{hypergeometric}: probability of the edge being included at least as many triangles if edges were random, given the size of the endpoints' neighborhoods (inverted, so that larger is more important)
#'
#' The \code{normalize} parameter determines whether edge scores are normalized.
#' * \code{none}: no normalization is performed
#' * \code{rank}: scores are normalized by neighborhood rank, such that the strongest edge in a node's neighborhood is ranked 1 (requires that \code{filter = degree})
#' * \code{embeddedness}: scores are normalized using the maximum Jaccard coefficient of the top k-ranked neighbors of each endpoint, for all k
#'
#' The \code{filter} parameter determines how edges are filtered based on their (normalized) edge scores.
#' * \code{threshold}: Edges with scores > `parameter` are retained in the backbone
#' * \code{proportion}: Specifies the approximate proportion of edges to retain in the backbone
#' * \code{degree}: Retains each node's d^`parameter` most important edges, where d is the node's degree (requires that \code{normalize = "rank"})
#' * \code{disparity}: Applies the disparity filter using [backbone_from_weighted()]
#' * \code{lans}: Applies locally adaptive network sparsification using [backbone_from_weighted()]
#' * \code{mlf}: Applies the marginal likelihood filter using [backbone_from_weighted()]
#'
#' @return If \code{return = "backbone"}, a backbone in the same class as \code{U}. If \code{return = "everything"}, then the backbone
#' is returned as an element in a list that also includes the original unweighted network, a narrative description, and the original
#' function call.
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references skeleton: {Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. *Mathematics of Operations Research, 24*, 383-413. \doi{10.1287/moor.24.2.383}}
#' @references gspar and lspar: {Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). \doi{10.1145/1989323.1989399}}
#' @references simmelian: {Nick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). \doi{10.1145/2492517.2492569}}
#' @references jaccard, meetmin, geometric, hyper: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @references degree: {Hamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. *Social Network Analysis and Mining, 6*, 22. \doi{10.1007/s13278-016-0332-2}}
#' @references quadrilateral: {Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. *Journal of Graph Algorithms and Applications, 19*, 595-618. \doi{10.7155/jgaa.00370}}
#' @export
#'
#' @examples
#' #A dense, unweighted network with three embedded communities
#' U <- igraph::sample_sbm(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U)  #Communities are not obvious
#'
#' #Extract backbone using the built-in "Local Sparsification" model
#' bb <- backbone_from_unweighted(U, model = "lspar", parameter = 0.5)
#' plot(bb)  #Communities are clearly visible
#'
#' #Extract backbone using local sparification, but explicitly specifying the model steps
#' bb <- backbone_from_unweighted(U, model = "custom", escore = "jaccard",
#'                                normalize = "rank", filter = "degree",
#'                                umst = FALSE, parameter = 0.5)
#' plot(bb)  #Communities are clearly visible
backbone_from_unweighted <- function(U,
                                     model = "lspar",
                                     parameter = 0.5,
                                     escore,
                                     normalize,
                                     filter,
                                     umst,
                                     narrative = FALSE,
                                     return = "backbone") {

  call <- match.call()

  #### Check parameters ####
  #All models
  if (!(model %in% c("custom", "skeleton", "gspar", "lspar", "simmelian", "jaccard", "meetmin", "geometric", "hyper", "degree", "quadrilateral"))) {stop("`model` must be one of: \"custom\", \"skeleton\", \"gspar\", \"lspar\", \"simmelian\", \"jaccard\", \"meetmin\", \"geometric\", \"hyper\", \"degree\", \"quadrilateral\"")}
  if (!is.numeric(parameter)) {stop("`parameter` must be a numeric value")}
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
    if (normalize=="rank" & filter!="degree") {stop("Using normalize=\"rank\" requires that filter=\"degree\"")}
    if (normalize!="rank" & filter=="degree") {stop("Using filter=\"degree\" requires that normalize=\"rank\"")}
  }

  #### Check and format input ####
  #Check that input is a weighted adjacency matrix or weighted unipartite igraph
  if (!methods::is(U,"matrix") & !methods::is(U,"Matrix") & !methods::is(U,"igraph")) {stop("`U` must be an adjacency matrix or Matrix, or an igraph object")}

  if (methods::is(U,"matrix")) {
    if (dim(as.matrix(U))[1] != dim(as.matrix(U))[2]) {stop("`U` must be a symmetric adjacency matrix")}
    if (!all(as.matrix(U) %in% c(0,1))) {stop("The entries of `U` must be either 0 or 1")}
    if (!isSymmetric(as.matrix(U))) {stop("`U` must be a symmetric adjacency matrix")}
  }

  if (methods::is(U,"igraph")) {
    if (igraph::is_bipartite(U)) {stop("`U` must be an undirected unipartite igraph object")}
    if (igraph::is_directed(U)) {stop("`U` must be an undirected unipartite igraph object")}
    if ("weight" %in% igraph::edge_attr_names(U)) {stop("An edge weight attribute is present in `U`, but will be ignored")}
  }

  #Convert input to adjacency matrix
  if (methods::is(U,"matrix")) {A <- U}  #matrix --> matrix
  if (methods::is(U,"Matrix")) {A <- as.matrix(U)}  #Matrix --> matrix
  if (methods::is(U,"igraph")) {A <- igraph::as_adjacency_matrix(U, names = FALSE, sparse = FALSE)}

  #### Compute edge scores ####
  G <- .escore(A, escore = escore)

  #### Apply edge score normalization ####
  G <- .normalize(G, normalize = normalize)

  #### Apply filter ####
  backbone <- .filter(G, filter = filter, parameter = parameter)

  #### Symmetrize ####
  backbone <- pmax(backbone, t(backbone))

  #### Add UMST ####
  if (umst) {
    tree <- igraph::graph_from_adjacency_matrix(G, mode = "max", weighted = TRUE)  #Convert weighted matrix to undirected igraph
    if (normalize!="rank") {igraph::E(tree)$weight <- igraph::E(tree)$weight*-1}  #If not using rank normalization, reverse-score weights so that mst() returns *maximum* spanning tree
    tree <- igraph::mst(tree)  #Find the (union of) maximum spanning trees
    tree <- igraph::as_adjacency_matrix(tree, sparse = FALSE)  #Convert tree to matrix
    backbone <- (backbone | tree)*1  #Include an edge if it is in either the backbone or tree
  }

  #### Construct narrative ####
  # First sentence (descriptive)
  text <- paste0("We used the backbone package for R (v", utils::packageVersion("backbone"), "; Neal, 2025) to extract the unweighted backbone of an unweighted network containing ", nrow(A), " nodes.")

  # Second sentence (model and outcome)
  if (model == "skeleton") {desc <- "Karger's (1999) Skeleton backbone"}
  if (model == "gspar") {desc <- "Satuluri et al's (2011) Global Sparsification backbone model"}
  if (model == "lspar") {desc <- "Satuluri et al's (2011) Local Sparsification backbone model"}
  if (model == "simmelian") {desc <- "Nick et al's (2013) Simmelian backbone model"}
  if (model == "jaccard") {desc <- "Goldberg and Roth's (2003) Jaccard backbone model"}
  if (model == "meetmin") {desc <- "Goldberg and Roth's (2003) MeetMin backbone model"}
  if (model == "geometric") {desc <- "Goldberg and Roth's (2003) Geometric backbone model"}
  if (model == "hyper") {desc <- "Goldberg and Roth's (2003) Hypergeometric backbone model"}
  if (model == "degree") {desc <- "Hamann et al.'s (2016) Local Degree backbone model"}
  if (model == "quadrilateral") {desc <- "Nocaj et al.'s (2015) Quadrilateral Simmelian backbone model"}
  if (model == "custom") {desc <- "a custom backbone model specification"}

  old <- sum(A!=0, na.rm=TRUE)  #Number of edges in original network
  new <- sum(backbone!=0)  #Number of edges in backbone
  reduced_edges <- round(((old - new) / old)*100,2)

  text <- paste0(text, " Edges were selected for retention in the backbone using ", desc, " (filtering parameter = ", parameter,"), which reduced the number of edges by ", reduced_edges, "%.")

  #References
  text <- paste0(text, "\n\nNeal, Z. P. 2025. backbone: An R Package to Extract Network Backbones. CRAN. https://doi.org/10.32614/CRAN.package.backbone")
  if (model == "skeleton") {text <- paste0(text, "\n\nKarger, D. R. (1999). Random sampling in cut, flow, and network design problems. Mathematics of Operations Research, 24, 383-413. https://doi/org/10.1287/moor.24.2.383")}
  if (model == "gspar" | model == "lspar") {text <- paste0(text, "\n\nSatuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). https://doi.org/10.1145/1989323.1989399")}
  if (model == "simmelian") {text <- paste0(text, "\n\nNick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). https://doi.org/10.1145/2492517.2492569")}
  if (model == "jaccard" | model == "meetmin" | model == "geometric" | model == "hyper") {text <- paste0(text, "\n\nGoldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. Proceedings of the National Academy of Sciences, 100, 4372-4376. https://doi.org/10.1073/pnas.0735871100")}
  if (model == "degree") {text <- paste0(text, "\n\nHamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. Social Network Analysis and Mining, 6, 22. https://doi.org/10.1007/s13278-016-0332-2")}
  if (model == "quadrilateral") {text <- paste0(text, "\n\nNocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. Journal of Graph Algorithms and Applications, 19, 595-618. https://doi.org/10.7155/jgaa.00370")}

  #### Display narrative ####
  if (narrative) {message(text)}

  #### Prepare backbone ####
  if (methods::is(U,"matrix")) {
    rownames(backbone) <- rownames(U)
    colnames(backbone) <- rownames(U)
  }

  if (methods::is(U,"Matrix")) {
    rownames(backbone) <- rownames(U)
    colnames(backbone) <- rownames(U)
    backbone <- Matrix::Matrix(backbone)
  }

  if (methods::is(U,"igraph")) {
    temp <- U  #Placeholder for backbone
    temp <- igraph::set_edge_attr(temp, "keep", value = backbone[igraph::as_edgelist(temp, names = FALSE)])  #Insert edge retention marker as attribute
    temp <- igraph::delete_edges(temp, which(igraph::E(temp)$keep==0))  #Delete any edges that should not be retained
    temp <- igraph::delete_edge_attr(temp, "keep")  #Delete edge returntion marker
    backbone <- temp
  }

  #### Return ####
  if (return == "backbone") {return(backbone)}
  if (return == "everything") {return(list(original = U, backbone = backbone, narrative = text, call = call))}
}
