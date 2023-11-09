#' Extract the backbone from a network using a sparsification model
#'
#' @description
#' A generic function to extract the backbone of an undirected, unipartite
#' network using a sparsification model described by a combination of an edge scoring metric, a
#' edge score normalization, and an edge score filter.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsification parameter
#' @param escore string: Method for scoring edges' importance
#' @param normalize string: Method for normalizing edge scores
#' @param filter string: Type of filter to apply
#' @param symmetrize boolean: TRUE if the result should be symmetrized
#' @param umst boolean: TRUE if the backbone should include the union of minimum spanning trees, to ensure connectivity
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `escore` parameter determines how an unweighted edge's importance is calculated.
#' Unless noted below, scores are symmetric and larger values represent more important edges.
#' There are 12 options for assigning an edge's score; when `escore = `
#' * `random`: a random number drawn from a uniform distribution
#' * `betweenness`: edge betweenness
#' * `triangles`: number of triangles that include the edge
#' * `jaccard`: jaccard similarity coefficient of the neighborhoods of an edge's endpoints, or alternatively, triangles normalized by the size of the union of the endpoints neighborhoods
#' * `dice`: dice similarity coefficient of the neighborhoods of an edge's endpoints
#' * `invlogweighted`: inverse log weighted similarity coefficient of the neighborhoods of an edge's endpoints
#' * `quadrangles`: number of quadrangles that include the edge
#' * `quadrilateral embeddedness`: geometric mean normalization of quadrangles
#' * `degree`: degree of neighbor to which an edge is adjacent (asymmetric)
#' * `meetmin`: triangles normalized by the smaller of the endpoints' neighborhoods' sizes
#' * `geometric`: triangles normalized by the product of the endpoints' neighborhoods' sizes
#' * `hypergeometric`: probability of the edge being included at least as many triangles if edges were random, given the size of the endpoints' neighborhoods (smaller is more important)
#'
#' The `normalize` parameter determines whether edge scores are normalized.
#' There are three options; when `normalize = `
#' * `none`: no normalization is performed
#' * `rank`: scores are normalized by neighborhood rank, such that the strongest edge in a node's neighborhood is ranked 1 (asymmetric)
#' * `embeddedness`: scores are normalized using the maximum Jaccard coefficient of the top k-ranked neighbors of each endpoint, for all k
#'
#' Using `escore == "degree"` or `normalize == "rank"` can yield an assymmetric network. When `symmetrize == TRUE`, the network is
#'  symmetrized before applying a filter by defining an edge between i and j if either i->j or i<-j.
#'
#' The `filter` parameter determines how edges are filtered based on their (normalized) edge scores.
#' There are three options; when `filter = `
#' * `threshold`: Edges with scores more important than `s` are retained in the backbone
#' * `proportion`: Specifies the proportion of most important edges to retain in the backbone
#' * `degree`: Retains each node's d^`s` most important edges, where d is the node's degree (requires that `normalize = "rank"`)
#'
#' Specific combinations of `escore`, `normalize`, `filter`, and `umst` correspond to specific sparsification models in the literature, and are available via the following wrapper functions:
#' [sparsify.with.skeleton()], [sparsify.with.gspar()], [sparsify.with.lspar()], [sparsify.with.simmelian()], [sparsify.with.jaccard()], [sparsify.with.meetmin()], [sparsify.with.geometric()], [sparsify.with.hypergeometric()], [sparsify.with.localdegree()], [sparsify.with.quadrilateral()].
#' See the documentation for these wrapper functions for more details and the associated citation.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @export
#'
#' @references {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify(U, s = 0.6, escore = "jaccard", normalize = "rank",
#' filter = "degree", narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify <- function(U, s, escore, normalize, filter, symmetrize = TRUE, umst = FALSE, class = "original", narrative = FALSE) {

  #### Helper Function: Edge score ranking ####
  nhood.rank <- function(x) {
    if (max(x)==0) {return(x)} else {  #Do nothing for isolate nodes
      old <- sort(unique(x))  #Find unique values
      new <- c((length(old)):1)  #Rank them 1 = highest, 2 = second highest, etc
      if (min(old)==0) {new[which(new==max(new))] <- 0}  #If zero was one of the values, rank them as 0
      x <- new[match(x, old)]  #Replace original values with corresponding ranks
      return(x)
    }
  }

  #### Convert supplied object to matrix ####
  G <- tomatrix(U)
  if (G$summary$bipartite==TRUE | G$summary$symmetric==FALSE | G$summary$weighted==TRUE) {stop("G must be an undirected, unweighted, unipartite network")}
  if (class == "original") {class <- G$summary$class}
  attribs <- G$attribs
  original <- G$G  #Original graph
  G <- G$G  #Copy to be manipulated

  #### Sparsification model checks ####
  if (is.null(s)) {stop("A sparsification parameter `s` must be specified")}
  if (escore != "random" & escore != "betweenness" & escore != "triangles" & escore != "jaccard" & escore != "dice" & escore != "invlogsquared" &
      escore != "quadrangles" & escore != "quadrilateral embeddedness" & escore != "degree" & escore != "meetmin" &
      escore != "geometric" & escore != "hypergeometric") {stop("escore must be one of: random, betweenness, triangles, jaccard, dice, invlogsquared, quadrangles, quadrilateral embeddedness, degree, meetmin, geometric, hypergeometric")}
  if (normalize != "none" & normalize != "rank" & normalize != "embeddedness") {stop("normalize must be one of: none, rank, embeddedness")}
  if (filter != "threshold" & filter != "proportion" & filter != "degree") {stop("filter must be one of: threshold, proportion, degree")}
  if (filter == "degree" & normalize != "rank") {stop("The degree filter requires that normalize = \"rank\"")}  #Degree filter assumes edge scores are integer ranks

  #### Compute edge scores, if requested ####
  #Random, from Karger (1994)
  if (escore == "random") {
    G <- G*stats::runif(length(G))  #Assign each edge a random weight
    G[lower.tri(G)] <- t(G)[lower.tri(G)]  #Make symmetric
    }

  #Edge betweenness, from Melancon & Sallaberry (2008)
  if (escore == "betweenness") {
    G <- igraph::graph_from_adjacency_matrix(G,mode="undirected")
    igraph::E(G)$weight <- igraph::edge_betweenness(G, directed = FALSE)
    G <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
    }

  #Number of triangles, from Nick et al. (2013)
  if (escore == "triangles") {
    G <- tcrossprod(G)
    G <- G * original
  }

  #Jaccard coefficient (aka Neighborhood-normalized number of triangles), from Satuluri et al. (2011)
  if (escore == "jaccard") {
    N <- tcrossprod(G)  #Count triangles (union of neighborhoods, numerator of jaccard)
    D <- nrow(G) - tcrossprod((!G)*1)  #Intersection of neighborhoods, denominator of jaccard
    G <- N/D  #Jaccard coefficient
    G[is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * original
  }

  #Dice coefficient
  if (escore == "dice") {
    N <- tcrossprod(G)  #Count triangles
    D <- matrix(1, nrow(G), ncol(G))  #Matrix of sum of degrees
    D[lower.tri(D)] <- utils::combn(rowSums(G), 2, FUN = sum)
    D[upper.tri(D)] <- t(D)[upper.tri(D)]
    G <- (2*N)/D  #Dice coefficient
    G[is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * original
  }

  #Inverse log-weighted
  if (escore == "invlogweighted") {
    G <- igraph::graph_from_adjacency_matrix(G,mode="undirected")
    G <- igraph::similarity(G, mode = "all", method = "invlogweighted", loops = FALSE)
    G <- G * original
  }

  #Number of maximal 4-cliques (i.e., quadrangles), from Nocaj et al. (2015)
  if (escore == "quadrangles" | escore == "quadrilateral embeddedness") {
    G <- igraph::graph_from_adjacency_matrix(G,mode="undirected")
    quads <- matrix(unlist(igraph::cliques(G, min=4, max=4)), nrow = 4) #Value can be replaced to count an edge's number of k-clique
    quads <- as.data.frame(table(data.frame(do.call(rbind,unlist(apply(quads, 2, function(x) utils::combn(sort(x), 2, simplify = FALSE)),recursive = FALSE)))))
    quads <- subset(quads, quads$Freq > 0)
    quads$edgeid <- igraph::get.edge.ids(G, as.numeric(as.vector(unlist(t(quads[,1:2])))))
    igraph::E(G)$weight <- 0
    igraph::E(G)$weight[quads$edge] <- quads$Freq[which(quads$edgeid==quads$edge)]
    G <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
  }

  #Neighborhood-normalized quadrangle count, from Nocaj et al. (2015)
  if (escore == "quadrilateral embeddedness") {  #G already contains the number of quadrangles per edge
    denominator <- sqrt(rowSums(G)%*%t(colSums(G)))
    G <- (G / denominator) * original
    G[is.nan(G)] <- 0  #Fix any divide-by-zero
  }

  #Degree of alter, from Hamann et al. (2016)
  if (escore == "degree") {
    G <- t(rowSums(G)*G)
    G <- G * original
  }

  #Meet/min, from Goldberg & Roth (2003)
  if (escore == "meetmin") {
    N <- tcrossprod(G)  #Shared neighbors
    D <- pmin(G*rowSums(G), t(G*rowSums(G)))  #Minimum of i's and j's degree
    G <- N/D  #Meet-min score
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
  }

  #Geometric, from Goldberg & Roth (2003)
  if (escore == "geometric") {
    N <- tcrossprod(G)^2  #Shared neighbors, squared
    D <- rowSums(G)%*%t(rowSums(G))
    G <- N/D  #Geometric score
    G[is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * original
  }

  #Hypergeometric, from Goldberg & Roth (2003)
  if (escore == "hypergeometric") {
    triangles <- tcrossprod(G)
    G <- outer(1:nrow(G),1:ncol(G), FUN = Vectorize( function(i,j) stats::phyper(triangles[i,j]-1, sum(G[i,])-1, (nrow(G)-2)-(sum(G[i,])-1), sum(G[j,])-1, lower.tail=FALSE) ))
    G <- G * original
  }

  #### Normalize edge scores ####
  #Neighborhood rank, from Satuluri et al. (2011)
  if (normalize == "rank" | normalize == "embeddedness") {for (i in 1:nrow(G)) {G[i,] <- nhood.rank(G[i,])}}  #Rank edges by row (i.e., from perspective of each node)

  #Embeddedness, from Nick et al. (2013) and Nocaj et al. (2015)
  if (normalize == "embeddedness") {  #Scores will already be transformed as neighborhood ranks
    scores <- matrix(0, nrow(G), ncol(G))  #Initialize matrix to hold embeddedness scores
    for (row1 in 1:(nrow(G)-1)) {
      for (row2 in (row1+1):nrow(G)) {  #Loop over each pair of rows
        list1 <- G[row1,]  #Vector of ranked edges for row1
        list2 <- G[row2,]  #Vector of ranked edges for row2

        #Find overlap between neighborhoods using non-parametric variant
        k <- max(list1,list2)
        if (k==0 | ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k))))==0) {  #If jaccard for max(k) is zero, stop
          scores[row1,row2] <- 0
        } else {  #Otherwise, compute jaccard for each k, use maximum
          j <- NULL
          for (k in 1:max(list1,list2)) {j <- c(j, ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k)))))}
          scores[row1,row2] <- max(j)
        }
      }
    }
    G <- scores
    G[lower.tri(G)] <- t(G)[lower.tri(G)]  #Make symmetric
  }

  #### Symmetrize if requested ####
  if (symmetrize) {
    G[lower.tri(G)] <- pmax(G[lower.tri(G)],t(G)[lower.tri(t(G))])
    G[upper.tri(G)] <- t(G)[upper.tri(G)]
  }

  #### Apply filter ####
  #Threshold
  if (filter == "threshold") {
    if (escore != "hypergeometric" & normalize != "rank") {G <- (G >= s)*1}  #Cases where large edge scores are stronger
    if (escore == "hypergeometric" | normalize == "rank") {G <- (G<=s & G!=0)*1}  #Cases where small non-zero edge scores are stronger
    }

  #Proportion
  if (filter == "proportion") {
    scores <- G[lower.tri(G)][which(G[lower.tri(G)]!=0)]  #Vector of non-zero edge scores
    tokeep <- ceiling(s*length(scores))  #Number of edges to keep
    if (escore != "hypergeometric" & normalize != "rank") {  #Cases where large edge scores are stronger
      keep.score <- sort(scores, decreasing = TRUE)[tokeep]  #Value of the tokeep^th edge score, starting from largest value
      G <- (G >= keep.score)*1  #Keep edges with scores at least as large
    }
    if (escore == "hypergeometric" | normalize == "rank") {  #Cases where small edge scores are stronger
      keep.score <- sort(scores, decreasing = FALSE)[tokeep]  #Value of the tokeep^th edge score, starting from smallest value
      G <- (G <= keep.score)*1  #Keep edges with scores at least as small
      G[which(original==0)] <- 0  #But, don't count edges with score = 0, which should be missing
    }
  }

  #Degree exponent, from Satuluri et al. (2011)
  if (filter == "degree") {
    G <- (G <= (floor(rowSums(original)^s)))*1  #Keep edges with scores at least as small as degree^s
    G[which(original==0)] <- 0  #But, don't count edges with score = 0, which should be missing
    }

  #### Add UMST if requested ####
  if (umst) {
    tree <- igraph::graph_from_adjacency_matrix(original, mode = "undirected")  #Convert original to igraph
    tree <- igraph::mst(tree)  #Find the UMST
    tree <- igraph::as_adjacency_matrix(tree, sparse = FALSE)  #Convert back to matrix
    G <- (G | tree)*1  #Include an edge if it is in either the sparsified graph or the tree
  }

  #### Display narrative if requested ####
  model <- ""
  if (escore=="random" & normalize=="none" & filter=="proportion" & umst==FALSE) {model <- "skeleton"}
  if (escore=="jaccard" & normalize=="none" & filter=="proportion" & umst==FALSE) {model <- "gspar"}
  if (escore=="jaccard" & normalize=="rank" & filter=="degree" & umst==FALSE) {model <- "lspar"}
  if (escore=="triangles" & normalize=="embeddedness" & filter=="threshold" & umst==FALSE) {model <- "simmelian"}
  if (escore=="jaccard" & normalize=="none" & filter=="threshold" & umst==FALSE) {model <- "jaccard"}
  if (escore=="meetmin" & normalize=="none" & filter=="threshold" & umst==FALSE) {model <- "meetmin"}
  if (escore=="geometric" & normalize=="none" & filter=="threshold" & umst==FALSE) {model <- "geometric"}
  if (escore=="hypergeometric" & normalize=="none" & filter=="threshold" & umst==FALSE) {model <- "hypergeometric"}
  if (escore=="degree" & normalize=="rank" & filter=="degree" & umst==FALSE) {text <- model <- "degree"}
  if (escore=="quadrilateral embeddedness" & normalize=="embeddedness" & filter=="threshold" & umst==TRUE) {model <- "quadrilateral"}
  if (model=="") {model <- "sparify"}
  reduced_edges <- round((sum(original!=0) - sum(G!=0)) / sum(original!=0),3)*100  #Percent decrease in number of edges
  reduced_nodes <- round((max(sum(rowSums(original)!=0),sum(colSums(original)!=0)) - max(sum(rowSums(G)!=0),sum(colSums(G)!=0))) / max(sum(rowSums(original)!=0),sum(colSums(original)!=0)),3) * 100  #Percent decrease in number of connected nodes
  if (narrative == TRUE) {write.narrative(agents = nrow(original), artifacts = NULL, weighted = FALSE, bipartite = FALSE, symmetric = TRUE,
                                          signed = FALSE, mtc = "none", alpha = NULL, s = s, ut = NULL, lt = NULL, trials = NULL, model = model,
                                          reduced_edges = reduced_edges, reduced_nodes = reduced_nodes)}

  #### Return backbone in desired class ####
  rownames(G) <- rownames(original)  #Restore labels if they were lost
  colnames(G) <- colnames(original)
  backbone <- frommatrix(G, attribs, convert = class)  #Convert to desired class
  return(backbone)
}

#### Wrappers ####
#' Extract Karger's (1999) skeleton backbone
#'
#' @description
#' `sparsify.with.skeleton` is a wrapper for [sparsify()] that extracts the skeleton backbone described by Karger (1999),
#' which preserves a specified proportion of random edges. It is equivalent to `sparsify(escore = "random", normalize = "none", filter = "proportion", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Proportion of edges to retain, 0 < s < 1; smaller values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. *Mathematics of Operations Research, 24*, 383-413. \doi{10.1287/moor.24.2.383}}
#' @export
#'
#' @examples
#' U <- igraph::erdos.renyi.game(60, .5)
#' plot(U) #A dense graph
#' sparse <- sparsify.with.skeleton(U, s = 0.25, narrative = TRUE)
#' plot(sparse) #A sparser graph
sparsify.with.skeleton <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "random", normalize = "none", filter = "proportion", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Satuluri et al's (2011) G-spar backbone
#'
#' @description
#' `sparsify.with.gspar` is a wrapper for [sparsify()] that extracts the G-spar backbone described by Satuluri et al. (2011).
#' It is equivalent to `sparsify(escore = "jaccard", normalize = "none", filter = "proportion", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Proportion of edges to retain, 0 < s < 1; smaller values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). \doi{10.1145/1989323.1989399}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.gspar(U, s = 0.4, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.gspar <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "jaccard", normalize = "none", filter = "proportion", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Satuluri et al's (2011) L-spar backbone
#'
#' @description
#' `sparsify.with.lspar` is a wrapper for [sparsify()] that extracts the L-spar backbone described by Satuluri et al. (2011).
#' It is equivalent to `sparsify(escore = "jaccard", normalize = "rank", filter = "degree", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsification exponent, 0 < s < 1; smaller values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). \doi{10.1145/1989323.1989399}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.lspar(U, s = 0.6, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.lspar <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "jaccard", normalize = "rank", filter = "degree", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Nick et al's (2013) Simmelian backbone
#'
#' @description
#' `sparsify.with.simmelian` is a wrapper for [sparsify()] that extracts the simmelian backbone described by Nick et al. (2013).
#' It is equivalent to `sparsify(escore = "triangles", normalize = "embeddedness", filter = "threshold", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsificiation threshold, 0 < s < 1; larger values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Nick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). \doi{10.1145/2492517.2492569}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.simmelian(U, s = 0.5, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.simmelian <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "triangles", normalize = "embeddedness", filter = "threshold", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Goldberg and Roth's (2003) Jaccard backbone
#'
#' @description
#' `sparsify.with.jaccard` is a wrapper for [sparsify()] that extracts the jaccard backbone described by Goldberg and Roth (2003).
#' It is equivalent to `sparsify(escore = "jaccard", normalize = "none", filter = "threshold", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsificiation threshold, 0 < s < 1; larger values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.jaccard(U, s = 0.3, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.jaccard <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "jaccard", normalize = "none", filter = "threshold", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Goldberg and Roth's (2003) MeetMin backbone
#'
#' @description
#' `sparsify.with.meetmin` is a wrapper for [sparsify()] that extracts the meetmin backbone described by Goldberg and Roth (2003).
#' It is equivalent to `sparsify(escore = "meetmin", normalize = "none", filter = "threshold", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsificiation threshold, 0 < s < 1; larger values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.meetmin(U, s = 0.5, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.meetmin <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "meetmin", normalize = "none", filter = "threshold", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Goldberg and Roth's (2003) Geometric backbone
#'
#' @description
#' `sparsify.with.geometric` is a wrapper for [sparsify()] that extracts the geometric backbone described by Goldberg and Roth (2003).
#' It is equivalent to `sparsify(escore = "geometric", normalize = "none", filter = "threshold", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsificiation threshold, 0 < s < 1; larger values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.geometric(U, s = 0.25, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.geometric <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "geometric", normalize = "none", filter = "threshold", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Goldberg and Roth's (2003) Hypergeometric backbone
#'
#' @description
#' `sparsify.with.hypergeometric` is a wrapper for [sparsify()] that extracts the hypergeometric backbone described by Goldberg and Roth (2003).
#' It is equivalent to `sparsify(escore = "hypergeometric", normalize = "none", filter = "threshold", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsificiation threshold, 0 < s < 1; smaller values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. *Proceedings of the National Academy of Sciences, 100*, 4372-4376. \doi{10.1073/pnas.0735871100}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.hypergeometric(U, s = 0.3, narrative = TRUE)
#' plot(sparse) #Clearly visible communities
sparsify.with.hypergeometric <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "hypergeometric", normalize = "none", filter = "threshold", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Hamann et al.'s (2016) Local Degree backbone
#'
#' @description
#' `sparsify.with.localdegree` is a wrapper for [sparsify()] that extracts the local degree backbone described by Hamann et al. (2016).
#' It is equivalent to `sparsify(escore = "degree", normalize = "rank", filter = "degree", umst = FALSE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsification exponent, 0 < s < 1; smaller values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Hamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. *Social Network Analysis and Mining, 6*, 22. \doi{10.1007/s13278-016-0332-2}}
#' @export
#'
#' @examples
#' U <- igraph::as.undirected(igraph::sample_pa(60, m = 3), mode = "collapse")
#' plot(U) #A hairball
#' sparse <- sparsify.with.localdegree(U, s = 0.3, narrative = TRUE)
#' plot(sparse) #Clearly visible hubs
sparsify.with.localdegree <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "degree", normalize = "rank", filter = "degree", s = s, umst = FALSE, class = class, narrative = narrative)
}

#' Extract Nocaj et al.'s (2015) Quadrilateral Simmelian backbone
#'
#' @description
#' `sparsify.with.quadrilateral` is a wrapper for [sparsify()] that extracts the quadrilateral Simmelian backbone described by Nocaj et al. (2015).
#' It is equivalent to `sparsify(escore = "quadrilateral embeddedness", normalize = "embeddedness", filter = "threshold", umst = TRUE)`.
#'
#' @param U An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param s numeric: Sparsification exponent, 0 < s < 1; larger values yield sparser graphs
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `U`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @return An unweighted, undirected, unipartite graph of class `class`.
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references model: {Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. *Journal of Graph Algorithms and Applications, 19*, 595-618. \doi{10.7155/jgaa.00370}}
#' @export
#'
#' @examples
#' U <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
#' plot(U) #A hairball
#' sparse <- sparsify.with.quadrilateral(U, s = 0.5, narrative = TRUE)
#' plot(sparse) #Clearly visible communities in a connected graph
sparsify.with.quadrilateral <- function(U, s, class = "original", narrative = FALSE) {
  sparsify(U, escore = "quadrilateral embeddedness", normalize = "embeddedness", filter = "threshold", s = s, umst = TRUE, class = class, narrative = narrative)
}
