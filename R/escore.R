#' Compute edge scores in an unweighted network
#'
#' @param A a binary adjacency matrix
#' @param escore string: type of edge score to compute
#'
#' @return A weighted adjacency matrix
#'
#' @noRd
.escore <- function(A, escore){

  W <- A  #Weighted graph to be generated

  #### Random, from Karger (1994) ###
  if (escore == "random") {
    W <- W*stats::runif(length(W))  #Assign each edge a random weight
    W[lower.tri(W)] <- t(W)[lower.tri(W)]  #Make symmetric
    }

  #### Edge betweenness, from Melancon & Sallaberry (2008) ####
  if (escore == "betweenness") {
    W <- igraph::graph_from_adjacency_matrix(W,mode="undirected")
    igraph::E(W)$weight <- igraph::edge_betweenness(W, directed = FALSE)
    W <- igraph::as_adjacency_matrix(W, attr = "weight", sparse = FALSE)
    }

  #### Number of triangles, from Nick et al. (2013) ####
  if (escore == "triangles") {
    W <- tcrossprod(W)
    W <- W * A
  }

  #### Jaccard coefficient (aka Neighborhood-normalized number of triangles), from Satuluri et al. (2011) ####
  if (escore == "jaccard") {
    N <- tcrossprod(W)  #Union of neighborhoods, excluding focal nodes
    D <- nrow(W) - tcrossprod((!W)*1)  #Intersection of neighborhoods
    D <- D - 2  #Exclude focal nodes from denominator
    W <- N/D  #Jaccard coefficient
    W[W==Inf | is.nan(W)] <- 0  #Fix any divide-by-zero
    W <- W * A  #Keep coefficient only for present edges
  }

  #### Dice coefficient ####
  if (escore == "dice") {
    N <- tcrossprod(W)  #Count triangles
    D <- matrix(1, nrow(W), ncol(W))  #Matrix of sum of degrees
    D[lower.tri(D)] <- utils::combn(rowSums(W), 2, FUN = sum)
    D[upper.tri(D)] <- t(D)[upper.tri(D)]
    D <- D - 2  #Exclude focal nodes from denominator
    W <- (2*N)/D  #Dice coefficient
    W[W==Inf | is.nan(W)] <- 0  #Fix any divide-by-zero
    W <- W * A  #Keep coefficient only for present edges
  }

  #### Number of 4-cliques (i.e., quadrangles), from Nocaj et al. (2015) ####
  if (escore == "quadrangles" | escore == "quadrilateral") {
    W <- igraph::graph_from_adjacency_matrix(W,mode="undirected")
    quads <- matrix(unlist(igraph::cliques(W, min=4, max=4)), nrow = 4) #Value can be replaced to count an edge's number of k-clique
    quads <- as.data.frame(table(data.frame(do.call(rbind,unlist(apply(quads, 2, function(x) utils::combn(sort(x), 2, simplify = FALSE)),recursive = FALSE)))))
    quads <- subset(quads, quads$Freq > 0)
    quads$edgeid <- igraph::get.edge.ids(W, as.numeric(as.vector(unlist(t(quads[,1:2])))))
    igraph::E(W)$weight <- 0
    igraph::E(W)$weight[quads$edge] <- quads$Freq[which(quads$edgeid==quads$edge)]
    W <- igraph::as_adjacency_matrix(W, attr = "weight", sparse = FALSE)
  }

  #### Neighborhood-normalized quadrangle count, from Nocaj et al. (2015) ####
  if (escore == "quadrilateral") {  #W already contains the number of quadrangles per edge
    denominator <- sqrt(rowSums(W)%*%t(colSums(W)))
    W <- (W / denominator) * A
    W[W==Inf | is.nan(W)] <- 0  #Fix any divide-by-zero
  }

  #### Degree of alter, from Hamann et al. (2016) ####
  if (escore == "degree") {
    W <- t(rowSums(W)*W)
    W <- W * A
  }

  #### Meet/min, from Goldberg & Roth (2003) ####
  if (escore == "meetmin") {
    N <- tcrossprod(W)  #Shared neighbors
    D <- pmin(W*rowSums(W), t(W*rowSums(W)))  #Minimum of i's and j's degree
    W <- N/D  #Meet-min score
    W[W==Inf | is.nan(W)] <- 0  #Fix any divide-by-zero
  }

  #### Geometric, from Goldberg & Roth (2003) ####
  if (escore == "geometric") {
    N <- tcrossprod(W)^2  #Shared neighbors, squared
    D <- rowSums(W)%*%t(rowSums(W))
    W <- N/D  #Geometric score
    W[W==Inf | is.nan(W)] <- 0  #Fix any divide-by-zero
    W <- W * A
  }

  #### Hypergeometric, from Goldberg & Roth (2003) ####
  if (escore == "hypergeometric") {
    triangles <- tcrossprod(W)
    W <- outer(1:nrow(W),1:ncol(W), FUN = Vectorize( function(i,j) stats::phyper(triangles[i,j]-1, sum(W[i,])-1, (nrow(W)-2)-(sum(W[i,])-1), sum(W[j,])-1, lower.tail=FALSE) ))
    W <- (1-W) * A  #Reverse-score so that larger weights are assigned to edges more worth keeping
  }

  return(W)
}
