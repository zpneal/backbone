#' Compute edge scores in an unweighted network
#'
#' @param A a binary adjacency matrix
#' @param escore string: type of edge score to compute
#'
#' @return A weighted adjacency matrix
#'
#' @noRd
.escore <- function(A, escore){

  G <- A  #Graph to be manipulated

  #### Random, from Karger (1994) ###
  if (escore == "random") {
    G <- G*stats::runif(length(G))  #Assign each edge a random weight
    G[lower.tri(G)] <- t(G)[lower.tri(G)]  #Make symmetric
    }

  #### Edge betweenness, from Melancon & Sallaberry (2008) ####
  if (escore == "betweenness") {
    G <- igraph::graph_from_adjacency_matrix(G,mode="undirected")
    igraph::E(G)$weight <- igraph::edge_betweenness(G, directed = FALSE)
    G <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
    }

  #### Number of triangles, from Nick et al. (2013) ####
  if (escore == "triangles") {
    G <- tcrossprod(G)
    G <- G * A
  }

  #### Jaccard coefficient (aka Neighborhood-normalized number of triangles), from Satuluri et al. (2011) ####
  if (escore == "jaccard") {
    N <- tcrossprod(G)  #Union of neighborhoods, excluding focal nodes
    D <- nrow(G) - tcrossprod((!G)*1)  #Intersection of neighborhoods
    D <- D - 2  #Exclude focal nodes from denominator
    G <- N/D  #Jaccard coefficient
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * A  #Keep coefficient only for present edges
  }

  #### Dice coefficient ####
  if (escore == "dice") {
    N <- tcrossprod(G)  #Count triangles
    D <- matrix(1, nrow(G), ncol(G))  #Matrix of sum of degrees
    D[lower.tri(D)] <- utils::combn(rowSums(G), 2, FUN = sum)
    D[upper.tri(D)] <- t(D)[upper.tri(D)]
    D <- D - 2  #Exclude focal nodes from denominator
    G <- (2*N)/D  #Dice coefficient
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * A  #Keep coefficient only for present edges
  }

  #### Number of 4-cliques (i.e., quadrangles), from Nocaj et al. (2015) ####
  if (escore == "quadrangles" | escore == "quadrilateral") {
    G <- igraph::graph_from_adjacency_matrix(G,mode="undirected")
    quads <- matrix(unlist(igraph::cliques(G, min=4, max=4)), nrow = 4) #Value can be replaced to count an edge's number of k-clique
    quads <- as.data.frame(table(data.frame(do.call(rbind,unlist(apply(quads, 2, function(x) utils::combn(sort(x), 2, simplify = FALSE)),recursive = FALSE)))))
    quads <- subset(quads, quads$Freq > 0)
    quads$edgeid <- igraph::get.edge.ids(G, as.numeric(as.vector(unlist(t(quads[,1:2])))))
    igraph::E(G)$weight <- 0
    igraph::E(G)$weight[quads$edge] <- quads$Freq[which(quads$edgeid==quads$edge)]
    G <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
  }

  #### Neighborhood-normalized quadrangle count, from Nocaj et al. (2015) ####
  if (escore == "quadrilateral") {  #G already contains the number of quadrangles per edge
    denominator <- sqrt(rowSums(G)%*%t(colSums(G)))
    G <- (G / denominator) * A
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
  }

  #### Degree of alter, from Hamann et al. (2016) ####
  if (escore == "degree") {
    G <- t(rowSums(G)*G)
    G <- G * A
  }

  #### Meet/min, from Goldberg & Roth (2003) ####
  if (escore == "meetmin") {
    N <- tcrossprod(G)  #Shared neighbors
    D <- pmin(G*rowSums(G), t(G*rowSums(G)))  #Minimum of i's and j's degree
    G <- N/D  #Meet-min score
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
  }

  #### Geometric, from Goldberg & Roth (2003) ####
  if (escore == "geometric") {
    N <- tcrossprod(G)^2  #Shared neighbors, squared
    D <- rowSums(G)%*%t(rowSums(G))
    G <- N/D  #Geometric score
    G[G==Inf | is.nan(G)] <- 0  #Fix any divide-by-zero
    G <- G * A
  }

  #### Hypergeometric, from Goldberg & Roth (2003) ####
  if (escore == "hypergeometric") {
    triangles <- tcrossprod(G)
    G <- outer(1:nrow(G),1:ncol(G), FUN = Vectorize( function(i,j) stats::phyper(triangles[i,j]-1, sum(G[i,])-1, (nrow(G)-2)-(sum(G[i,])-1), sum(G[j,])-1, lower.tail=FALSE) ))
    G <- G * A
  }

  return(G)
}
