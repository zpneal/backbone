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
    W <- W * A
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
    W <- W * A
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

#' Normalize edge scores
#'
#' @param W a weighted adjacency matrix
#' @param normalize string: type of normalization
#'
#' @return A weighted adjacency matrix
#'
#' @noRd
.normalize <- function(W, normalize) {

  original <- W

  #### Neighborhood rank, from Satuluri et al. (2011) ####
  if (normalize == "rank" | normalize == "embeddedness") {
    for (i in 1:nrow(W)) {  #For each row (i.e., from the perspective of each node)
      x <- W[i,]  #Vector of values from this row
      old <- sort(unique(x))  #Find unique values
      new <- c((length(old)):1)  #Rank them 1 = highest, 2 = second highest, etc
      if (min(old)==0) {new[which(new==max(new))] <- 0}  #If zero was one of the values, rank them as 0
      x <- new[match(x, old)]  #Replace original values with corresponding ranks
      W[i,] <- x  #Put ranks into row
    }
  }

  #### Embeddedness, from Nick et al. (2013) ####
  if (normalize == "embeddedness") {  #Scores will already be transformed as neighborhood ranks
    scores <- matrix(0, nrow(W), ncol(W))  #Initialize matrix to hold embeddedness scores
    for (row1 in 1:(nrow(W)-1)) {
      for (row2 in (row1+1):nrow(W)) {  #Loop over each pair of rows
        list1 <- W[row1,-c(row1,row2)]  #Vector of ranked edges for row1, excluding row1 and row2
        list2 <- W[row2,-c(row1,row2)]  #Vector of ranked edges for row2, excluding row1 and row2

        #Find overlap between neighborhoods using non-parametric variant
        k <- max(list1,list2)
        if (k==0 | ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k))))==0) {  #If jaccard for max(k) is zero, stop
          scores[row1,row2] <- 0
        } else {  #Otherwise, compute jaccard for each k, use maximum
          j <- NULL
          for (k in 1:max(list1,list2)) {j <- c(j, ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k)))))}
          scores[row1,row2] <- max(j, na.rm = TRUE)
        }
      }
    }
    W[lower.tri(W)] <- t(W)[lower.tri(W)]  #Fill in rest of matrix
    W <- scores * (original!=0)*1  #Only keep scores for present edges
  }

  return(W)
}

#' Filter based on edge scores
#'
#' @param W a weighted adjacency matrix
#' @param filter string: filter method
#' @param parameter numeric: filtering parameter
#'
#' @return A binary adjacency matrix representing the backbone
#'
#' @noRd
.filter <- function(W, filter, parameter){

  #### Threshold ####
  #Keep edges with weights greater than `parameter`; Depends on edge weight scaling, but larger values keep fewer edges
  if (filter == "threshold") {B <- (W > parameter)*1}

  #### Proportion ####
  #Keep strongest `parameter` proportion of edges; 0 = keep 0% of edges, 1 = keep 100% of edges
  if (filter == "proportion") {
    scores <- W[which(W!=0)]  #Vector of non-zero edge scores
    B <- (W >= stats::quantile(scores, probs = (1 - parameter)))*1
  }

  #### Degree exponent, from Satuluri et al. (2011) ####
  #Keep edges with neighborhood rank scores at least as small as degree^`parameter`; 0 = keep one edge per node, 1 = keep all edges per node
  if (filter == "degree") {
    A <- (W != 0)*1  #Adjacency matrix
    B <- (W <= (floor(rowSums(A)^parameter)) & W!=0)*1
  }

  #### Weighted backbone models ####
  #Parameter functions as alpha: 0 = keep no edges, 1 = keep all edges
  if (filter %in% c("disparity", "lans", "mlf")) {B <- backbone_from_weighted(W, model = filter, alpha = parameter)}

  return(B)

}
