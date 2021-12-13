#' Generates a bipartite network with given edge probability
#'
#' `bipartite.from.probability` returns a bipartite graph, as an object of the requested class,
#'  in which each edge has a given probability and where no node is an isolate or maximally connected.
#'
#' @param R integer: number of rows
#' @param C integer: number of columns
#' @param P numeric: probability of an edge; if P = 0 a probability will be chosen randomly
#' @param class string: the class of the returned backbone graph, one of c("matrix", "Matrix", "sparseMatrix", "igraph", "network").
#'
#' @export
#'
#' @examples
#' B <- bipartite.from.probability(R = 10, C = 10)
#' B <- bipartite.from.probability(R = 10, C = 10, P = .5)
#' B <- bipartite.from.probability(R = 10, C = 10, P = .5, class = "igraph")
bipartite.from.probability <- function(R, C, P=0, class="matrix") {

  #Parameter check
  if (!is.numeric(R) | !is.numeric(C) | !is.numeric(P)) {stop("R, C, and P must be numeric")}
  if (R<0 | C<0 | R%%1!=0 | C%%1!=0) {stop("R and C must be positive integers")}
  if (P<0 | P>1) {stop("P must be between 0 and 1")}

  #Compute and check allowable probability
  minP <- max(R,C) / (R*C)  #Minimum probability
  maxP <- ((R*C) - max(R,C)) / (R*C)  #Maximum probability
  if (P == 0) {P <- sample(seq(from=minP, to=maxP, by=0.001),1)}  #If unspecified, pick a probability
  if (P < minP | P > maxP) {stop("P is outside the allowed range for this size bipartite network")}

  #Create an empty bipartite matrix
  B <- matrix(stats::rbinom(R*C,1,0),R,C)

  #Randomly create bipartite matrices until no rows or columns contain only 0s or 1s
  while (min(rowSums(B))==0 | max(rowSums(B))==C | min(colSums(B))==0 | max(colSums(B))==R ) {B <- matrix(stats::rbinom(R*C,1,P),R,C)}

  #Convert to desired class and return
  if (class == "Matrix"){B <- Matrix::Matrix(B)}
  if (class == "sparseMatrix"){B <- Matrix::Matrix(B, sparse = TRUE)}
  if (class == "network"){B <- network::network(B, bipartite = TRUE)}
  if (class == "igraph"){B <- igraph::graph_from_incidence_matrix(B)}
  return(B)
}
