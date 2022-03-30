#' Generates a bipartite network with given row and column degree distributions
#'
#' `bipartite.from.distribution` returns a bipartite graph, as an object of the requested class,
#' with row and column degree distributions that approximately follow beta distributions with given parameters.
#'
#' @param R integer: number of rows
#' @param C integer: number of columns
#' @param P numeric: probability of an edge
#' @param rowdist vector length 2: Row degrees will approximately follow a Beta(a,b) distribution
#' @param coldist vector length 2: Column degrees will approximately follow a Beta(a,b) distribution
#' @param class string: the class of the returned backbone graph, one of c("matrix", "Matrix", "sparseMatrix", "igraph", "network")
#'
#' @export
#'
#' @examples
#' B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1,
#' rowdist = c(1,1), coldist = c(1,1))  #Uniform
#' B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1,
#' rowdist = c(1,10), coldist = c(1,10))  #Right-tailed
#' B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1,
#' rowdist = c(10,1), coldist = c(10,1))  #Left-tailed
#' B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1,
#' rowdist = c(10,10), coldist = c(10,10))  #Normal
#' B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1,
#' rowdist = c(10000,10000), coldist = c(10000,10000))  #Constant
bipartite.from.distribution <- function(R, C, P, rowdist = c(1,1), coldist = c(1,1), class="matrix") {

  # Vector of N integers with sum S and approximately distributed as beta(a,b) ####
  integers <- function(N, S, a, b){
    for (i in 1:500000) {  #Attempt up to 500000 times
      if (a < b) {  #For right-tailed distributions
        rand <- rep(1,N)  #Start with vector of 1s, to set a floor
        plus <- stats::rbeta(N,a,b)  #Distribution of amount to add to the floor
        Sstar <- S-N  #Total value that can be added to the floor
        rand <- rand + round(plus * (Sstar / sum(plus)))  #Add extra value to floor
        }
      if (a >= b) {  #For symmetric and left-tailed distributions
        rand <- stats::rbeta(N,a,b)  #Distribution of values
        rand <- round(rand * (S / sum(rand)))  #Normalize to match desired total
        }
      if (sum(rand) == S & min(rand) != 0) {break}  #Stop when an allowable vector is generated
      }
    if (sum(rand) == S | min(rand) != 0) {return(rand)} else {stop("These degree distributions are not possible")}
    }

  #Parameter check
  if (!is.numeric(R) | !is.numeric(C) | !is.numeric(P)) {stop("R, C, and P must be numeric")}
  if (R<0 | C<0 | R%%1!=0 | C%%1!=0) {stop("R and C must be positive integers")}
  if (P<0 | P>1) {stop("P must be between 0 and 1")}

  # Generate bipartite
  ones <- round(R * C * P)  #Number of 1s in matrix, given dimensions and density
  for (i in 1:10000) {
    R <- integers(R,ones,rowdist[1],rowdist[2])  #Create a vector of agent degrees; if not possible, error
    C <- integers(C,ones,coldist[1],coldist[2])  #Create a vector of artifact degrees; if not possible, error
    B <- bipartite.from.sequence(R,C)  #Create bipartite
    if (all.equal(R,rowSums(B)) == TRUE | all.equal(C,colSums(B)) == TRUE) {break}  #Repeat until conditions are achieved
  }

  #Verify and return
  if (all.equal(R,rowSums(B)) == TRUE | all.equal(C,colSums(B)) == TRUE) {
    if (class == "Matrix"){B <- Matrix::Matrix(B)}
    if (class == "sparseMatrix"){B <- Matrix::Matrix(B, sparse = TRUE)}
    if (class == "network"){B <- network::network(B, bipartite = TRUE)}
    if (class == "igraph"){B <- igraph::graph_from_incidence_matrix(B)}
    return(B)
  } else {stop("A bipartite network with these degree distributions is not possible")}
}
