#' Generates a bipartite graph from row and column degree sequences
#' 
#' `bipartite.from.sequence` returns a bipartite graph, as an object of the requested class,
#'  that has the given row and column degree sequences.
#'
#' @param R numeric vector: requested row degree sequence of positive integers
#' @param C numeric vector: requested column degree sequence of positive integers
#' @param class string: the class of the returned backbone graph, one of c("matrix", "Matrix", "sparseMatrix", "igraph", "network")
#' 
#' @return
#' @export
#'
#' @examples
#' B <- bipartite.from.sequence(R = c(1,1,2), C = c(1,1,2))
#' B <- bipartite.from.sequence(R = c(1,1,2), C = c(1,1,2), class = "igraph")
#' B <- bipartite.from.sequence(R = c(1,1,2), C = c(1,1,2), class = "network")
bipartite.from.sequence <- function(R,C,class="matrix"){  
  #Replacement for sample() function so that if length(x)=1, it simply returns x
  newsample <- function(x) {if (length(x) <= 1) {return(x)} else {return(sample(x,1,replace = FALSE, prob = NULL))}}
  
  #Parameter check
  if (sum(R)!=sum(C)) {stop("sum(R) must equal sum(C)")}
  if (!is.numeric(R) | !is.numeric(C)) {stop("R and C must be numeric")}
  if (any(R%%1!=0) | (any(C%%1!=0)) | any(R<1) | any(C<1)) {stop("R and C must only contain positive integers")}
  
  #Set initial values
  r <- length(R) #number of rows
  row_s <- R #row sum vector to change
  col_s <- C #col sum vector to change
  B <- matrix(0, nrow = length(R), ncol = length(C)) #matrix of all zeros to change
  
  #Start inserting 1's to matrix
  for (counter in 1:length(C))
  repeat {
    #Loop stopping conditions
    if ((sum(row_s) == 0) & (sum(col_s) == 0)) {break}
    #Choose a random column index that corresponds to a nonzero entry
    col_index <- newsample(which(col_s > 0))
    #Find the column sum of that column
    c_sum <- col_s[col_index]
    #Rank the entries of the row sum sum vector row_s, breaking ties at random
    ranks <- rank(row_s, ties.method = "random")
    #Select c_sum largest entries, find their indices
    indices <- which(ranks > (r-c_sum))
    #Replace entries in matrix with a 1 if row index in indices, col index as chosen
    B[indices, col_index] <- 1
    #Update row and col sums
    row_s <- (R - rowSums(B))
    col_s <- (C - colSums(B))
  } #end while loop
  
  if ((sum(row_s) == 0) & (sum(col_s) == 0)) {
    B <- curveball(B)
    if (class == "Matrix"){B <- Matrix::Matrix(B)}
    if (class == "sparseMatrix"){B <- Matrix::Matrix(B, sparse = TRUE)}
    if (class == "network"){B <- network::network(B, bipartite = TRUE)}
    if (class == "igraph"){B <- igraph::graph_from_incidence_matrix(B)}
    return(B)  
  } else {stop("A bipartite matrix with these degree sequences does not exist.")}
  
}

