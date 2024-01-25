#' Randomize a binary matrix using the fastball algorithm
#'
#' `fastball` randomizes a binary matrix, preserving the row and column sums
#'
#' @param M matrix: a binary matrix (see details)
#' @param trades integer: number of trades; the default is 5R trades (approx. mixing time)
#'
#' @return
#' matrix: A random binary matrix with same row sums and column sums as M.
#'
#' @details
#' Given a matrix `M`, `fastball` randomly samples a new matrix from the space of all matrices with the same row and column sums as `M`.
#'
#' @references fastball: {Godard, Karl and Neal, Zachary P. 2022. fastball: A fast algorithm to sample bipartite graphs with fixed degree sequences. *Journal of Complex Networks* \doi{10.1093/comnet/cnac049}}
#'
#' @export
#' @examples
#' M <- matrix(rbinom(200,1,0.5),10,20)  #A random 10x20 binary matrix
#' Mrand <- fastball(M)  #Random matrix with same row and column sums
fastball <- function(M, trades = 5 * nrow(M)) {
  if (methods::is(M, "matrix")) {
    
    #Convert matrix to adjacency list
    if (as.numeric(R.Version()$major)>=4 & as.numeric(R.Version()$minor)>=1) {
      L <- apply(M==1, 1, which, simplify = FALSE)  #Slightly faster, requires R 4.1.0
    } else {
      L <- lapply(asplit(M == 1, 1), which)  #Slightly slower, works for earlier version of R
      }
    
    Lrand <- fastball_cpp(L, trades)
    Mrand <- matrix(0,nrow(M),ncol(M))
    for (row in 1:nrow(Mrand)) {Mrand[row,Lrand[[row]]] <- 1L}
    return(Mrand)
  } else {return(fastball_cpp(M, 5 * length(M)))}
}