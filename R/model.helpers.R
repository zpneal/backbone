#' Convert graph object to adjacency matrix
#'
#' @param graph matrix, sparse matrix, \link{igraph}, edgelist, or \link[network]{network} object
#' @param convert class to convert to, one of "matrix", "sparseMatrix", "igraph", "edgelist", or "network"
#' @param extract Boolean, TRUE if using function within backbone.extract, FALSE if not.
#' @details An object is considered an edgelist if it is (1) a matrix or sparse matrix, and (2) has only two columns.
#'     Each column is understood as a bipartite set, with edges only going between members of column 1 and members of column 2.
#' @return list(class, adjacency), a list containing the class of parameter graph, and the adjacency matrix of the graph
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{davis.sp <- as(davis, "sparseMatrix")}
#' \dontrun{davis.graph <- igraph::graph.incidence(davis)}
#' \dontrun{davis.nw <- network::network(davis, ignore.eval = FALSE,
#'     names.eval = "weight", loops = TRUE)}
#' \dontrun{backbone:::class.convert(davis, "matrix")}
#' \dontrun{backbone:::class.convert(davis.sp, "matrix")}
#' \dontrun{backbone:::class.convert(davis.graph, "matrix")}
#' \dontrun{backbone:::class.convert(davis.nw, "matrix")}
#' \dontrun{bb.sdsm <- sdsm(davis)}
#' \dontrun{bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)}
#' \dontrun{backbone:::class.convert(bb, "matrix")}
#' \dontrun{backbone:::class.convert(bb, "sparseMatrix")}
#' \dontrun{backbone:::class.convert(bb, "igraph")}
#' \dontrun{backbone:::class.convert(bb, "network")}
class.convert <- function(graph, convert = "matrix", extract = FALSE){
  class <- class(graph)

  #### Converting to "matrix" ####
  if (convert == "matrix"){
    if ((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix")) | (methods::is(graph, "Matrix"))) {
      if (dim(graph)[2] == 2){ # for edgelists, assumes bipartite!
        g <- igraph::graph.data.frame(graph, directed = F)
        igraph::V(g)$type <- igraph::V(g)$name %in% graph[,2] #second column of edges is TRUE type
        G <- igraph::get.incidence(g)
        class <- "edgelist"} else {
          if (is.numeric(graph[1]) != TRUE){
            class(graph)<-"numeric"
          }
          G <- graph}
    }
    if (methods::is(graph, "igraph")) {
      if (igraph::is.bipartite(graph)){
        G <- igraph::get.incidence(graph)
      } else {
        if (length(igraph::edge.attributes(graph)) > 0){
          G <- igraph::get.adjacency(graph, attr = "weight")
          } else{G <- igraph::get.adjacency(graph)}
      }
    }
    if (methods::is(graph, "network")) {
      if ("weight" %in% network::list.edge.attributes(graph)){
        G <- as.matrix(graph, attr = "weight")
      } else {G <- as.matrix(graph)}
    }
  ### Remove rows/cols with zero sums ###
    if (extract == FALSE){
      R <- Matrix::rowSums(G)
      C <- Matrix::colSums(G)
      r <- which(R == 0)
      c <- which(C == 0)
      if (length(r)>0){
        G <- G[-r,]
        warning("Rows with a zero sum have been removed from the data. The rows removed were: ",
        paste0(r, " "))
      }
      if (length(c)>0){
        G <- G[,-c]
        warning("Columns with a zero sum have been removed from the data. The columns removed were: ", paste0(c, " "))
      }
    }
  }#end convert to matrix

  #### Converting to "sparseMatrix ####
  if (convert == "sparseMatrix"){
    G <- methods::as(graph, "sparseMatrix")
  }

  #### Converting to "igraph" ####
  if (convert == "igraph"){
    if (-1 %in% graph){
      G <- igraph::graph.adjacency(graph,mode = "undirected", weighted = TRUE)
      } else {G <- igraph::graph.adjacency(graph, mode = "undirected")}
  }

  #### Converting to "network" ####
  if (convert == "network"){
    G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = FALSE)
  }

  #### Converting to "edgelist" ####
  if (convert == "edgelist"){
    g <- igraph::graph.adjacency(graph, weighted = TRUE)
    G <- igraph::get.data.frame(g)
  }
  return(list(class, G))
}

#' Polytope method for finding a matrix that maximizes entropy function
#'
#' @param G matrix, an adjacency matrix representing a graph
#' @details Uses convex optimization via the \link[CVXR]{CVXR-package} to find a matrix \eqn{M}{M} that maximizes the entropy function
#'     where \eqn{M}{M} satisfies the following constraints:
#'     (1) the values of \eqn{M}{M} are between 0 & 1, (2) the row sums of the matrix
#'     are equal to the row sums of the original matrix, (3) the column sums of the matrix
#'     are equal to the column sums of the original matrix.
#' @details This method is utilized in the function \link{sdsm} to compute probabilities of an edge existing in a graph.
#'    Method is called polytope as it is optimizing over the convex hull of the set of matrices (thought of as vectors) with
#'    the same row and column sums as the input.
#' @return matrix containing optimal solution to entropy function under constraints
#' @export
#'
#' @examples
#' polytope(davis)
polytope <- function(G){

  #### Define Variable to solve for ####
  matrix <- CVXR::Variable(dim(G)[1], dim(G)[2])

  #### Define row & column sums ####
  mat1 <- matrix(1, dim(G)[2], 1)
  mat2 <- matrix(1, dim(G)[1], 1)

  #### Define Constraints ####
  constraint1 <- matrix >= 0
  constraint2 <- matrix <= 1
  if (methods::is(G, "sparseMatrix")) {constraint3 <- (matrix%*%mat1) == Matrix::rowSums(G)
  } else {constraint3 <- (matrix%*%mat1) == rowSums(G)}
  if (methods::is(G, "sparseMatrix")) {constraint4 <- t(matrix)%*%mat2 == Matrix::colSums(G)
  } else {constraint4 <- t(matrix)%*%mat2 == colSums(G)}
  constraints <- list(constraint1, constraint2, constraint3, constraint4)

  #### Define Objective, the function to solve ####
  objective <- CVXR::Maximize(sum(CVXR::entr(matrix)+CVXR::entr(1-matrix)))

  #### Define Problem: objective with the constrants ####
  problem <- CVXR::Problem(objective, constraints)

  #### Solve the problem ####
  result <- suppressWarnings(CVXR::psolve(problem))

  ### Warning/Stop if not optimal ###
  if (result$status == "optimal_inaccurate") {warning("polytope result not optimal")}
  if (result$status != "optimal" & result$status != "optimal_inaccurate") {stop("unable to compute SDSM-Polytope")}

  #### Results ####
  new_matrix <- result$getValue(matrix)

  ### Restrict values between 0 and 1 ###
  gr <- which(new_matrix>1)
  new_matrix[gr] <- 1
  le <- which(new_matrix<0)
  new_matrix[le] <- 0

  #### Return Matrix of Probabilities ####
  return(new_matrix)
}


#' curveball algorithm
#'
#' @param M matrix
#'
#' @return rm, a matrix with same row sums and column sums as M, but randomized 0/1 entries.
#' @export
#'
#' @references Algorithm and R implementation: \href{https://www.nature.com/articles/ncomms5114}{Strona, Giovanni, Domenico Nappo, Francesco Boccacci, Simone Fattorini, and Jesus San-Miguel-Ayanz. 2014. “A Fast and Unbiased Procedure to Randomize Ecological Binary Matrices with Fixed Row and Column Totals.” Nature Communications 5 (June). Nature Publishing Group: 4114. DOI:10.1038/ncomms5114.}
#' @examples
#' curveball(davis)
curveball<-function(M){
  #### Define Variables ####
  RC=dim(M)
  R=RC[1]
  C=RC[2]
  hp=list()

  #### Mark Locations of One's ####
  for (row in 1:dim(M)[1]) {hp[[row]]=(which(M[row,]==1))}
  l_hp=length(hp)

  #### Curveball Swaps ####
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
  }

  #### Define and Return Random Matrix ####
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

#' Poisson Binomial distribution computed with Refined Normal Approximation
#'
#' @param kk values where the cdf is to be computed
#' @param pp vector of success probabilities for indicators
#' @param wts the weights for each probability
#'
#' @return cdf, cumulative distribution function
#'
#' @references Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. Computational Statistics & Data Analysis, Vol. 59, pp. 41-51.
#' @details These values are approximated using the Refined Normal Approximation (RNA method).
#'     These functions are originally described by \link[poibin]{ppoibin} and used here under GPL-2 License.
#' @keywords internal
#' @examples
#' \dontrun{probs <- polytope(davis)}
#' \dontrun{P <- davis %*% t(davis)}
#' \dontrun{prob.mat <- matrix(probs, nrow = nrow(davis), ncol = ncol(davis))}
#' \dontrun{prob.imat <- sweep(prob.mat, MARGIN = 2, prob.mat[1,], `*`)}
#' \dontrun{mapply(backbone:::rna, kk= as.data.frame(t(P[1,])), pp = as.data.frame(t(prob.imat)))}
rna <-function(kk,pp,wts=NULL)
  #### Check Arguments ####
{
  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }
  if(is.null(wts))
  {
    wts=rep(1,length(pp))
  }

  #### Define Variables ####
  pp=rep(pp,wts)
  muk=sum(pp)
  sigmak=sqrt(sum(pp*(1-pp)))
  gammak=sum(pp*(1-pp)*(1-2*pp))
  ind=gammak/(6*sigmak^3)
  kk1=(kk+.5-muk)/sigmak

  #### Compute Statistic and Return ####
  vkk.r=stats::pnorm(kk1)+gammak/(6*sigmak^3)*(1-kk1^2)*stats::dnorm(kk1)
  vkk.r[vkk.r<0]=0
  vkk.r[vkk.r>1]=1
  res=vkk.r
  return(res)
}

