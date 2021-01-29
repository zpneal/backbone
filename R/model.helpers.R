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
    if (igraph::is.weighted(G) == T){
      igraph::E(G)$sign <- igraph::E(G)$weight
    }
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
rna <-function(kk,pp,wts=NULL){
  #### Check Arguments ####
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


#' bipartite.null: generates a backbone object from a bipartite matrix using a null model defined by constraining row and/or column sums.
#' @param B graph: Bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#' @param rows boolean: TRUE if the row sums should be constrained by the null model, FALSE if not.
#' @param cols boolean: TRUE if the column sums should be constrained by the null model, FALSE if not.
#' @param trials integer: number of monte carlo trials used to estimate the \link{fdsm} null model (rows = TRUE, cols = TRUE)
#' @param ... optional arguments
#' @details When only rows are constrained, the hypergeometric null model (\link{hyperg}) is used.
#'     When rows and columns are constrained, the stochastic degree sequence model (\link{sdsm}) is used.
#'     When rows and columns are constrained and trials are specified, the fixed degree sequence model (\link{fdsm}) is used.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @export
#'
#' @examples bipartite.null(davis, rows = TRUE, cols = FALSE) #runs hyperg on davis data
bipartite.null <- function(B,
                      rows = TRUE,
                      cols = TRUE,
                      trials = NULL,
                      ...){
 if ((rows==TRUE)&(cols==TRUE)){
   if (is.null(trials)){
     return(sdsm(B,progress = TRUE,...))
   } else {
     return(fdsm(B,trials = trials, progress = TRUE))
   } #end else
 } #end if r/c T
 else if ((rows == TRUE)&(cols == FALSE)){
   return(hyperg(B))
 }
 else if ((rows == FALSE) & (cols == TRUE)){
   stop("This null model is not currently implemented.")
 } else if ((rows == FALSE) & (cols == FALSE)){
   stop("This null model does not exist.")
 }
} #end bipartite.null


