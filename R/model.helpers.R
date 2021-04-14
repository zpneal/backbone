#' Converts an input graph object to an adjacency matrix and identifies its characteristics
#'
#' @param graph A graph object of class "matrix", "sparseMatrix", \link{igraph}, matrix or dataframe edgelist, or \link[network]{network}
#'
#' @return a list(summary, G)
#'    `summary` is a dataframe containing characteristics of the supplied object
#'    `G` is an adjacency matrix
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(rbinom(5*5,1,.5),5,5)
#' test <- backbone:::tomatrix(M)
tomatrix <- function(graph){
  class <- class(graph)[1]
  isbipartite <- FALSE

  if (!(methods::is(graph, "matrix")) & !(methods::is(graph, "sparseMatrix")) & !(methods::is(graph, "Matrix")) & !(methods::is(graph, "igraph")) & !(methods::is(graph, "network")) & !(methods::is(graph, "data.frame"))) {stop("input bipartite data must be a matrix, igraph, or network object.")}

  #### Convert matrix-like object ####
  if (((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix")) | (methods::is(graph, "Matrix")))) {
    if (dim(graph)[2] > 3) {
      G <- as.matrix(graph)  #Coerce to matrix
      class(G) <- "numeric"  #Coerce to numeric
      if (any(is.na(G))) {stop("The object contains non-numeric entries")}

      if (dim(G)[1]!=dim(G)[2]) {isbipartite <- TRUE}  #A rectangular matrix are treated as bipartite
      if (dim(G)[1]==dim(G)[2] & !is.null(rownames(G)) & !is.null(colnames(G))) { #A labeled square matrix is treated as bipartite IFF
        if (!identical(rownames(G),colnames(G)) &                                 #the row and column labels differ, and
            !isSymmetric(G)) {                                                    #it is not symmetric
          isbipartite <- TRUE
        }
      }
    }
  }

  #### Convert edge list ####
  if (((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix")) | (methods::is(graph, "Matrix")) | methods::is(graph, "data.frame"))) {
    if (dim(graph)[2] == 2 | dim(graph)[2] == 3) {
      class <- "edgelist"  #Update starting class as edgelist
      if ((methods::is(graph, "data.frame")) == FALSE) {G <- as.data.frame(as.matrix(graph))} else {G <- graph} #Coerce to dataframe if necessary
      colnames(G) <- LETTERS[1:dim(graph)[2]]  #Name columns A, B, and (if present) C
      isbipartite <- length(intersect(G[,1],G[,2])) == 0  #Treat as bipartite if there is no overlap in node lists A and B

      if (isbipartite == TRUE) { #Bipartite
        G <- igraph::graph_from_data_frame(G, directed = F)
        igraph::V(G)$type <- igraph::V(G)$name %in% graph[,2] #second column of edges is TRUE type
        if (dim(graph)[2] == 2) {G <- igraph::as_incidence_matrix(G, sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_incidence_matrix(G, attr = "C", sparse = FALSE)} #Weighted
      }

      if (isbipartite == FALSE) { #Unipartite
        G <- igraph::graph_from_data_frame(G, directed = T)
        if (dim(graph)[2] == 2) {G <- igraph::as_adjacency_matrix(G, type = "both", sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_adjacency_matrix(G, type = "both", attr = "C", sparse = FALSE)} #Weighted
      }
    }
  }

  #### Convert igraph ####
  if (methods::is(graph, "igraph")) {

    if (igraph::is.bipartite(graph) == TRUE){ #Bipartite
      isbipartite <- TRUE
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_incidence_matrix(graph, sparse = FALSE) #Unweighted
      } else {G <- igraph::as_incidence_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted
    }

    if (igraph::is.bipartite(graph) == FALSE){ #Unipartite
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_adjacency_matrix(graph, sparse = FALSE) #Unweighted
      } else {G <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted
    }
  }

  #### Convert statnet ####
  if (methods::is(graph, "network")) {

    if (network::is.bipartite(graph) == TRUE) { #Bipartite
      isbipartite <- TRUE
      if ("weight" %in% network::list.edge.attributes(graph)) {
        G <- network::as.matrix.network(graph, type = "incidence", attrname = "weight")} else { #Weighted
          G <- network::as.matrix.network(graph, type = "incidence")} #Unweighted
    }

    if (network::is.bipartite(graph) == FALSE) { #Unipartite
      if ("weight" %in% network::list.edge.attributes(graph)) {
        G <- network::as.matrix.network(graph, type = "adjacency", attrname = "weight")} else { #Weighted
          G <- network::as.matrix.network(graph, type = "adjacency")} #Unweighted
    }
  }

  #### If the graph is bipartite, remove rows/columns with zero sums ####
  if (isbipartite){
    R <- Matrix::rowSums(G)
    C <- Matrix::colSums(G)
    r <- which(R == 0)
    c <- which(C == 0)
    if (length(r)>0){G <- G[-r,]}
    if (length(c)>0){G <- G[,-c]}
  }

  #### Summary dataframe ####
  if (isbipartite) {issymmetric <- FALSE} else {issymmetric <- isSymmetric(G)}
  isweighted <- any(!(G%in%c(0,1)))
  summary <- data.frame(
    class = class,
    bipartite = isbipartite,
    symmetric = issymmetric,
    weighted = isweighted)

  #### Report input type and modifications ####
  if (issymmetric) {dir <- "undirected"} else {dir <- "directed"}
  if (isweighted) {weigh <- "a weighted"} else {weigh <- "an unweighted"}
  if (isbipartite) {
    message(paste0("This ", class, " object looks like ", weigh, " bipartite network of ", nrow(G), " agents and ", ncol(G), " artifacts."))
    if (length(r)>0) {message("These zero-sum rows have been removed from the data: ", paste0(r, " "))}
    if (length(c)>0) {message("These zero-sum columns have been removed from the data: ", paste0(c, " "))}
  }
  if (!isbipartite) {message(paste0("This ", class, " object looks like ", weigh, " ", dir, " network containing ", nrow(G), " nodes."))}

  return(list(summary = summary, G = G))
}

#' Converts a backbone adjacency matrix to an object of specified class
#'
#' @param graph a matrix
#' @param convert class to convert to, one of "matrix", "sparseMatrix", "igraph", "edgelist", or "network"
#'
#' @return backbone graph: Binary or signed backbone graph of class given in parameter `convert`.
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
#' test <- backbone:::frommatrix(M, "Matrix")
frommatrix <- function(graph, convert = "matrix"){

  if (convert == "matrix"){G <- graph}

  if (convert == "Matrix"){G <- Matrix::Matrix(graph)}

  if (convert == "sparseMatrix"){G <- methods::as(graph, "sparseMatrix")}

  if (convert == "network"){
    if (isSymmetric(graph)) {G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = FALSE)}
    if (!isSymmetric(graph)) {G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = TRUE)}
  }

  if (convert == "edgelist"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
    G <- igraph::as_data_frame(G, what = "edges")
  }

  if (convert == "igraph"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
    if (gsize(G)!=0) {igraph::E(G)$sign <- igraph::E(G)$weight}  #To facilitate use with library(signnet)
  }

  return(G)
}



