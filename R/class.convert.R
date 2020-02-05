#' Convert graph object to adjacency matrix
#'
#' @param graph, matrix, sparse matrix, igraph, edgelist, or network object
#'
#' @return list(class, adjacency), a list containing the class of parameter graph, and the adjacency matrix of the graph
#'
#' @examples
graph_to_adjacency <- function(graph){
  class <- class(graph)
  if ((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix"))) {adjacency <- graph}
  #if (((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix"))) & (dim(graph)[2]==2)) {adjacency <- igraph::get.incidence(igraph::graph_from_data_frame(graph))}
  if (methods::is(graph, "igraph")) {
    if (igraph::is.bipartite(graph)){
      adjacency <- igraph::get.incidence(graph)
    }
    else { adjacency <- igraph::get.adjacency(graph)}
  }
  if (methods::is(graph, "network")) {adjacency <- as.matrix(graph, attr = "weight")}
  return(list(class, adjacency))
}

### Notes about edgelist & graph objects
# edgelist is not it's own class, could be identified by 2 columns
# if we put in a dim argument, the igraph objects are throwing errors, but shouldn't be? something wrong in the logic flow
# d3 is sparse but b3 isn't?


#' Converts adjacency matrix into a graph object
#'
#' @param adj, an adjacency matrix of a graph
#' @param class, character, class of a graph object
#'
#' @return graph object of class `class`
#'
#' @examples
adjacency_to_graph <- function(adj, class){
  if (class == "matrix"){graph <- adj}
  if (class == "dgCMatrix"){graph <- as(adj, "sparseMatrix")}
  if (class == "igraph"){graph <- igraph::graph.adjacency(adj, mode = "undirected", weighted = TRUE)}
  if (class == "network"){graph <- network::network(adj, ignore.eval = FALSE, names.eval = "weight", directed = FALSE, loops = TRUE)}
  return(graph)
}



