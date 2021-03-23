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

  if (convert == "network"){G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = TRUE)}

  if (convert == "edgelist"){
    G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)
    G <- igraph::as_data_frame(G, what = "edges")
  }

  if (convert == "igraph"){
    G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)
    igraph::E(G)$sign <- igraph::E(G)$weight  #To facilitate use with library(signnet)
  }

  return(G)
}
