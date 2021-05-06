#' Compute fixed row sums / hypergeometric backbone probabilities
#'
#' `fixedrow` computes the probability of observing
#'     a higher or lower edge weight using the hypergeometric distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the row vertex degrees are fixed but the column vertex degrees are allowed to vary.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: class, model name, number of rows, number of columns, and running time.
#' @references {Tumminello, Michele and Miccichè, Salvatore and Lillo, Fabrizio and Piilo, Jyrki and Mantegna, Rosario N. 2011. "Statistically Validated Networks in Bipartite Complex Systems." PLOS ONE, 6(3), \doi{10.1371/journal.pone.0017994}}
#' @references {Neal, Zachary. 2013. “Identifying Statistically Significant Edges in One-Mode Projections.” Social Network Analysis and Mining 3 (4). Springer: 915–24. \doi{10.1007/s13278-013-0107-y}}
#' @export
#'
#' @examples
#' fixedrow_probs <- fixedrow(davis)

fixedrow <- function(B){

  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(B)
  class <- convert$summary[[1]]
  B <- convert$G
  if (convert$summary[[4]]==TRUE){stop("Graph must be unweighted.")}
  if (convert$summary[[2]]==FALSE){warning("This object is being treated as a bipartite network.")}


  #### Bipartite Projection ####
  P <- tcrossprod(B)
  rs <- rowSums(B)

  #### Hypergeometric Distribution ####

  ### Set up df for values ###
  df <- data.frame(as.vector(P))
  names(df)[names(df)=="as.vector.P."] <- "projvalue"

  ### Compute row sums ###
  df$row_sum_i <- rep(rs, times = nrow(B))

  ### Match each row sum i with each row sum j and their Pij value ###
  df$row_sum_j <- rep(rs, each = nrow(B))

  ### Compute different in number of artifacts and row sum ###
  df$diff <- ncol(B)-df$row_sum_i

  ### Probability of Pij or less ###
  df$hgl <- stats::phyper(df$projvalue, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)

  ### Probability of Pij or more ###
  df$hgu <- stats::phyper(df$projvalue-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)

  #### Create Positive and Negative Probability Matrices ####
  Positive <- matrix(as.numeric(df$hgu), nrow = nrow(B), ncol = nrow(B))
  Negative <- matrix(as.numeric(df$hgl), nrow = nrow(B), ncol = nrow(B))

  ### Add back in rownames ###
  rownames(Positive) <- rownames(B)
  colnames(Positive) <- rownames(B)
  rownames(Negative) <- rownames(B)
  colnames(Negative) <- rownames(B)

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- rowSums(B)
  c <- colSums(B)

  a <- c("Model", "Input Class", "Bipartite", "Symmetric", "Weighted", "Number of Rows", "Number of Columns", "Running Time (secs)")
  b <- c("Hypergeometric Model", class[1], convert$summary$bipartite, convert$summary$symmetric, convert$summary$weighted, dim(B)[1], dim(B)[2], as.numeric(total.time))

  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}


#' A wrapper for the \link{fixedrow} function
#'
#' `hyperg` computes the probability of observing
#'     a higher or lower edge weight using the hypergeometric distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the row vertex degrees are fixed but the column vertex degrees are allowed to vary.
#' @return \link{fixedrow}
#' @references {Tumminello, Michele and Miccichè, Salvatore and Lillo, Fabrizio and Piilo, Jyrki and Mantegna, Rosario N. 2011. "Statistically Validated Networks in Bipartite Complex Systems." PLOS ONE, 6(3), \doi{10.1371/journal.pone.0017994}}
#' @references {Neal, Zachary. 2013. “Identifying Statistically Significant Edges in One-Mode Projections.” Social Network Analysis and Mining 3 (4). Springer: 915–24. \doi{10.1007/s13278-013-0107-y}}
#' @export
#'
#' @examples
#' hyperg_probs <- hyperg(davis)

hyperg <- function(B){
  return(fixedrow(B))
}


