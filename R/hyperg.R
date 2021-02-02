#' Compute hypergeometric backbone probabilities
#'
#' `hyperg` computes the probability of observing
#'     a higher or lower edge weight using the hypergeometric distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: Bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the row vertex degrees are fixed but the column vertex degrees are allowed to vary.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @references {Tumminello, Michele and Miccichè, Salvatore and Lillo, Fabrizio and Piilo, Jyrki and Mantegna, Rosario N. 2011. "Statistically Validated Networks in Bipartite Complex Systems." PLOS ONE, 6(3), \doi{10.1371/journal.pone.0017994}}
#' @references {Neal, Zachary. 2013. “Identifying Statistically Significant Edges in One-Mode Projections.” Social Network Analysis and Mining 3 (4). Springer: 915–24. \doi{10.1007/s13278-013-0107-y}}
#' @export
#'
#' @examples
#' hyperg_probs <- hyperg(davis)

hyperg <- function(B){

  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- class.convert(B)
  class <- convert[[1]]
  B <- convert[[2]]

  #### Argument Checks ####
  if (!(methods::is(B, "matrix")) & !(methods::is(B, "sparseMatrix"))) {stop("input bipartite data must be a matrix")}
  message("Finding the distribution using hypergeometric distribution")

  #### Bipartite Projection ####
  if (methods::is(B, "sparseMatrix")) {
    P <- Matrix::tcrossprod(B)
    rs <- Matrix::rowSums(B)
  } else {
    P <- tcrossprod(B)
    rs <- rowSums(B)
  }

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
  if (methods::is(B, "sparseMatrix")) {
    r <- Matrix::rowSums(B)
    c <- Matrix::colSums(B)
  } else {
    r <- rowSums(B)
    c <- colSums(B)
  }
  a <- c("Input Class", "Model", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Hypergeometric Model", dim(B)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(B)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}

