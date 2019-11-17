#' Compute hypergeometric backbone
#'
#' `hyperg` computes the probability of observing
#'     a higher or lower edge weight using the hypergeometric distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B Matrix: Bipartite network
#'
#' @return list(positive, negative, summary).
#' positive gives matrix of probability of ties above the observed value.
#' negative gives matrix of probability of ties below the observed value.
#' summary: a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#'
#' @references \href{https://doi.org/10.1007/s13278-013-0107-y}{Neal, Zachary. 2013. “Identifying Statistically Significant Edges in One-Mode Projections.” Social Network Analysis and Mining 3 (4). Springer: 915–24. DOI:10.1007/s13278-013-0107-y.}
#' @export
#'
#' @examples
#' hypergeometric_bb <- hyperg(davis)

hyperg <- function(B){

  #Run Time
  run.time.start <- Sys.time()

  #Argument Checks
  if (class(B)!="matrix") {stop("input bipartite data must be a matrix")}
  message("Finding the Backbone using Hypergeometric Distribution")

  P <- tcrossprod(B)
  df <- data.frame(as.vector(P))
  names(df)[names(df)=="as.vector.P."] <- "projvalue"

  #Compute row sums
  df$row_sum_i <- rep(rowSums(B), times = nrow(B))

  #Match each row sum i with each row sum j and their Pij value
  df$row_sum_j <- rep(rowSums(B), each = nrow(B))

  #Compute different in number of artifacts and row sum
  df$diff <- ncol(B)-df$row_sum_i

  #Probability of Pij or less
  df$hgl <- stats::phyper(df$projvalue, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)

  #Probability of Pij or more
  df$hgu <- stats::phyper(df$projvalue-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)

  Positive <- matrix(as.numeric(df$hgu), nrow = nrow(B), ncol = nrow(B))
  Negative <- matrix(as.numeric(df$hgl), nrow = nrow(B), ncol = nrow(B))

  #Add back in rownames
  rownames(Positive) <- rownames(B)
  colnames(Positive) <- rownames(B)
  rownames(Negative) <- rownames(B)
  colnames(Negative) <- rownames(B)

  #Run Time
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start), 2))

  #Compile Summary
  a <- c("Model", "Number of Rows", "Skew of Row Sums", "Number of Columns", "Skew of Column Sums", "Running Time")
  b <- c("Hypergeometric Model", dim(B)[1], round((sum((r-mean(r))**3))/((length(r))*((sd(r))**3)), 5), dim(B)[2], round((sum((c-mean(c))**3))/((length(c))*((sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  return(list(positive = Positive, negative = Negative, summary = model.summary))
}

