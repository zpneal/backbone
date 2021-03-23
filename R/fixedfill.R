#' Compute fixed fill backbone probabilities
#'
#' `fixedfill` computes the probability of observing
#'     a higher or lower edge weight.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the number of edges present is equal to the number of edges in B.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @export
#'
#' @examples
#' fixed_probs <- fixedfill(davis)

fixedfill <- function(B){

  #### Argument Checks ####
  if (!(methods::is(B, "matrix")) & !(methods::is(B, "sparseMatrix")) & !(methods::is(B, "igraph")) & !(methods::is(B, "network"))) {stop("input bipartite data must be a matrix, igraph, or network object.")}

  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(B)
  class <- convert$summary[[1]]
  B <- convert$G

  if (any(!B%in%c(0,1))){stop("Graph must be unweighted.")}

  #### Bipartite Projection ####
  P <- tcrossprod(B)
  rs <- rowSums(B)

  #### Compute Probabilities ####

  m <- dim(B)[1]
  n <- dim(B)[2]
  f <- sum(B)

  ### Diagonal Values ###
  diagonal <- diag(P)
  diagn <- stats::phyper(diagonal, n, (m-1)*n, f-diagonal, lower.tail = TRUE)
  diagp <- stats::phyper(diagonal-1, n, (m-1)*n, f-diagonal, lower.tail=FALSE)

  ### Off Diagonal Values ###
  prob <- function(k){
    lb <- max(0,n+k-f)
    ub <- min(n-k, (m-1)*n+k-f)
    range <- lb:ub
    return((choose(n,k))/(choose(m*n, f))*sum(2^(n-k-range)*choose(n-k,range)*choose(((m-2)*n), f-k-n+range)))
  }
  max <- max(P)
  probs <- sapply(0:max, FUN = prob)

  #### Create Positive and Negative Probability Matrices ####
  Negative <- apply(P, c(1,2), FUN = function(k)sum(probs[1:(k+1)]))
  Positive <- apply(P, c(1,2), FUN = function(k) 1- sum(probs[1:k]))
  diag(Negative) <- diagn
  diag(Positive) <- diagp

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- rowSums(B)
  c <- colSums(B)
  a <- c("Input Class", "Model", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Fixed Fill Model", dim(B)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(B)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}
