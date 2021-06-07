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
#' @details The fixedfill function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the number of edges present is equal to the number of edges in B. When B is large, this function may be impractically
#'     slow and may return a backbone object that contains NaN values.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: class, model name, number of rows, number of columns, and running time.
#' @references {Neal, Domagalski, and Sagan. 2021. "Comparing Models for Extracting the Backbone of Bipartite Projections."} \href{https://arxiv.org/abs/2105.13396}{arXiv:2105.13396 cs.SI})
#' @export
#'
#' @examples
#' fixed_probs <- fixedfill(davis)

fixedfill <- function(B){

  #### Class Conversion ####
  convert <- tomatrix(B)
  class <- convert$summary$class
  B <- convert$G
  if (convert$summary$weighted==TRUE){stop("Graph must be unweighted.")}
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
  }

  #### Bipartite Projection ####
  P <- tcrossprod(B)
  rs <- rowSums(B)

  #### Compute Probabilities ####
  m <- dim(B)[1]
  n <- dim(B)[2]
  f <- sum(B)

  ### Off Diagonal Values ###
  ## This computes log of k! ##
  logsum <- function(k){
    if (k==0){
      return(0)
    }
    return(sum(log(1:k)))
  }

  ## This computes log of (n choose k) ##
  logbinom <- function(n,k){
    if (k == 0){
      return(0)
    }
    else if (k == 1){
      return(log(n))
    }
    else {
      x <- sum(log(n:(n-k+1)))
      y <- sum(log(k:1))
      return(x-y)
    }
  }

  prob_log <- function(k) {
    lb <- max(0, n + k - f)
    ub <- min(n - k, (m - 1) * n + k - f)
    range <- lb:ub
    logvalues <- matrix(0, nrow = 1, ncol = length(range))
    i = 1
    for (r in range){
      logvalues[i] <- (log(2^(n-k-r))+logsum(n)-logsum(k)-logsum(r)-logsum(n-k-r)+logbinom((m-2)*n,f-n-k+r)-logbinom(m*n,f))
      i <- i+1
    }
    return(sum((exp(logvalues))))
  }

  maxk <- max(P)  #Largest observed k
  probs <- sapply(0:maxk, FUN = prob_log)  #Probability of observing each k, for 0 <= k <= maxk
  probs <- c(probs, 1 - sum(probs))  #Add one more entry for probability of observing any k > maxk (upper tail of PMF)

  #### Create Positive and Negative Probability Matrices ####
  Positive <- apply(P, c(1,2), FUN = function(k)sum(probs[(k+1):(maxk+2)]))  #Sum of probabilities Pij <= k <= maxk and beyond
  Negative <- apply(P, c(1,2), FUN = function(k)sum(probs[1:(k+1)]))  #Sum of probabilities 0 <= k <= Pij

  ### Insert NAs for p-values along diagonal
  #diagonal <- diag(P)
  #diagn <- stats::phyper(diagonal, n, (m-1)*n, f-diagonal, lower.tail = TRUE)
  #diagp <- stats::phyper(diagonal-1, n, (m-1)*n, f-diagonal, lower.tail=FALSE)
  diag(Positive) <- NA
  diag(Negative) <- NA

  #### Compile Summary ####
  r <- rowSums(B)
  c <- colSums(B)

  a <- c("Model", "Input Class", "Bipartite", "Symmetric", "Weighted", "Number of Rows", "Number of Columns")
  b <- c("Fixed Fill Model", convert$summary$class, convert$summary$bipartite, convert$summary$symmetric, convert$summary$weighted, dim(B)[1], dim(B)[2])

  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}
