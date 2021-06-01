#' Compute fixed column sums / Poisson binomial backbone probabilities
#'
#' `fixedcol` computes the probability of observing
#'     a higher or lower edge weight using the Poisson binomial distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param method string: Specifies the method of the Poisson Binomial distribution computation used by the "ppbinom" function in \link[PoissonBinomial]{PoissonBinomial-Distribution}.
#'     "RefinedNormal" gives quick, very accurate approximations, while "DivideFFT" gives the quickest exact computations.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the column vertex degrees are fixed but the row vertex degrees are allowed to vary.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: class, model name, number of rows, number of columns, and running time.
#' @export
#' @examples
#' fixedcol(davis)
fixedcol <- function(B,
                     method = "RefinedNormal"){

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
  cs <- colSums(B)
  m = dim(B)[1]

  #### Poisson Binomial Distribution ####
  ### Compute parameters for the diagonal ###
  pd <- cs/m
  diagonaln <- PoissonBinomial::ppbinom(x = diag(P), probs = pd, method = method)
  diagonalp <- PoissonBinomial::ppbinom(x = diag(P)-1, probs = pd, method = method, lower.tail = FALSE)

  ### Compute parameters for the off diagonal###
  pt <- ((cs*(cs-1))/(m*(m-1)))

  ### Poisson binomial distribution probabilities ###
  Negative <- matrix(PoissonBinomial::ppbinom(x = P, probs = pt, method = method),nrow = m, ncol = m,byrow = TRUE)
  Positive <- matrix(PoissonBinomial::ppbinom(x = (P-1), probs = pt, method = method, lower.tail = FALSE),nrow = m, ncol = m,byrow = TRUE)
  diag(Negative) <- diagonaln
  diag(Positive) <- diagonalp

  ### Add back in rownames ###
  rownames(Positive) <- rownames(B)
  colnames(Positive) <- rownames(B)
  rownames(Negative) <- rownames(B)
  colnames(Negative) <- rownames(B)

  ### Insert NAs for p-values along diagonal
  diag(Positive) <- NA
  diag(Negative) <- NA

  #### Compile Summary ####
  r <- rowSums(B)
  c <- colSums(B)

  a <- c("Model", "Input Class", "Bipartite", "Symmetric", "Weighted", "Number of Rows", "Number of Columns")
  b <- c("Fixed Column Model", convert$summary$class, convert$summary$bipartite, convert$summary$symmetric, convert$summary$weighted, dim(B)[1], dim(B)[2])

  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}
