#' Compute fixed column sums / Poisson binomial backbone probabilities
#'
#' `fixedcol` computes the probability of observing
#'     a higher or lower edge weight using the Poisson binomial distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param method string: Specifies the method of the Poisson Binomial distribution computation used by the ``ppbinom" function in \link[PoissonBinomial]{PoissonBinomial-Distribution}.
#'     "RefinedNormal" gives quick, very accurate approximations, while "DivideFFT" gives the quickest exact computations.
#'
#' @details Specifically, this function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the column vertex degrees are fixed but the row vertex degrees are allowed to vary.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @export
#' @examples
#' fixedcol(davis)
fixedcol <- function(B,
                     method = "RefinedNormal"){

  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(B)
  class <- convert$summary[[1]]
  B <- convert$G
  if (convert$summary[[4]]==TRUE){stop("Graph must be unweighted.")}
  if (convert$summary[[2]]==FALSE){warning("This object is being treated as a bipartite network.")}


  #### Bipartite Projection ####
  ### If sparse matrix input, use sparse matrix operations ###
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

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- rowSums(B)
  c <- colSums(B)

  a <- c("Input Class", "Model", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Poisson Binomial Sequence Model", dim(B)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(B)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}
