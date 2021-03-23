#' Compute fixed column sums / Poisson binomial backbone probabilities
#'
#' `fixedcol` computes the probability of observing
#'     a higher or lower edge weight using the Poisson binomial distribution.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
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
fixedcol <- function(B){

  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(B)
  class <- convert$summary[[1]]
  B <- convert$G
  if (convert$summary[[2]]==FALSE){stop("Graph must be bipartite.")}
  if (convert$summary[[4]]==TRUE){stop("Graph must be unweighted.")}

  #### Bipartite Projection ####
  ### If sparse matrix input, use sparse matrix operations ###
  P <- tcrossprod(B)
  cs <- colSums(B)

  m = dim(B)[1]

  #### Poisson Binomial Distribution ####
  ### Compute parameters for the diagonal ###
  pd <- cs/m
  diagonaln <- rna(kk = diag(P), pp = pd)
  diagonalp <- 1-rna(kk = diag(P)-1, pp = pd)

  ### Compute parameters for the off diagonal###
  pt <- ((cs*(cs-1))/(m*(m-1)))

  ### Use RNA to get probabilities ###
  Negative <- as.array(rna(kk = P, pp = pt))
  Positive <- as.array(1-rna(kk = (P-1), pp = pt))
  diag(Negative) <- diagonaln
  diag(Positive) <- diagonalp

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
