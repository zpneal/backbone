#' The stochastic degree sequence model (sdsm) for backbone probabilities
#'
#' `sdsm` computes the probability of edge weights being
#'     above or below the observed edge weights in a bipartite projection
#'     using the stochastic degree sequence model.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param method string: Specifies the method of the Poisson Binomial distribution computation used by the "ppbinom" function in \link[PoissonBinomial]{PoissonBinomial-Distribution}.
#'     "RefinedNormal" gives quick, very accurate approximations, while "DivideFFT" gives the quickest exact computations.
#' @param ... optional arguments
#' @details The sdsm function compares an edge's observed weight in the projection \code{B*t(B)}
#'    to the distribution of weights expected in a projection obtained from a random bipartite network where
#'    both the row vertex degrees and column vertex degrees are approximately fixed. It uses the Bipartite Configuration Model \link{bicm} (Saracco et. al (2015, 2017)) to compute probabilities for the Poisson binomial distribution.
#'
#' @details The "backbone" S3 class object returned is composed of two matrices, and a summary dataframe.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: class, model name, number of rows, number of columns, and running time.
#' @references sdsm: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. Social Networks, 39, Elsevier: 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @references bicm: {Saracco, F., Straka, M. J., Clemente, R. D., Gabrielli, A., Caldarelli, G., & Squartini, T. (2017). Inferring monopartite projections of bipartite networks: An entropy-based approach. New Journal of Physics, 19(5), 053022. \doi{10.1088/1367-2630/aa6b38}}
#' @references bicm: {Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. Scientific Reports, 5(1), 10595. \doi{10.1038/srep10595}}
#' @export
#'
#' @examples
#'sdsm_probs <- sdsm(davis)
sdsm <- function(B,
                 method = "RefinedNormal",
                 ...){

  #### Argument Checks ####
  args <- match.call()
  exist <- ("model" %in% names(args))
  if (exist == TRUE){
    message("This model is deprecated. SDSM now uses the 'bicm' model.
             To run an older model, you must install a previous version of backbone.
             This can be done by using:
            ''require(devtools)'' and
            ''install_version(''backbone'', version = '1.2.2')")
  }

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

  ### Create Positive and Negative Matrices to hold backbone ###
  Positive <- matrix(0, nrow(P), ncol(P))
  Negative <- matrix(0, nrow(P), ncol(P))

  #### Compute Probabilities for SDSM ####
  prob.mat <- bicm(graph=B,...)

  #### Assemble and Probabilities ####
  rows <- dim(prob.mat)[1]

  #### Compute Null Edge Weight Distributions Using Poisson Binomial Distribution ####
  for (i in 1:rows){
    ### Compute prob.mat[i,]*prob.mat[j,] for each j ###
    prob.imat <- sweep(prob.mat, MARGIN = 2, prob.mat[i,], `*`)

    ### Find cdf, below or equal to value for negative, above or equal to value for positive ###
    negative <- as.array(mapply(PoissonBinomial::ppbinom, x= as.data.frame(t(P[i,])), probs = as.data.frame(t(prob.imat)), method = "RefinedNormal"))
    positive <- as.array(mapply(PoissonBinomial::ppbinom, x=(as.data.frame(t(P[i,])-1)), probs = as.data.frame(t(prob.imat)), method = "RefinedNormal", lower.tail = FALSE))

    ### Set values in Positive & Negative matrices ###
    Positive[i,] <- positive
    Negative[i,] <- negative
  } #end for i in rows

  ### Add back in rownames
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
  b <- c("Stochastic Degree Sequence Model", convert$summary$class, convert$summary$bipartite, convert$summary$symmetric, convert$summary$weighted, dim(B)[1], dim(B)[2])

  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)

} #end sdsm function
