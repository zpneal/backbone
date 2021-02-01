#' The stochastic degree sequence model (sdsm) for backbone probabilities
#'
#' `sdsm` computes the probability of edge weights being
#'     above or below the observed edge weights in a bipartite projection
#'     using the stochastic degree sequence model.
#'     Once computed, use \code{\link{backbone.extract}} to return
#'     the backbone matrix for a given alpha value.
#'
#' @param B graph: Bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#' @param progress Boolean: If \link[utils]{txtProgressBar} should be used to measure progress
#' @param ... optional arguments
#' @details Specifically, the sdsm function compares an edge's observed weight in the projection \code{B*t(B)}
#'    to the distribution of weights expected in a projection obtained from a random bipartite network where
#'    both the row vertex degrees and column vertex degrees are approximately fixed.
#' @details sdsm uses the Bipartite Configuration Model \link{bicm} (Saracco et. al (2015, 2017)) to compute probabilities for the Poisson binomial distribution.
#'
#' @details The "backbone" S3 class object returned is composed of two matrices, and a summary dataframe.
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value in the projection, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @references sdsm: \href{https://www.sciencedirect.com/science/article/abs/pii/S0378873314000343}{Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. Social Networks, 39, Elsevier: 84-97. DOI: 10.1016/j.socnet.2014.06.001}
#' @references bicm: \href{https://doi.org/10.1088/1367-2630/aa6b38}{Saracco, F., Straka, M. J., Clemente, R. D., Gabrielli, A., Caldarelli, G., & Squartini, T. (2017). Inferring monopartite projections of bipartite networks: An entropy-based approach. New Journal of Physics, 19(5), 053022. DOI: 10.1088/1367-2630/aa6b38}
#' @references bicm: \href{https://doi.org/10.1038/srep10595}{Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. Scientific Reports, 5(1), 10595. DOI: 10.1038/srep10595}
#' @export
#'
#' @examples
#'sdsm_probs <- sdsm(davis)
sdsm <- function(B,
                 progress = FALSE,
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

  if (!(methods::is(B, "matrix")) & !(methods::is(B, "sparseMatrix")) & !(methods::is(B, "igraph")) & !(methods::is(B, "network"))) {stop("input bipartite data must be a matrix, igraph, or network object.")}

  ### Run Time ###
  run.time.start <- Sys.time()
  message(paste0("Finding the distribution using SDSM with BiCM probabilities."))

  #### Class Conversion ####
  convert <- class.convert(B, "matrix")
  class <- convert[[1]]
  B <- convert[[2]]

  #### Bipartite Projection ####
  ### If sparse matrix input, use sparse matrix operations ###
  if (methods::is(B, "sparseMatrix")) {
    B <- Matrix::Matrix(B, sparse = T)
  }
  P <- Matrix::tcrossprod(B)

  ### Create Positive and Negative Matrices to hold backbone ###
  Positive <- matrix(0, nrow(P), ncol(P))
  Negative <- matrix(0, nrow(P), ncol(P))

  #### Compute Probabilities for SDSM ####
  prob.mat <- bicm(graph=B,progress=progress,...)

  #### Assemble and Probabilities ####
  rows <- dim(prob.mat)[1]

  #### Compute Null Edge Weight Distributions Using Poisson Binomial RNA ####

  for (i in 1:rows){
    ### Compute prob.mat[i,]*prob.mat[j,] for each j ###
    prob.imat <- sweep(prob.mat, MARGIN = 2, prob.mat[i,], `*`)

    ### Find cdf, below or equal to value for negative, above or equal to value for positive ###
    ### Using RNA approximation ###
    negative <- as.array(mapply(rna, kk= as.data.frame(t(P[i,])), pp = as.data.frame(t(prob.imat))))
    positive <- as.array((1- mapply(rna, kk=(as.data.frame(t(P[i,])-1)), pp = as.data.frame(t(prob.imat)))))

    ### Set values in Positive & Negative matrices ###
    Positive[i,] <- positive
    Negative[i,] <- negative
  } #end for i in rows
  rownames(Positive) <- rownames(B)
  colnames(Positive) <- rownames(B)
  rownames(Negative) <- rownames(B)
  colnames(Negative) <- rownames(B)

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- Matrix::rowSums(B)
  c <- Matrix::colSums(B)

  a <- c("Input Class", "Model", "Method", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Stochastic Degree Sequence Model", "BiCM", dim(B)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(B)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)

} #end sdsm function
