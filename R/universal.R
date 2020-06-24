#' Compute universal threshold backbone
#'
#' `universal` returns a backbone graph in which edge weights are set to
#'    1 if above the given upper parameter threshold,
#'    set to -1 if below the given lower parameter threshold, and are 0 otherwise.
#'
#' @param M graph: Bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#' @param upper Real or FUN: upper threshold value or function to be applied to the edge weights. Default is 0.
#' @param lower Real or FUN: lower threshold value or function to be applied to the edge weights. Default is NULL.
#' @param bipartite Boolean: TRUE if bipartite matrix, FALSE if weighted matrix. Default is FALSE.
#'
#' @return backbone, a list(backbone, summary). The `backbone` object is a graph object of the same class as M.
#'     The `summary` contains a data frame summary of the inputted matrix and the model used including:
#'     model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @export
#'
#' @examples
#' test <- universal(davis%*%t(davis), upper = function(x)mean(x)+sd(x), lower=function(x)mean(x))
#' test2 <- universal(davis, upper = function(x)mean(x)+2*sd(x), lower = 2, bipartite = TRUE)
#' test3 <- universal(davis, upper = 4, lower = 2, bipartite = TRUE)

universal <- function(M,
                      upper = 0,
                      lower = NULL,
                      bipartite = FALSE){
  #### Argument Checks ####
  if (!(methods::is(upper, "function")) & (!(methods::is(upper, "numeric")))) {stop("upper must be either function or numeric")}
  if (!(methods::is(lower, "function")) & (!(methods::is(lower, "numeric"))) & (!(methods::is(lower, "NULL")))) {stop("lower must be either function or numeric")}

  ### Run Time ###
  run.time.start <- Sys.time()

  ### Class Conversion ###
  convert <- class.convert(M)
  class <- convert[[1]]
  M <- convert[[2]]

  #### Bipartite Projection ####
  if (bipartite == TRUE){
    if (methods::is(M, "sparseMatrix")) {
      P <- Matrix::tcrossprod(M)
    } else {
      P <- tcrossprod(M)
    }
  } else {
    P <- M
  }

  #### Set Threshold Values ####
  if (class(upper) == "function"){
    ut <- upper(P)
  } else {ut <- upper}

  if (class(lower) == "function"){
    lt <- lower(P)
  } else {lt <- lower}

  #### Create Backbone Matrix ####
  backbone <- matrix(0, nrow(P), ncol(P))
  negative <- (P<lt)+0
  positive <- (P>ut)+0

  if (length(lower) > 0){
    backbone <- backbone - negative
  }

  backbone <- backbone + positive

  diag(backbone) <- 0

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  if (methods::is(M, "sparseMatrix")) {
    r <- Matrix::rowSums(M)
    c <- Matrix::colSums(M)
  } else {
    r <- rowSums(M)
    c <- colSums(M)
  }
  a <- c("Input Class", "Model", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time")
  b <- c(class[1], "Universal Threshold", dim(M)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(M)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Convert to Indicated Class Object ####
  backbone_converted <- class.convert(backbone, class[1])

  #### Return Backbone and Summary ####
  bb <- list(backbone = backbone_converted[[2]], summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}


