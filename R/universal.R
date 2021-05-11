#' Compute universal threshold backbone
#'
#' `universal` returns a backbone graph in which edge weights are set to
#'    1 if above the given upper parameter threshold,
#'    set to -1 if below the given lower parameter threshold, and are 0 otherwise.
#'
#' @param M graph: Graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param upper Real, FUN, or NULL: upper threshold value or function to be applied to the edge weights. Default is NULL.
#' @param lower Real, FUN, or NULL: lower threshold value or function to be applied to the edge weights. Default is NULL.
#' @param bipartite Boolean: TRUE if bipartite graph, FALSE if weighted graph. Default is NULL. If TRUE, input graph should be unweighted.
#' @param narrative Boolean: TRUE if suggested text for a manuscript is to be returned
#'
#' @details If both `upper` and `lower` are `NULL`, a weighted projection is returned.
#' @details If `bipartite` is `NULL`, the function tries to guess at whether the data is bipartite or unipartite based on its shape.
#' @return backbone, a list(backbone, summary). The `backbone` object is a graph object of the same class as M.
#'     The `summary` contains a data frame summary of the inputted matrix and the model used including:
#'     model name, number of rows, number of columns, and running time.
#' @export
#'
#' @examples
#' test <- universal(davis%*%t(davis), upper = function(x)mean(x)+sd(x), lower=function(x)mean(x))
#' test2 <- universal(davis, upper = function(x)mean(x)+2*sd(x), lower = 2, bipartite = TRUE)
#' test3 <- universal(davis, upper = 4, lower = 2, bipartite = TRUE)
universal <- function(M,
                      upper = NULL,
                      lower = NULL,
                      bipartite = NULL,
                      narrative = FALSE){
  #### Argument Checks ####
  if (!(methods::is(upper, "function")) & (!(methods::is(upper, "numeric"))) & (!(methods::is(upper, "NULL")))) {stop("upper must be either function, numeric, or NULL")}
  if (!(methods::is(lower, "function")) & (!(methods::is(lower, "numeric"))) & (!(methods::is(lower, "NULL")))) {stop("lower must be either function, numeric, or NULL")}

  #### Class Conversion ####
  convert <- tomatrix(M)
  class <- convert$summary[[1]]
  M <- convert$G

  #### If not specified, guess whether bipartite or unipartite ####
  if (nrow(M)==ncol(M) & length(bipartite)==0) {
    bipartite <- FALSE
    {warning("The input data is treated as unipartite")}
  }
  if (nrow(M)!=ncol(M) & length(bipartite)==0) {
    bipartite <- TRUE
    {warning("The input data is treated as bipartite")}
  }
  if (length(upper)==0 & length(lower)==0 & bipartite==FALSE) {stop("If the input data is a weighted unipartite graph, then upper and/or lower must be a function or numeric")}

  #### Bipartite Projection ####
  if (bipartite == TRUE){
    if (nrow(M)==ncol(M)) {warning("The input data is square, however you indicated that bipartite = TRUE")}
    if (methods::is(M, "sparseMatrix")) {
      P <- Matrix::tcrossprod(M)
    } else {
      P <- tcrossprod(M)
    }
  } else {
    if (nrow(M)!=ncol(M)) {warning("The input data is rectangular, however you indicated that bipartite = FALSE")}
    P <- M
  }

  #### If both NULL, return the weighted projection ####
  if ((is.null(upper)) & (is.null(lower))){
    backbone <- P
  } else {

    #### For Universal Thresholds ####
    #### Set Threshold Values ####
    if (class(upper) == "function"){
      ut <- upper(P)
    } else {ut <- upper}

    if (class(lower) == "function"){
      lt <- lower(P)
    } else {lt <- lower}

    #### Create Backbone Matrix ####
    backbone <- matrix(0, nrow(P), ncol(P))

    #### Identify negative edges ####
    if (length(lower) > 0){
      negative <- (P<lt)+0
      backbone <- backbone - negative
    }

    #### Identify positive edges ####
    if (length(upper) > 0) {
      positive <- (P>ut)+0
      backbone <- backbone + positive
    }

  } #end else
  diag(backbone) <- 0

  #### Compile Summary ####
  r <- rowSums(M)
  c <- colSums(M)

  a <- c("Model", "Input Class", "Bipartite", "Symmetric", "Weighted", "Number of Rows", "Number of Columns")
  b <- c("Universal Threshold", class[1], convert$summary$bipartite, convert$summary$symmetric, convert$summary$weighted, dim(M)[1], dim(M)[2])

  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Convert to Indicated Class Object ####
  backbone_converted <- frommatrix(backbone, class[1])

  #### Return Backbone and Summary ####
  bb <- list(backbone = backbone_converted, summary = model.summary)
  class(bb) <- "backbone"

  #### Display suggested manuscript text ####
  if (narrative == TRUE) {
    message("Suggested manuscript text and citations:")
    message(" ")
    if (bipartite == FALSE){text <- paste0("From a weighted graph containing ", nrow(M), " agents, we extracted its universal threshold backbone using the backbone package (Domagalski, Neal, & Sagan, 2020).")}
    if ((bipartite == TRUE) & (is.null(upper)) & (is.null(lower))){text <- paste0("From a bipartite graph containing ", nrow(M), " agents and ", ncol(M), " artifacts, we computed the weighted bipartite projection using the backbone package (Domagalski, Neal, & Sagan, 2020).")}
    if ((bipartite == TRUE) & !((is.null(upper)) & (is.null(lower)))){text <- paste0("From a bipartite graph containing ", nrow(M), " agents and ", ncol(M), " artifacts, we computed the weighted bipartite projection, then extracted its universal threshold backbone using the backbone package (Domagalski, Neal, & Sagan, 2020).")}
    if (length(upper) > 0 & length(lower) == 0) {text <- paste0(text, " Edges were retained as positive if their weight was above ", round(ut,3),".")}
    if (length(lower) > 0 & length(upper) == 0) {text <- paste0(text, " Edges were retained as negative if their weight was below ", round(lt,3),".")}
    if (length(lower) > 0 & length(upper) > 0) {text <- paste0(text, " Edges were retained as positive if their weight was above ", round(ut,3),", and as negative if their weight was below ", round(lt,3),".")}
    message(text)
    message("")
    message("Domagalski, R., Neal, Z. P., and Sagan, B. (2020). backbone: An R Package for Backbone Extraction of Weighted Graphs. arXiv:1912.12779 [cs.SI]")
  } #end narrative == TRUE

  return(backbone = bb)
} #end function


