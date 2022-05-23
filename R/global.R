#' Compute global threshold backbone
#'
#' `global` extracts the backbone of a weighted network using a global threshold
#'
#' @param W A weighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}, or dataframe; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object.
#' @param upper real, FUN, or NULL: upper threshold value or function that evaluates to an upper threshold value.
#' @param lower real, FUN, or NULL: lower threshold value or function that evaluates to a lower threshold value.
#' @param keepzeros boolean: TRUE if zero-weight edges in `W` should be excluded from (i.e. also be zero in) the backbone
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `W`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `global` function retains a edge in the backbone if its weight exceeds `upper`. If a `lower` threshold is also
#'    specified, it returns a signed backbone in which edge weights are set to 1 if above the given upper threshold,
#'    set to -1 if below the given lower threshold, and set to 0 otherwise. Prior to v2.0.0, the `global` function was
#'    called [universal()], which is now depricated.
#'
#' If `W` is an unweighted bipartite graph, any rows and columns that contain only zeros or only ones are removed, then
#'    the global threshold is applied to its weighted bipartite projection.
#'
#' @return Binary or signed backbone graph of class given in parameter `class`.
#'
#' @references {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *arXiv:2203.11055 \[cs.SI\]*. \doi{10.48550/arXiv.2203.11055}}
#' @references {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @export
#'
#' @examples
#' W <- matrix(sample(0:5, 100, replace = TRUE), 10) #Random weighted graph
#' diag(W) <- 0
#' W
#' global(W, narrative = TRUE)  #Keep all non-zero edges
#' global(W, upper = 4, lower = 2, narrative = TRUE)  #Signed with specified thresholds
#' global(W, upper = function(x)mean(x),  #Above-average --> positive edges
#'           lower = function(x)mean(x), narrative = TRUE)  #Below-average --> negative edges
global <- function(W, upper = 0, lower = NULL, keepzeros = TRUE, class = "original", narrative = FALSE){

  #### Argument Checks ####
  if (!(methods::is(upper, "function")) & (!(methods::is(upper, "numeric"))) & (!(methods::is(upper, "NULL")))) {stop("upper must be either function, numeric, or NULL")}
  if (!(methods::is(lower, "function")) & (!(methods::is(lower, "numeric"))) & (!(methods::is(lower, "NULL")))) {stop("lower must be either function, numeric, or NULL")}
  if (is.null(upper) & is.null(lower)) {stop("upper and lower cannot both be NULL")}

  #### Class Conversion ####
  convert <- tomatrix(W)
  if (class == "original") {class <- convert$summary$class}
  attribs <- convert$attribs
  M <- convert$G
  if (convert$summary$bipartite) {
    artifacts <- ncol(M)      #Record the number of artifacts
    M <- tcrossprod(M)  #Convert bipartite to weighted projection
    diag(M) <- 0        #Fill diagonal with zeroes
  }
  symmetric <- convert$summary$symmetric
  signed <- !is.null(lower)

  #### Set Threshold Values ####
  if (methods::is(upper, "function")){ut <- upper(M)} else {ut <- upper}
  if (methods::is(lower, "function")){lt <- lower(M)} else {lt <- lower}

  #### Apply Global Thresholds ####
  if (!is.null(upper)) {positive <- M > ut} else {positive <- 0}
  if (!is.null(lower)) {negative <- M < lt} else {negative <- 0}
  backbone <- positive + (negative * -1)

  #### Restore Zeroes ####
  if (keepzeros == TRUE) {backbone[which(M==0)] <- 0}

  #### Display suggested manuscript text ####
  if (narrative == TRUE) {
    reduced_edges <- round((sum(M!=0) - sum(backbone!=0)) / sum(M!=0),3)*100  #Percent decrease in number of edges
    reduced_nodes <- round((max(sum(rowSums(M)!=0),sum(colSums(M)!=0)) - max(sum(rowSums(backbone)!=0),sum(colSums(backbone)!=0))) / max(sum(rowSums(M)!=0),sum(colSums(M)!=0)),3) * 100  #Percent decrease in number of connected nodes
    if (narrative == TRUE) {write.narrative(agents = nrow(M), artifacts = NULL, weighted = TRUE, bipartite = FALSE, symmetric = TRUE,
                                            signed = signed, mtc = "none", alpha = NULL, s = NULL, ut = ut, lt = lt, trials = NULL, model = "global",
                                            reduced_edges = reduced_edges, reduced_nodes = reduced_nodes)}
  }

  backbone <- frommatrix(backbone, attribs, convert = class)
  return(backbone)
}

#' Wrapper for global()
#' @param W A weighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}, or dataframe; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object.
#' @param upper Real, FUN, or NULL: upper threshold value or function that evaluates to an upper threshold value.
#' @param lower Real, FUN, or NULL: lower threshold value or function that evaluates to a lower threshold value.
#' @param keepzeros Boolean: TRUE if zero-weight edges in `M` should be missing in the backbone
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative Boolean: TRUE if suggested text for a manuscript is to be returned
#' @export
universal <- function(W, upper = 0, lower = NULL, keepzeros = TRUE, class = "original", narrative = FALSE){
  warning("The universal() function is now called global(); please use global() instead.")
  global(W, upper = 0, lower = NULL, keepzeros = TRUE, class = "original", narrative = FALSE)
}
