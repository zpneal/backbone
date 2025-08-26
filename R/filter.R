#' Filter based on edge scores
#'
#' @param W a weighted adjacency matrix
#' @param filter string: filter method
#' @param parameter numeric: filtering parameter
#'
#' @return A binary adjacency matrix representing the backbone
#'
#' @noRd
.filter <- function(W, filter, parameter){

  #### Threshold ####
  #Keep edges with weights greater than `parameter`; Depends on edge weight scaling, but larger values keep fewer edges
  if (filter == "threshold") {B <- (W > parameter)*1}
  
  #### Proportion ####
  #Keep strongest `parameter` proportion of edges; 0 = keep 0% of edges, 1 = keep 100% of edges
  if (filter == "proportion") {
    scores <- W[which(W!=0)]  #Vector of non-zero edge scores
    B <- (W >= stats::quantile(scores, probs = (1 - parameter)))*1
  }
  
  #### Degree exponent, from Satuluri et al. (2011) ####
  #Keep edges with neighborhood rank scores at least as small as degree^`parameter`; 0 = keep one edge per node, 1 = keep all edges per node
  if (filter == "degree") {
    A <- (W != 0)*1  #Adjacency matrix
    B <- (W <= (floor(rowSums(A)^parameter)) & W!=0)*1
  }
  
  #### Weighted backbone models ####
  #Parameter functions as alpha: 0 = keep no edges, 1 = keep all edges
  if (filter %in% c("disparity", "lans", "mlf")) {B <- backbone_from_weighted(W, model = filter, alpha = parameter)}

  return(B)
  
}