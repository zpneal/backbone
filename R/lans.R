#' Compute edgewise p-values under the Locally Adaptive Network Sparsification
#'
#' @param A A weighted adjacency matrix
#' @param missing_as_zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone
#'
#' @return
#' If `signed = FALSE` a list containing a matrix of upper-tail p-values.
#'
#' If `signed = TRUE` a list containing a matrix of lower-tail and upper-tail p-values
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references lans: {Foti, N. J., Hughes, J. M., & Rockmore, D. N. (2011). Nonparametric sparsification of complex multiscale networks. *PLOS One, 6*, e16431. \doi{10.1371/journal.pone.0016431}}
#'
#' @noRd
.lans <- function(A, missing_as_zero, signed){

  #### Compute p-values ####
  upper <- matrix(NA, nrow(A), ncol(A))
  if (signed) {lower <- matrix(NA, nrow(A), ncol(A))}
  p_ij <- A / rowSums(A)  #Fractional edge weight from i to j
  for (row in 1:nrow(p_ij)) {upper[row,] <- 1 - unlist(lapply(p_ij[row,], function(i) sum(p_ij[row,] <= i & p_ij[row,]!=0))) / sum(p_ij[row,]!=0)}
  if (signed) {for (row in 1:nrow(p_ij)) {lower[row,] <- 1 - unlist(lapply(p_ij[row,], function(i) sum(p_ij[row,] >= i & p_ij[row,]!=0))) / sum(p_ij[row,]!=0)}}

  if (isSymmetric(A)) {  #If network started as symmetric, backbone should be symmetric
    upper <- pmin(upper,t(upper))  #Use smaller p-value from perspective of both nodes for a given edge
    if (signed) {lower <- pmin(lower,t(lower))}
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[A == 0] <- NA
    if (signed) {lower[A == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
}
