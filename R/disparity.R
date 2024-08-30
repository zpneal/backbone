#' Compute edgewise p-values under the Disparity Filter
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
#' @references disparity filter: {Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. *Proceedings of the National Academy of Sciences, 106*, 6483-6488. \doi{10.1073/pnas.0808904106}}
#'
#' @noRd
.disparity <- function(A, missing_as_zero, signed){

  #### Set Parameters ####
  strength <- rowSums(A)
  binary <- (A>0)+0
  degree <- rowSums(binary)

  #### Compute p-values ####
  if (isSymmetric(A)) {
    P <- A/strength
    pvalues <- (1-P)^(degree-1)
    upper <- as.matrix(pvalues)      #Asymmetric p-values, one from the perspective of each node
    upper <- pmin(upper,t(upper))    #From Serrano: "satisfy the above criterion for at least one of the two nodes"
    if (signed) {lower <- 1-upper}
  }

  if (!isSymmetric(A)) {
    outp <- A/strength
    outvalues <- (1-outp)^(degree-1)
    inp <- t(A)/(colSums(A))
    invalues <- t((1-inp)^(colSums(binary)-1))
    upper <- pmin(invalues,outvalues)
    if (signed) {lower <- 1-upper}
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[A == 0] <- NA
    if (signed) {lower[A == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
