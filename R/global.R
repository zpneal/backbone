#' Extract global threshold backbone
#'
#' @param A A weighted adjacency matrix
#' @param missing_as_zero boolean: treat missing edges as edges with zero weight and consider them for inclusion/exclusion in backbone
#' @param parameter numeric vector of length 1 or 2
#'
#' @return
#' If \code{length(parameter)==1}, an unweighted backbone as a matrix
#' If \code{length(parameter)==2}, a signed backbone as a matrix
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#'
#' @noRd
.global <- function(A, missing_as_zero, parameter){

  #### Apply Global Thresholds ####
  backbone <- matrix(NA, nrow(A), ncol(A))  #Start with empty adjacency matrix

  if (missing_as_zero) {  #Evaluate all edges
    backbone[A > max(parameter)] <- 1  #Preserve edges above upper threshold
    if (length(parameter)==2) {backbone[A < min(parameter)] <- -1}  #Optionally, preserve edges below lower threshold
    backbone[is.na(backbone)] <- 0  #Fill in missing edges
  }

  if (!missing_as_zero) {  #Evaluate non-zero edges
    backbone[A > max(parameter) & A!=0] <- 1  #Preserve edges above upper threshold
    if (length(parameter)==2) {backbone[A < min(parameter) & A!=0] <- -1}  #Optionally, preserve edges below lower threshold
    backbone[is.na(backbone)] <- 0  #Fill in missing edges
  }

  return(backbone)
  }
