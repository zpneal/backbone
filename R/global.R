#' Extract global threshold backbone
#'
#' @param A A weighted adjacency matrix
#' @param parameter numeric vector of length 1 or 2
#'
#' @return
#' If \code{length(parameter)==1}, an unweighted backbone as a matrix
#' If \code{length(parameter)==2}, a signed backbone as a matrix
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#'
#' @noRd
.global <- function(A, parameter){

  #### Check Parameter ####
  if (!is.numeric(parameter)) {stop("parameter must be a numeric vector of length 1 or 2")}
  if (length(parameter)<1 | length(parameter)>2) {stop("parameter must be a numeric vector of length 1 or 2")}
  
  #### Apply Global Thresholds ####
  backbone <- matrix(NA, nrow(A), ncol(A))  #Start with empty adjacency matrix
  backbone[which(A > max(parameter))] <- 1  #Preserve edges above upper threshold
  if (length(parameter)==2) {backbone[which(A < min(parameter))] <- -1}  #Optionally, preserve edges below lower threshold
  backbone[which(is.na(backbone))] <- 0  #Fill in missing edges
  return(backbone)
  }
