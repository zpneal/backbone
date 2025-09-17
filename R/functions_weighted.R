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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
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
    upper <- pmin(pvalues,t(pvalues))    #From Serrano: "satisfy the above criterion for at least one of the two nodes"
    if (signed) {lower <- 1-upper}
  } else {
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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
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

#' Compute edgewise p-values under the Marginal Likelihood Filter
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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references mlf: {Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. *Physical Review E, 93*, 012304. \doi{10.1103/PhysRevE.93.012304}}
#'
#' @noRd
.mlf <- function(A, missing_as_zero, signed){

  #### Compute p-values ####
  if (isSymmetric(A)) {
    upper <- matrix(NA, nrow(A), ncol(A))
    if (signed) {lower <- matrix(NA, nrow(A), ncol(A))}
    T <- sum(rowSums(A))/2
    p <- (rowSums(A) %*% t(rowSums(A))) / (2 * (T^2))
    for (col in 1:ncol(A)) {  #Loop over lower triangle
      for (row in col:nrow(A)) {

        if (missing_as_zero) {  #If missing edges should be treated as zero, test each one
          upper[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {lower[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "less")$p.value}
        }

        if (!missing_as_zero & A[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
          upper[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {lower[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "less")$p.value}
        }

      }
    }
    upper[upper.tri(upper)] <- t(upper)[upper.tri(upper)]  #Add upper triangle
    if (signed) {lower[upper.tri(lower)] <- t(lower)[upper.tri(lower)]}
  }

  if (!isSymmetric(A)) {
    upper <- matrix(NA, nrow(A), ncol(A))
    if (signed) {lower <- matrix(NA, nrow(A), ncol(A))}
    T <- sum(rowSums(A))
    p <- (rowSums(A) %*% t(colSums(A))) / (T^2)
    for (col in 1:ncol(A)) {  #Loop over full matrix
      for (row in 1:nrow(A)) {

        if (missing_as_zero) {  #If missing edges should be treated as zero, test each one
          upper[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {lower[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "less")$p.value}
        }

        if (!missing_as_zero & A[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
          upper[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {lower[row,col] <- stats::binom.test(A[row,col], T, p[row,col], alternative = "less")$p.value}
        }

      }
    }
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
}

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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
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
