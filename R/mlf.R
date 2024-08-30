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
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
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
