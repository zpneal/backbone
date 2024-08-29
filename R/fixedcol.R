#' Compute edgewise p-values under the Fixed Column Model
#'
#' @param I a binary incidence matrix
#' @param missing_as_zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone
#'
#' @return
#' If `signed = FALSE` a list containing a matrix of upper-tail p-values.
#'
#' If `signed = TRUE` a list containing a matrix of lower-tail and upper-tail p-values
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references fixedcol: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @noRd
.fixedcol <- function(I, missing_as_zero, signed){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Parameters for computing p-values ####
  cs <- colSums(I)
  m <- dim(I)[1]
  p <- ((cs*(cs-1))/(m*(m-1)))

  mu <- sum(p)
  pq <- p*(1-p)
  sigma <- sqrt(sum(pq))
  gamma <- sum(pq*(1-2*p))

  #### Compute p-values using RNA-approximation of poisson-binomial ####
  upper <- ((P-1)+.5-mu)/sigma
  upper <- 1 - (stats::pnorm(upper)+gamma/(6*sigma^3)*(1-upper^2)*stats::dnorm(upper))
  diag(upper) <- NA

  if (signed) {
    lower <- (P+.5-mu)/sigma
    lower <- stats::pnorm(lower)+gamma/(6*sigma^3)*(1-lower^2)*stats::dnorm(lower)
    diag(lower) <- NA
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
