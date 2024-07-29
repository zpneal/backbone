#' Compute edgewise p-values under the Fixed Fill Model
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
#' @references fixedfill: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @noRd
.fixedfill <- function(I, missing_as_zero, signed){

  #### Define helper functions ####
  # Compute log of k! ##
  logsum <- function(k){
    if (k==0){
      return(0)
    }
    return(sum(log(1:k)))
  }

  # Compute log of (n choose k) ##
  logbinom <- function(n,k){
    if (k == 0){
      return(0)
    }
    else if (k == 1){
      return(log(n))
    }
    else {
      x <- sum(log(n:(n-k+1)))
      y <- sum(log(k:1))
      return(x-y)
    }
  }

  prob_log <- function(k) {
    lb <- max(0, n + k - f)
    ub <- min(n - k, (m - 1) * n + k - f)
    range <- lb:ub
    logvalues <- matrix(0, nrow = 1, ncol = length(range))
    i = 1
    for (r in range){
      logvalues[i] <- (log(2^(n-k-r))+logsum(n)-logsum(k)-logsum(r)-logsum(n-k-r)+logbinom((m-2)*n,f-n-k+r)-logbinom(m*n,f))
      i <- i+1
    }
    return(sum((exp(logvalues))))
  }

  #### Compute bipartite parameters ####
  rs <- rowSums(I)  #Row sums (i.e., agent degrees)
  m <- dim(I)[1]  #Numer of rows (agents)
  n <- dim(I)[2]  #Number of columns (artifacts)
  f <- sum(I)  #Fill (number of edges)

  #### Compute projection parameters ####
  P <- tcrossprod(I)  #Weighted bipartite projection
  max_w <- max(P)  #Largest observed weight w
  probs <- sapply(0:max_w, FUN = prob_log)  #Probability of observing each w, for 0 <= w <= max_w
  probs <- c(probs, 1 - sum(probs))  #Add one more entry for probability of observing any w > max_w (upper tail of PMF)

  #### Compute P-values ####
  upper <- apply(P, c(1,2), FUN = function(k)sum(probs[(k+1):(max_w+2)]))  #Sum of probabilities Pij <= w <= max_w and beyond
  if (signed) {lower <- apply(P, c(1,2), FUN = function(k)sum(probs[1:(k+1)]))}  #Sum of probabilities 0 <= k <= Pij

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
