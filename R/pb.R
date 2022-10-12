#' Poisson binomial distribution function
#'
#' `pb` computes the poisson binomial distribution function using the refined normal approximation.
#'
#' @param k numeric: value where the pdf should be evaluated
#' @param p vector: vector of success probabilities
#' @param lowertail boolean: If TRUE return both upper & lower tail probabilities,
#'    if FALSE return only upper tail probability
#'
#' @details
#' The Refined Normal Approximation (RNA) offers a close approximation when `length(p)` is
#'    large (Hong, 2013). This function is based on `ppoibin()` from the `poibin` package.
#'
#' @return vector, length 2: The first value (if lower = TRUE) is the lower tail probability, the
#'    probability of observing `k` or fewer successes when each trial has probability `p` of success.
#'    The second value is the upper tail probability, the probability of observing `k` or more
#'    successes when each trial has probability `p` of success.
#'
#' @references
#' {Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. *Computational Statistics and Data Analysis, 59*, 41-51. \doi{10.1016/j.csda.2012.10.006}}
#'
#' @export
#'
#' @examples
#' pb(50,runif(100))
pb <-function(k, p, lowertail=TRUE) {

  #Compute parameters
  mu <- sum(p)
  sigma <- sqrt(sum(p*(1-p)))
  gamma <- sum(p*(1-p)*(1-2*p))

  #Lower tail p-value, if requested
  if (lowertail) {
    x <- (k+.5-mu)/sigma
    lower <- stats::pnorm(x)+gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x)
  } else {lower <- NA}

  #Upper tail p-value
  x <- (k-1+.5-mu)/sigma
  upper <- 1 - (stats::pnorm(x)+gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x))

  #Combine, truncate, return
  prob <- c(lower,upper)
  prob[prob<0] <- 0
  prob[prob>1] <- 1
  return(prob)
}