#' Poisson Binomial distribution computed with Refined Normal Approximation
#'
#' @param kk values where the cdf is to be computed
#' @param pp vector of success probabilities for indicators
#' @param wts the weights for each probability
#'
#' @return cdf, cumulative distribution function
#'
#' @references Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. Computational Statistics & Data Analysis, Vol. 59, pp. 41-51.
#' @details These values are approximated using the Refined Normal Approximation (RNA method).
#'     These functions are originally described by \link[poibin]{ppoibin} and used here under GPL-2 License.
#' @keywords internal
#' @examples
#' \dontrun{probs <- polytope(davis)}
#' \dontrun{P <- davis %*% t(davis)}
#' \dontrun{prob.mat <- matrix(probs, nrow = nrow(davis), ncol = ncol(davis))}
#' \dontrun{prob.imat <- sweep(prob.mat, MARGIN = 2, prob.mat[1,], `*`)}
#' \dontrun{mapply(backbone:::rna, kk= as.data.frame(t(P[1,])), pp = as.data.frame(t(prob.imat)))}
rna <-function(kk,pp,wts=NULL){
  #### Check Arguments ####
  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }
  if(is.null(wts))
  {
    wts=rep(1,length(pp))
  }

  #### Define Variables ####
  pp=rep(pp,wts)
  muk=sum(pp)
  sigmak=sqrt(sum(pp*(1-pp)))
  gammak=sum(pp*(1-pp)*(1-2*pp))
  ind=gammak/(6*sigmak^3)
  kk1=(kk+.5-muk)/sigmak

  #### Compute Statistic and Return ####
  vkk.r=stats::pnorm(kk1)+gammak/(6*sigmak^3)*(1-kk1^2)*stats::dnorm(kk1)
  vkk.r[vkk.r<0]=0
  vkk.r[vkk.r>1]=1
  res=vkk.r
  return(res)
}



