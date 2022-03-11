#' Extract backbone using the Fixed Fill Model
#'
#' `fixedfill` extracts the backbone of a bipartite projection using the Fixed Fill Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros or only ones are automatically removed before computations.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `fixedfill` function compares an edge's observed weight in the projection \eqn{B*t(B)} to the distribution
#'     of weights expected in a projection obtained from a random bipartite graph where the number of edges present
#'     (i.e., the number of cells *filled* with a 1) is equal to the number of edges in B. When B is large, this function
#'     may be impractically slow and may return a backbone object that contains `NaN` values.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected under the null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*
#'
#' @return
#' If `alpha` != NULL: Binary or signed backbone graph of class `class`.
#'
#' If `alpha` == NULL: An S3 backbone object containing three matrices (the weighted graph, edges' upper-tail p-values,
#'    edges' lower-tail p-values), and a string indicating the null model used to compute p-values, from which a backbone
#'    can subsequently be extracted using [backbone.extract()]. The `signed`, `mtc`, `class`, and `narrative` parameters
#'    are ignored.
#'
#' @references {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#' @export
#'
#' @examples
#' #A binary bipartite network of 30 agents & 75 artifacts; agents form three communities
#' B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
#'                  matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.2),10)),
#'            cbind(matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.8),10),
#'                  matrix(rbinom(250,1,.2),10)),
#'            cbind(matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.8),10)))
#'
#' P <- B%*%t(B) #An ordinary weighted projection...
#' plot(igraph::graph_from_adjacency_matrix(P, mode = "undirected",
#'                                          weighted = TRUE, diag = FALSE)) #...is a dense hairball
#'
#' bb <- fixedfill(B, alpha = 0.05, narrative = TRUE, class = "igraph") #A fixedfill backbone...
#' plot(bb) #...is sparse with clear communities

fixedfill <- function(B, alpha = 0.05, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

  #### Argument Checks ####
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Class Conversion ####
  convert <- tomatrix(B)
  if (class == "original") {class <- convert$summary$class}
  B <- convert$G
  if (convert$summary$weighted==TRUE){stop("Graph must be unweighted.")}
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
  }

  #### Bipartite Projection ####
  P <- tcrossprod(B)
  rs <- rowSums(B)

  #### Compute Probabilities ####
  m <- dim(B)[1]
  n <- dim(B)[2]
  f <- sum(B)

  ### Off Diagonal Values ###
  ## This computes log of k! ##
  logsum <- function(k){
    if (k==0){
      return(0)
    }
    return(sum(log(1:k)))
  }

  ## This computes log of (n choose k) ##
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

  maxk <- max(P)  #Largest observed k
  probs <- sapply(0:maxk, FUN = prob_log)  #Probability of observing each k, for 0 <= k <= maxk
  probs <- c(probs, 1 - sum(probs))  #Add one more entry for probability of observing any k > maxk (upper tail of PMF)

  #### Create Backbone Object ####
  Pupper <- apply(P, c(1,2), FUN = function(k)sum(probs[(k+1):(maxk+2)]))  #Sum of probabilities Pij <= k <= maxk and beyond
  Plower <- apply(P, c(1,2), FUN = function(k)sum(probs[1:(k+1)]))  #Sum of probabilities 0 <= k <= Pij
  bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "fixedfill")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = "matrix")
    retained <- round((sum((backbone!=0)*1)) / (sum((P!=0)*1) - nrow(P)),3)*100
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = FALSE, bipartite = TRUE, symmetric = TRUE,
                            signed = signed, mtc = mtc, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "fixedfill", retained = retained)}
    backbone <- frommatrix(backbone, convert = class)
    return(backbone)
  }
}
