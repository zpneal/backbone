#' Extract backbone using the Ordinal Stochastic Degree Sequence Model
#'
#' `osdsm` extracts the backbone of a bipartite projection using the Ordinal Stochastic Degree Sequence Model.
#'
#' @param B An ordinally weighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros or only ones are automatically removed before computations.
#' @param trials integer: the number of bipartite graphs generated to approximate the edge weight distribution. If NULL, the number of trials is selected based on `alpha` (see details)
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `osdsm` function compares an edge's observed weight in the projection `B*t(B)` to the distribution of weights
#'    expected in a projection obtained from a random bipartite network where both the rows and the columns contain
#'    approximately the same number of each value. The edges in `B` must be integers, and are assumed to represent an
#'    ordinal-level measure such as a Likert scale.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*
#'
#' The p-values used to evaluate the statistical significance of each edge are computed using Monte Carlo methods. The number of
#'    `trials` performed affects the precision of these p-values, and the confidence that a given p-value is less than the
#'    desired `alpha` level. Because these p-values are proportions (i.e., the proportion of times an edge is weaker/stronger
#'    in the projection of a random bipartite graphs), evaluating the statistical significance of an edge is equivalent to
#'    comparing a proportion (the p-value) to a known proportion (alpha). When `trials = NULL`, the `power.prop.test` function
#'    is used to estimate the required number of trials to make such a comparison with a `alpha` type-I error rate, (1-`alpha`) power,
#'    and when the riskiest p-value being evaluated is at least 5% smaller than `alpha`. When any `mtc` correction is applied,
#'    for simplicity this estimation is based on a conservative Bonferroni correction.
#'
#' @return
#' If `alpha` != NULL: Binary or signed backbone graph of class `class`.
#'
#' If `alpha` == NULL: An S3 backbone object containing three matrices (the weighted graph, edges' upper-tail p-values,
#'    edges' lower-tail p-values), and a string indicating the null model used to compute p-values, from which a backbone
#'    can subsequently be extracted using [backbone.extract()]. The `signed`, `mtc`, `class`, and `narrative` parameters
#'    are ignored.
#'
#' @references {Neal, Z. P. (2017). Well connected compared to what? Rethinking frames of reference in world city network research. *Environment and Planning A, 49*, 2859-2877. \doi{10.1177/0308518X16631339}}
#'
#' @export
#'
#' @examples
#' #A weighted binary bipartite network of 20 agents & 50 artifacts; agents form two communities
#' B <- rbind(cbind(matrix(sample(0:3, 250, replace = TRUE, prob = ((1:4)^2)),10),
#'                  matrix(sample(0:3, 250, replace = TRUE, prob = ((4:1)^2)),10)),
#'            cbind(matrix(sample(0:3, 250, replace = TRUE, prob = ((4:1)^2)),10),
#'                  matrix(sample(0:3, 250, replace = TRUE, prob = ((1:4)^2)),10)))
#'
#' P <- B%*%t(B) #An ordinary weighted projection...
#' plot(igraph::graph_from_adjacency_matrix(P, mode = "undirected",
#'                                          weighted = TRUE, diag = FALSE)) #...is a dense hairball
#'
#' bb <- osdsm(B, alpha = 0.05, narrative = TRUE, class = "igraph", trials = 1000) #An oSDSM backbone...
#' plot(bb) #...is sparse with clear communities

osdsm <- function(B, alpha = 0.05, trials = NULL, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

  #### Class Conversion and Argument Checks ####
  convert <- tomatrix(B)
  if (class == "original") {class <- convert$summary$class}
  B <- convert$G
  if (convert$summary$weighted==FALSE){stop("Graph must be weighted.")}
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
  }
  if (any(B!=as.integer(B)) | any(B < 0)) {stop("Edge weights must be positive integers.")}
  if (!is.null(trials)) {if (trials < 1 | trials%%1!=0) {stop("trials must be a positive integer")}}
  if (is.null(trials) & is.null(alpha)) {stop("If trials = NULL, then alpha must be specified")}
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Bipartite Projection ####
  P <- tcrossprod(B)

  #### Compute number of trials needed ####
  if (is.null(trials)) {
    trials.alpha <- alpha
    if (signed == TRUE) {trials.alpha <- trials.alpha / 2}  #Two-tailed test
    if (mtc != "none") {  #Adjust trial.alpha using Bonferroni
      if (signed == TRUE) {trials.alpha <- trials.alpha / ((nrow(B)*(nrow(B)-1))/2)}  #Every edge must be tested
      if (signed == FALSE) {trials.alpha <- trials.alpha / (sum(P>0)/2)}  #Every non-zero edge in the projection must be tested
    }
    trials <- ceiling((stats::power.prop.test(p1 = trials.alpha * 0.95, p2 = trials.alpha, sig.level = alpha, power = (1-alpha), alternative = "one.sided")$n)/2)
  }

  ### Create Positive and Negative Matrices to hold backbone ###
  Pupper <- matrix(0, nrow(P), ncol(P))
  Plower <- matrix(0, nrow(P), ncol(P))

  #### Compute probabilities for SDSM ####
  #Vectorize the bipartite data
  A <- as.vector(B)                             #Cell values in column 1
  A <- cbind(A, rep(1:nrow(B), times=ncol(B)))  #Row IDs in column 2
  A <- cbind(A, rep(1:ncol(B), times=nrow(B)))  #Column IDs in column 3

  #For each ordinal value, compute and attach row and column sums
  for (i in min(A[,1]):max(A[,1])) {
    if (dim(A)[2]>3) {A[,4] <- ((A[,1]==i)*1)}       #Does cell value equal i, put in column 4
    if (dim(A)[2]==3) {A <- cbind(A, (A[,1]==i)*1)}  #Does cell value equal i, put in column 4
    A <- cbind(A, stats::ave(A[,4],A[,2],FUN=sum))   #Number of times row contains i, put in column 5,7,9...
    A <- cbind(A, stats::ave(A[,4],A[,3],FUN=sum))   #Number of times column contains i, put in column 6,8,10...
  }

  #Compute probabilities using a Proportional Odds Logistic Regression
  A <- as.data.frame(A[,c(1,5:ncol(A))])  #Data frame with cell values and row/column sums
  A[,1] <- factor(A[,1], ordered = TRUE)  #Cell values as ordered factor
  fitted <- suppressWarnings(MASS::polr(A ~ ., data = A)$fitted.values)  #Fitted probabilities (model will be rank-deficient, this is OK)
  probs <- matrix(0, nrow=nrow(fitted), ncol=ncol(fitted))  #Matrix to hold cumulative probabilities
  for (i in 1:ncol(fitted)) {  #Compute cumulative probabilities
    if (i == 1) {probs[,i] <- fitted[,i]}
    if (i > 1) {probs[,i] <- rowSums(fitted[,1:i])}
  }

  #### Build null models ####
  message(paste0("Constructing empirical edgewise p-values using ", trials, " trials -" ))
  pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)  #Start progress bar
  for (i in 1:trials){

    #Start estimation timer; print message
    if (i == 1) {
      start.time <- Sys.time()
      message("Finding the Backbone using Ordinal SDSM")
    }

    #Use probabilities to create an SDSM Bstar
    Bstar <- stats::runif(nrow(probs))                  #Random number
    Bstar <- rowSums((Bstar > probs)*1)                 #Compare to probabilities
    Bstar <- matrix(Bstar, nrow=nrow(B), ncol=ncol(B))  #Convert to matrix

    #Construct Pstar from Bstar, check whether Pstar edge is larger/smaller than P edge
    Pstar <- tcrossprod(Bstar)
    Pupper <- Pupper + (Pstar > P)+0
    Plower <- Plower + (Pstar < P)+0

    ### Increment progress bar ###
    utils::setTxtProgressBar(pb, i)

  } #end for loop
  close(pb) #End progress bar

  #### Create Backbone Object ####
  Pupper <- (Pupper/trials)
  Plower <- (Plower/trials)
  bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "osdsm")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = "matrix")
    retained <- round((sum((backbone!=0)*1)) / (sum((P!=0)*1) - nrow(P)),3)*100
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = TRUE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, mtc = mtc, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = trials, model = "osdsm", retained = retained)}
    backbone <- frommatrix(backbone, convert = class)
    return(backbone)
  }
}

