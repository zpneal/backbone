#' Extract backbone using the Fixed Degree Sequence Model
#'
#' `fdsm` extracts the backbone of a bipartite projection using the Fixed Degree Sequence Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param trials numeric: the number of bipartite graphs generated to approximate the edge weight distribution. If NULL, the number of trials is selected based on `alpha` (see details)
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' @param progress boolean: TRUE if the progress of Monte Carlo trials should be displayed.
#' @param ... optional arguments
#'
#' @details
#' The `fdsm` function compares an edge's observed weight in the projection \code{B*t(B)} to the distribution of weights
#'    expected in a projection obtained from a random bipartite network where both the row vertex degrees and column
#'    vertex degrees are *exactly* fixed at their values in `B`. It uses the [fastball()] algorithm to generate random
#'    bipartite matrices with give row and column vertex degrees.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger?) is performed for each edge. The resulting backbone
#'    contains edges whose weights are significantly *stronger* than expected in the null model. When `signed = TRUE`, a
#'    two-tailed test (is the weight stronger or weaker?) is performed for each edge. The resulting backbone contains
#'    positive edges for those whose weights are significantly *stronger*, and negative edges for those whose weights are
#'    significantly *weaker*, than expected in the null model.
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
#' If `alpha` == NULL: An S3 backbone object containing (1) the weighted graph as a matrix, (2) upper-tail p-values as a
#'    matrix, (3, if `signed = TRUE`) lower-tail p-values as a matrix, (4, if present) node attributes as a dataframe, and
#'    (5) several properties of the original graph and backbone model, from which a backbone can subsequently be extracted
#'    using [backbone.extract()].
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references fdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports*. \doi{10.1038/s41598-021-03238-3}}
#' @references fastball: {Godard, Karl and Neal, Zachary P. 2022. fastball: A fast algorithm to sample bipartite graphs with fixed degree sequences. *Journal of Complex Networks* \doi{10.1093/comnet/cnac049}}
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
#' bb <- fdsm(B, alpha = 0.05, trials = 1000, narrative = TRUE, class = "igraph") #An FDSM backbone...
#' plot(bb) #...is sparse with clear communities

fdsm <- function(B, alpha = 0.05, trials = NULL, missing.as.zero = FALSE, signed = FALSE, mtc = "none", class = "original", narrative = FALSE, progress = TRUE, ...){

  #### Argument Checks ####
  if (!is.null(trials)) {if (trials < 1 | trials%%1!=0) {stop("trials must be a positive integer")}}
  if (is.null(trials) & is.null(alpha)) {stop("If trials = NULL, then alpha must be specified")}
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Class Conversion ####
  convert <- tomatrix(B)
  if (class == "original") {class <- convert$summary$class}
  attribs <- convert$attribs
  B <- convert$G
  if (convert$summary$weighted==TRUE){stop("Graph must be unweighted.")}
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
  }

  #### Bipartite Projection ####
  P <- tcrossprod(B)

  #### Prepare for randomization loop ####
  ### Create Positive and Negative Matrices to hold backbone ###
  rotate <- FALSE  #initialize
  Pupper <- matrix(0, nrow(P), ncol(P))  #Create positive matrix to hold number of times null co-occurence >= P
  if (signed) {Plower <- matrix(0, nrow(P), ncol(P))}  #Create negative matrix to hold number of times null co-occurence <= P
  if (nrow(B) > ncol(B)) {  #If B is long, make it wide before randomizing so that randomization is faster
    rotate <- TRUE
    B <- t(B)
  }

  #Convert matrix to adjacency list
  if (as.numeric(R.Version()$major)>=4 & as.numeric(R.Version()$minor)>=1) {
    L <- apply(B==1, 1, which, simplify = FALSE)  #Slightly faster, requires R 4.1.0
  } else {
    L <- lapply(asplit(B == 1, 1), which)  #Slightly slower, works for earlier version of R
  }

  #### Build Null Models ####
  if (is.null(trials)) {trials <- trials.needed(M = P, alpha = alpha, signed = signed, missing.as.zero = missing.as.zero, mtc = mtc, ...)}
  if (progress) {message(paste0("Constructing empirical edgewise p-values using ", trials, " trials -" ))}
  if (progress) {pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)}  #Start progress bar
  for (i in 1:trials){

    ### Generate an FDSM Bstar ###
    Lstar <- fastball(L)
    Bstar <- matrix(0,nrow(B),ncol(B))
    for (row in 1:nrow(Bstar)) {Bstar[row,Lstar[[row]]] <- 1L}
    if (rotate) {Bstar <- t(Bstar)}  #If B got rotated from long to wide for randomization, rotate Bstar back from wide to long

    ### Construct Pstar from Bstar ###
    Pstar <- tcrossprod(Bstar)

    ### Check whether Pstar edge is larger/smaller than P edge ###
    Pupper <- Pupper + (Pstar >= P)+0
    if (signed) {Plower <- Plower + (Pstar <= P)+0}

    ### Increment progress bar ###
    if (progress) {utils::setTxtProgressBar(pb, i)}

  } #end for loop
  if (progress) {close(pb)} #End progress bar

  #### Compute p-values ####
  if (rotate) {B <- t(B)}  #If B got rotated from long to wide, rotate B back from wide to long
  Pupper <- (Pupper/trials)
  if (signed) {Plower <- (Plower/trials)}

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing.as.zero) {
    Pupper[P == 0] <- NA
    if (signed) {Plower[P == 0] <- NA}
  }

  #### Create Backbone Object ####
  bb <- list(G = P,  #Preliminary backbone object
             Pupper = Pupper,
             model = "fdsm",
             agents = nrow(B),
             artifacts = ncol(B),
             weighted = FALSE,
             bipartite = TRUE,
             symmetric = TRUE,
             class = class,
             trials = trials)
  if (signed) {bb <- append(bb, list(Plower = Plower))}  #Add lower-tail values, if requested
  if (!is.null(attribs)) {bb <- append(bb, list(attribs = attribs))}  #Add node attributes, if present
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = class, narrative = narrative)
    return(backbone)
  }
} #end fdsm function
