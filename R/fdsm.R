#' Compute edgewise p-values under the Fixed Degree Sequence Model
#'
#' @param I a binary incidence matrix
#' @param missing_as_zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone
#' @param alpha real: significance level of hypothesis test(s)
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param trials numeric: number of bipartite graphs generated using fastball to approximate the edge weight distribution
#'
#' @return
#' If `signed = FALSE` a list containing a matrix of upper-tail p-values.
#'
#' If `signed = TRUE` a list containing a matrix of lower-tail and upper-tail p-values
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references fdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#' @references fastball: {Godard, K. and Neal, Z. P. (2022). fastball: A fast algorithm to randomly sample bipartite graphs with fixed degree sequences. *Journal of Complex Networks, 10*, cnac049. \doi{10.1093/comnet/cnac049}}
#'
#' @noRd
.fdsm <- function(I, missing_as_zero, signed, alpha, mtc, trials = NULL){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Compute number of trials required ####
  if (trials == 0) {
    if (signed == TRUE) {test.alpha <- alpha / 2} else {test.alpha <- alpha}  #Adjust alpha if two-tailed test
    if (mtc != "none") {  #If multiple test correction is requested, conservatively adjust alpha using Bonferroni
      if (!missing_as_zero) {tests <- sum(lower.tri(P) & P!=0)}  #Non-zero entries in lower triangle
      if (missing_as_zero) {tests <- sum(lower.tri(P))}  #Entries in lower triangle
      test.alpha <- test.alpha / tests
    }
    #p1 = A hypothetical empirical monte carlo p-value we want to evaluate that is close to (within alpha percent of) the alpha level
    #p2 = The alpha level against which we are evaluating p1, with any two-tailed or mtc adjustments
    #Because type-I errors (a false edge is included in the backbone) is as bad as type-II errors (a true edge is omitted from the backbone), therefore power = alpha
    trials <- ceiling((stats::power.prop.test(p1 = test.alpha * (1 - alpha), p2 = test.alpha, sig.level = alpha, power = (1-alpha), alternative = "one.sided")$n)/2)
  }

  #### Prepare for randomization loop ####
  ### Create Positive and Negative Matrices to hold backbone ###
  rotate <- FALSE  #initialize
  upper <- matrix(0, nrow(P), ncol(P))  #Create positive matrix to hold number of times null co-occurence >= P
  if (signed) {lower <- matrix(0, nrow(P), ncol(P))}  #Create negative matrix to hold number of times null co-occurence <= P
  if (nrow(I) > ncol(I)) {  #If I is long, make it wide before randomizing so that randomization is faster
    rotate <- TRUE
    I <- t(I)
  }

  #Convert matrix to adjacency list
  if (as.numeric(R.Version()$major)>=4 & as.numeric(R.Version()$minor)>=1) {
    L <- apply(I == 1, 1, which, simplify = FALSE)  #Slightly faster, requires R 4.1.0
  } else {
    L <- lapply(asplit(I == 1, 1), which)  #Slightly slower, works for earlier version of R
  }

  #### Build Null Models ####
  message(paste0("Constructing edges' Monte Carlo p-values" ))
  pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)  #Start progress bar
  for (i in 1:trials){

    ### Generate an FDSM Bstar ###
    Lstar <- fastball(L)
    Istar <- matrix(0,nrow(I),ncol(I))
    for (row in 1:nrow(Istar)) {Istar[row,Lstar[[row]]] <- 1L}

    ### Construct Pstar from Istar ###
    if (rotate) {Pstar <- crossprod(Istar)}  #If I *was* rotated, generate projection on columns
    if (!rotate) {Pstar <- tcrossprod(Istar)}  #If I *was* not rotated, generate projection on rows

    ### Check whether Pstar edge is larger/smaller than P edge ###
    upper <- upper + (Pstar >= P) + 0
    if (signed) {lower <- lower + (Pstar <= P) + 0}

    ### Increment progress bar ###
    utils::setTxtProgressBar(pb, i)

  } #end for loop
  close(pb) #End progress bar

  #### Compute p-values ####
  upper <- (upper / trials)
  if (signed) {lower <- (lower / trials)}

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
