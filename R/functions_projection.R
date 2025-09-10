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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
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
  diag(upper) <- NA

  if (signed) {
    lower <- apply(P, c(1,2), FUN = function(k)sum(probs[1:(k+1)]))  #Sum of probabilities 0 <= k <= Pij
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

#' Compute edgewise p-values under the Fixed Row Model
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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references fixedrow: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @noRd
.fixedrow <- function(I, missing_as_zero, signed){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Prepare dyad list ####
  df <- data.frame(row = row(P)[upper.tri(P)],            #Dataframe of dyads in upper triangle
                   col = col(P)[upper.tri(P)],
                   weight = as.vector(P[upper.tri(P)]))

  rs <- rowSums(I)  #Find row sums in bipartite (agent degrees)
  df$row_sum_i <- rs[df$row]
  df$row_sum_j <- rs[df$col]
  df$diff <- ncol(I)-df$row_sum_i  #Difference in total number of artifacts and i's degree

  if (!missing_as_zero) {df$weight[which(df$weight==0)] <- NA}  #If missing edges should not be tested, replace weight with NA

  #### Compute p-values ####
  df$upper <- stats::phyper(df$weight-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)
  upper <- matrix(NA, nrow = nrow(P), ncol = nrow(P))  #Start with empty matrix of upper-tail p-values
  upper[upper.tri(upper)] <- df$upper  #Insert upper-tail p-values
  upper[lower.tri(upper)] = t(upper)[lower.tri(upper)]  #Make symmetric

  if (signed) {
    df$lower <- stats::phyper(df$weight, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)
    lower <- matrix(NA, nrow = nrow(P), ncol = nrow(P))  #Start with empty matrix of upper-tail p-values
    lower[upper.tri(lower)] <- df$lower
    lower[lower.tri(lower)] = t(lower)[lower.tri(lower)]
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }

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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
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

#' Compute edgewise p-values under the Stochastic Degree Sequence Model
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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references sdsm: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @references sdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @noRd
.sdsm <- function(I, missing_as_zero, signed){

  P <- tcrossprod(I)  #Weighted bipartite projection

  probs <- bicm(I)  #Bipartite edge probabilities under BiCM
  probs <- lapply(seq_len(nrow(probs)), function(i) probs[i,])  #Store probabilities as list

  #### Compute p-values ####
  upper <- matrix(NA, nrow(P), ncol(P))                #Set upper-tail p-value to NA initially, untested edges have p = NA
  if (signed) {lower <- matrix(NA, nrow(P), ncol(P))}  #If signed, set lower-tail p-value to NA initially

  for (col in 1:(ncol(P)-1)) {  #Loop over lower triangle of projection
    for (row in (col+1):nrow(P)) {

      if (missing_as_zero) {  #If missing edges should be treated as zero, test each one
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

      if (!missing_as_zero & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

    }
  }
  upper[upper.tri(upper)] <- t(upper)[upper.tri(upper)]  #Make upper-tail p-value matrix symmetric
  if (signed) {lower[upper.tri(lower)] <- t(lower)[upper.tri(lower)]}  #Make lower-tail p-value matrix symmetric

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }

#' Compute edgewise p-values under the Stochastic Degree Sequence Model with Edge Constraints
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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references sdsm-ec model: {Neal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. *International Conference on Complex Networks and Their Applications, 12*, 127-136. \doi{10.1007/978-3-031-53468-3_11}}
#'
#' @noRd
.sdsm_ec <- function(I, missing_as_zero, signed){

  #### Construct weighted projection ####
  I_unweighted <- I
  I_unweighted[I_unweighted==10] <- 0  #Make structural 0s ordinary 0
  I_unweighted[I_unweighted==11] <- 1  #Make structural 1s ordinary 1
  P <- tcrossprod(I_unweighted)  #Projection, not considering any structural 0s or 1s

  #### Compute probabilities with edge constraints using Logit ####
  # Prepare dyad list
  A <- data.frame(edge = as.vector(I),     #Data frame of bipartite dyads
                  row = as.vector(row(I)),
                  col = as.vector(col(I)))
  A$edge2 <- A$edge  #Copy of edges
  A$edge2[which(A$edge>1)] <- 0  #Set structural edges to 0 so they're not considered in marginals
  A$rowmarg <- stats::ave(A$edge2,A$row,FUN=sum,na.rm=TRUE)  #Compute rowsums (agent degree), excluding structural edges
  A$colmarg <- stats::ave(A$edge2,A$col,FUN=sum,na.rm=TRUE)  #Compute colsums (artifact degree), excluding structural edges

  #Compute probabilities on non-structural edges using logit
  A$edge2[which(A$edge>1)] <- NA  #Set structural edges to NA so they're not considered in marginals
  model.estimates <- suppressWarnings(stats::glm(formula = edge2 ~ rowmarg + colmarg, family = stats::binomial(link="logit"), data=A))
  A$probs <- as.vector(suppressWarnings(stats::predict(model.estimates, newdata = A, type = "response")))

  #Insert structural probabilities
  A$probs[which(A$edge==10)] <- 0  #Structural zeros have probability = 0
  A$probs[which(A$edge==11)] <- 1  #Structural ones have probability = 1

  #Probability matrix
  probs <- matrix(A$probs, nrow = nrow(I), ncol = ncol(I))  #Probability matrix
  probs <- lapply(seq_len(nrow(probs)), function(i) probs[i,])  #Store probabilities as list

  #### Compute p-values ####
  upper <- matrix(NA, nrow(P), ncol(P))                #Set upper-tail p-value to NA initially, untested edges have p = NA
  if (signed) {lower <- matrix(NA, nrow(P), ncol(P))}  #If signed, set lower-tail p-value to NA initially

  for (col in 1:(ncol(P)-1)) {  #Loop over lower triangle of projection
    for (row in (col+1):nrow(P)) {

      if (missing_as_zero) {  #If missing edges should be treated as zero, test each one
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

      if (!missing_as_zero & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

    }
  }
  upper[upper.tri(upper)] <- t(upper)[upper.tri(upper)]  #Make upper-tail p-value matrix symmetric
  if (signed) {lower[upper.tri(lower)] <- t(lower)[upper.tri(lower)]}  #Make lower-tail p-value matrix symmetric

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }

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
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references fdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#' @references fastball: {Godard, K. and Neal, Z. P. (2022). fastball: A fast algorithm to randomly sample bipartite graphs with fixed degree sequences. *Journal of Complex Networks, 10*, cnac049. \doi{10.1093/comnet/cnac049}}
#'
#' @noRd
.fdsm <- function(I, missing_as_zero, signed, alpha, mtc, trials){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Compute number of trials required ####
  if (is.null(trials)) {
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
  diag(upper) <- NA

  if (signed) {
    lower <- (lower / trials)
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

#' Compute Poisson binomial distribution function using the refined normal approximation
#'
#' @param k numeric: value where the pdf should be evaluated
#' @param p vector: vector of success probabilities
#' @param lowertail boolean: If TRUE return both upper & lower tail probabilities,
#'    if FALSE return only upper tail probability
#'
#' @return vector, length 2: The first value (if lower = TRUE) is the lower tail probability, the
#'    probability of observing `k` or fewer successes when each trial has probability `p` of success.
#'    The second value is the upper tail probability, the probability of observing `k` or more
#'    successes when each trial has probability `p` of success.
#'
#' @references
#' {Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. *Computational Statistics and Data Analysis, 59*, 41-51. \doi{10.1016/j.csda.2012.10.006}}
#'
#' @noRd
.pb <-function(k, p, lowertail) {
  #Compute parameters
  mu <- sum(p)
  pq <- p*(1-p)
  sigma <- sqrt(sum(pq))
  gamma <- sum(pq*(1-2*p))

  #Lower tail p-value, if requested
  if (lowertail) {
    x <- (k+.5-mu)/sigma
    lower <- stats::pnorm(x)+gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x)
  } else {lower <- NA}

  #Upper tail p-value
  x <- ((k-1)+.5-mu)/sigma
  upper <- stats::pnorm(x,lower.tail=F)-gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x)

  #Combine, truncate, return
  prob <- c(lower,upper)
  prob[prob<0] <- 0
  prob[prob>1] <- 1
  return(prob)
}
