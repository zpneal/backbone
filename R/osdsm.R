#' Extract backbone using the Ordinal Stochastic Degree Sequence Model
#'
#' `osdsm` extracts the backbone of a bipartite projection using the Ordinal Stochastic Degree Sequence Model.
#'
#' @param B An ordinally weighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object.
#' @param trials integer: the number of bipartite graphs generated to approximate the edge weight distribution. If NULL, the number of trials is selected based on `alpha` (see details)
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `osdsm` function compares an edge's observed weight in the projection `B*t(B)` to the distribution of weights
#'    expected in a projection obtained from a random bipartite network where both the rows and the columns contain
#'    approximately the same number of each value. The edges in `B` must be integers, and are assumed to represent an
#'    ordinal-level measure such as a Likert scale that starts at 0.
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
#' If `alpha` == NULL: An S3 backbone object containing (1) the weighted graph as a matrix, (2) upper-tail p-values as a
#'    matrix, (3, if `signed = TRUE`) lower-tail p-values as a matrix, and (4) a string indicating the null model used to
#'    compute p-values, from which a backbone can subsequently be extracted using [backbone.extract()]. The `mtc`, `class`,
#'    and `narrative` parameters are ignored.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references osdsm: {Neal, Z. P. (2017). Well connected compared to what? Rethinking frames of reference in world city network research. *Environment and Planning A, 49*, 2859-2877. \doi{10.1177/0308518X16631339}}
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
#' bb <- osdsm(B, alpha = 0.05, narrative = TRUE,  #An oSDSM backbone...
#'             class = "igraph", trials = 100)
#' plot(bb) #...is sparse with clear communities

osdsm <- function(B, alpha = 0.05, trials = NULL, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

  #### Class Conversion and Argument Checks ####
  convert <- tomatrix(B)
  if (class == "original") {class <- convert$summary$class}
  attribs <- convert$attribs
  B <- convert$G
  if (convert$summary$weighted==FALSE){stop("Graph must be weighted.")}
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
  }
  if (!is.null(trials)) {if (trials < 1 | trials%%1!=0) {stop("trials must be a positive integer")}}
  if (is.null(trials) & is.null(alpha)) {stop("If trials = NULL, then alpha must be specified")}
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}
  weights <- sort(unique(as.vector(B)))  #Unique edge weights
  if (sum(weights%%1)!=0) {stop("Edge weights must be integers that reflect values on an ordinal scale.")}
  if (any(weights < 0)) {stop("Edge weights must be positive integers that reflect values on an ordinal scale.")}
  if (!(0 %in% weights)) {warning("Although no edges have weight = 0, the weight scale is assumed to start at 0.")}
  weights <- c(0:max(B))  #If weights are valid, this is the full scale

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
  if (signed) {Plower <- matrix(0, nrow(P), ncol(P))}

  #### Compute probabilities for SDSM ####
  #Vectorize the bipartite data
  A <- data.frame(value = as.vector(B))     #Edge weight
  A$rowid <- rep(1:nrow(B), times=ncol(B))  #Row index
  A$colid <- rep(1:ncol(B), each=nrow(B))   #Column index

  #Compute conditional probabilities using logistic regression (see Neal, 2017)
  for (value in 1:max(weights)) {  #For each edge weight > 0
    dat <- data.frame(y = (A$value>=value)*1, x1 = stats::ave(A$value>=value,A$rowid,FUN=sum), x2 = stats::ave(A$value>=value,A$colid, FUN=sum))
    fitted <- suppressWarnings(stats::glm(y ~ x1 + x2, data = dat[which(A$value>=(value-1)),], family = "binomial"))
    A <- cbind(A, stats::predict(fitted, newdata = dat, type = "response"))
  }

  #Transform into unconditional probabilities (see Neal, 2017)
  for (value in c(1:max(weights),0)) {  #For each edge weight value, doing weight = 0 last
    if (value == 1) {ucondp <- A[,4]}
    if (value > 1) {ucondp <- apply(A[,4:(value+3)], 1, prod)}
    if (value != 0 & value != max(weights)) {ucondp <- ucondp * (1 - A[,(value+4)])}
    if (value == 0) {ucondp <- 1 - rowSums(A[,(4+max(weights)):ncol(A)])}
    A <- cbind(A, ucondp)
  }
  A <- A[,c(1:3,(4+max(weights)):ncol(A))]
  colnames(A) <- c("value", "rowid", "colid", paste0("p", c(1:max(weights),0)))
  A$rand <- NA

  #### Build null models ####
  message(paste0("Constructing empirical edgewise p-values using ", trials, " trials -" ))
  pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)  #Start progress bar
  for (i in 1:trials){

    #Use probabilities to create an SDSM Bstar
    A$rand <- apply(X = A[,4:(max(weights)+4)], MARGIN = 1, FUN = function(x) sample(c(1:max(weights),0), size = 1, replace= TRUE, prob = x))
    Bstar <- matrix(A$rand, nrow=nrow(B), ncol=ncol(B))  #Convert to matrix

    #Construct Pstar from Bstar, check whether Pstar edge is larger/smaller than P edge
    Pstar <- tcrossprod(Bstar)
    Pupper <- Pupper + (Pstar > P)+0
    if (signed) {Plower <- Plower + (Pstar < P)+0}

    #Increment progress bar
    utils::setTxtProgressBar(pb, i)

  } #end for loop
  close(pb) #End progress bar

  #### Create Backbone Object ####
  Pupper <- (Pupper/trials)
  if (signed) {Plower <- (Plower/trials)}
  if (signed) {bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "osdsm")} else {bb <- list(G = P, Pupper = Pupper, model = "osdsm")}
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = "matrix")
    reduced_edges <- round(((sum(P!=0)-nrow(P)) - sum(backbone!=0)) / (sum(P!=0)-nrow(P)),3)*100  #Percent decrease in number of edges
    reduced_nodes <- round((max(sum(rowSums(P)!=0),sum(colSums(P)!=0)) - max(sum(rowSums(backbone)!=0),sum(colSums(backbone)!=0))) / max(sum(rowSums(P)!=0),sum(colSums(P)!=0)),3) * 100  #Percent decrease in number of connected nodes
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = TRUE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, mtc = mtc, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "osdsm",
                                            reduced_edges = reduced_edges, reduced_nodes = reduced_nodes)}
    backbone <- frommatrix(backbone, attribs, convert = class)
    return(backbone)
  }
}

