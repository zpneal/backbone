#' Extract backbone using the Fixed Degree Sequence Model
#'
#' `fdsm` extracts the backbone of a bipartite projection using the Fixed Degree Sequence Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param trials numeric: The number of bipartite graphs generated to approximate the edge weight distribution. If NULL, the number of trials is selected to ensure p-values are computed with sufficient power given `alpha`.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param fwer string: type of familywise error rate correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' @param ... optional arguments
#'
#' @details
#' The `fdsm` function compares an edge's observed weight in the projection \code{B*t(B)} to the distribution of weights
#'    expected in a projection obtained from a random bipartite network where both the row vertex degrees and column
#'    vertex degrees are *exactly* fixed at their values in `B`. It uses the [fastball()] algorithm to generate random
#'    bipartite matrices with give row and column vertex degrees. The [fdsm.trials()] function can be used to estimate
#'    the number of random bipartite matrices that must be generated to obtain stable edge *p*-values.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
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
#'    can subsequently be extracted using [backbone.extract()]. The `signed`, `fwer`, `class`, and `narrative` parameters
#'    are ignored.
#'
#' @references fdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports*. \doi{10.1038/s41598-021-03238-3}}
#' @references fastball: {Godard, Karl and Neal, Zachary P. 2022. fastball: A fast algorithm to sample bipartite graphs with fixed degree sequences. \href{https://arxiv.org/abs/2112.04017}{*arXiv:2112.04017*}}#'
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

fdsm <- function(B, trials = NULL, method = "fastball",
                 alpha = NULL, signed = FALSE, fwer = "none", class = "original", narrative = FALSE,
                 ...){

  #### Argument Checks ####
  if (!is.null(trials)) {if (trials < 1 | trials%%1!=0) {stop("trials must be a positive integer")}}
  if (is.null(trials) & is.null(alpha)) {stop("If trials = NULL, then alpha must be specified")}

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

  #### Compute number of trials needed ####
  if (is.null(trials)) {
    trials.alpha <- alpha
    if (signed == TRUE) {trials.alpha <- trials.alpha / 2}  #Two-tailed test
    if (fwer != "none") {  #Adjust trial.alpha using Bonferroni
      if (signed == TRUE) {trials.alpha <- trials.alpha / ((nrow(B)*(nrow(B)-1))/2)}  #Every edge must be tested
      if (signed == FALSE) {trials.alpha <- trials.alpha / (sum(P>0)/2)}  #Every non-zero edge in the projection must be tested
    }
    trials <- ceiling((power.prop.test(p1 = trials.alpha * 0.95, p2 = trials.alpha, sig.level = alpha, power = (1-alpha), alternative = "one.sided")$n)/2)
  }

  #### Prepare for randomization loop ####
  ### Create Positive and Negative Matrices to hold backbone ###
  rotate <- FALSE  #initialize
  Pupper <- matrix(0, nrow(P), ncol(P))  #Create positive matrix to hold number of times null co-occurence >= P
  Plower <- matrix(0, nrow(P), ncol(P))  #Create negative matrix to hold number of times null co-occurence <] P
  if (nrow(B) > ncol(B)) {  #If B is long, make it wide before randomizing so that randomization is faster
    rotate <- TRUE
    B <- t(B)
  }
  Blist <- apply(B==1, 1, which)  #Convert B to an adjacency list
  message(paste0("Constructing empirical edgewise p-values using ", trials, " trials -" ))
  pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)  #Start progress bar

  #### Build Null Models ####
  for (i in 1:trials){

    ### Generate an FDSM Bstar ###
    Bstar <- fastball(Blist, nrow(B), ncol(B))
    if (rotate) {Bstar <- t(Bstar)}  #If B got rotated from long to wide for randomization, rotate Bstar back from wide to long

    ### Construct Pstar from Bstar ###
    Pstar <- tcrossprod(Bstar)

    ### Check whether Pstar edge is larger/smaller than P edge ###
    Pupper <- Pupper + (Pstar >= P)+0
    Plower <- Plower + (Pstar <= P)+0

    ### Increment progress bar ###
    utils::setTxtProgressBar(pb, i)

  } #end for loop
  close(pb) #End progress bar

  #### Create Backbone Object ####
  if (rotate) {B <- t(B)}  #If B got rotated from long to wide, rotate B back from wide to long
  Pupper <- (Pupper/trials)
  Plower <- (Plower/trials)
  bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "fdsm")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, fwer = fwer, class = "matrix")
    retained <- round((sum((backbone!=0)*1)) / (sum((P!=0)*1) - nrow(P)),3)*100
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = FALSE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, fwer = fwer, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = trials, model = "fdsm", retained = retained)}
    backbone <- frommatrix(backbone, convert = class)
    return(backbone)
  }
} #end fdsm function