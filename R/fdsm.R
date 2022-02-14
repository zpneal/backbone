#' Extract backbone using the Fixed Degree Sequence Model
#'
#' `fdsm` extracts the backbone of a bipartite projection using the Fixed Degree Sequence Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param trials numeric: The number of bipartite graphs generated to approximate the edge weight distribution.
#' @param method string: The method used to generate random bipartite graphs, one of c("fastball", "curveball")
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
#' @references curveball: {Strona, Giovanni, Domenico Nappo, Francesco Boccacci, Simone Fattorini, and Jesus San-Miguel-Ayanz. 2014. A Fast and Unbiased Procedure to Randomize Ecological Binary Matrices with Fixed Row and Column Totals. *Nature Communications, 5*, 4114. \doi{10.1038/ncomms5114}}
#'
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
#' bb <- fdsm(B, alpha = 0.05, narrative = TRUE, class = "igraph") #An FDSM backbone...
#' plot(bb) #...is sparse with clear communities

fdsm <- function(B, trials = 1000, method = "fastball",
                 alpha = NULL, signed = FALSE, fwer = "none", class = "original", narrative = FALSE,
                 ...){

  #### Argument Checks ####
  if (trials < 0) {stop("trials must be a positive integer")}
  if ((trials > 1) & (trials%%1!=0)) {stop("trials must be decimal < 1, or a positive integer")}

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

  #### Prepare for randomization loop ####
  ### Create Positive and Negative Matrices to hold backbone ###
  rotate <- FALSE  #initialize
  Pupper <- matrix(0, nrow(P), ncol(P))  #Create positive matrix to hold number of times null co-occurence >= P
  Plower <- matrix(0, nrow(P), ncol(P))  #Create negative matrix to hold number of times null co-occurence <] P
  if (nrow(B) > ncol(B)) {  #If B is long, make it wide before randomizing so that randomization is faster
    rotate <- TRUE
    B <- t(B)
  }
  if (method == "fastball") {Bindex <- apply(B==1, 1, which)}  #If using fastball, create an indexed list of 1s
  message("Constructing empirical edgewise p-values -")
  pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)  #Start progress bar

  #### Build Null Models ####
  for (i in 1:trials){

    ### Generate an FDSM Bstar ###
    if (method == "fastball") {Bstar <- fastball(Bindex, nrow(B), ncol(B))}
    if (method == "curveball") {Bstar <- curveball(B)}
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

#' Randomize a binary matrix using the curveball algorithm
#'
#' `curveball` randomizes a binary matrix, preserving the row and column sums
#'
#' @param M matrix: a binary matrix or adjacency list
#' @param R integer: number of rows in `M`
#' @param C integer: number of columns in `M`
#' @param trades integer: number of trades; the default is 5R trades (approx. mixing time)
#'
#' @return
#' matrix: A random binary matrix with same row sums and column sums as M.
#'
#' @details
#' `curveball` is a slightly modified version of Strona et al.'s (2014) R implementation of the curveball
#'    algorithm for generating random binary matrices with fixed row and column sums. For maximum efficiency,
#'    the input should be oriented wide (more columns than rows) rather than long (more rows than columns).
#'    It is also more efficient is the input is supplied not as a matrix but as an adjacency list. Given a
#'    binary matrix `M`, it can be converted into an adjacency list using `apply(M==1, 1, which, simplify = FALSE)`.
#'    When `M` is supplied as an indexed list, `R` and `C` must be specified.
#'
#' @export
#' @references {Strona, Giovanni, Domenico Nappo, Francesco Boccacci, Simone Fattorini, and Jesus San-Miguel-Ayanz. 2014. A Fast and Unbiased Procedure to Randomize Ecological Binary Matrices with Fixed Row and Column Totals. *Nature Communications, 5*, 4114. \doi{10.1038/ncomms5114}}
#'
#' @examples
#' M <- matrix(rbinom(200,1,0.5),10,20)  #A random bipartite graph
#' curveball(M)  #Generation of a random matrix
#' Mlist <- apply(M==1, 1, which, simplify = FALSE)  #Converting M to an adjacency list
#' curveball(Mlist, R = 10, C = 20)  #Faster generation of a random matrix
curveball<-function(M, R = nrow(M), C = ncol(M), trades = 5*R){

  #### Convert to adjacency list (if necessary); Define variables ####
  force(R)  #Evaluate R and C now, setting defaults before matrix gets indexed
  force(C)
  if (methods::is(M, "matrix")) {M <- apply(M==1, 1, which, simplify = FALSE)}
  hp=M
  l_hp=length(hp)

  #### Curveball Swaps ####
  #From here down, verbatim https://static-content.springer.com/esm/art%3A10.1038%2Fncomms5114/MediaObjects/41467_2014_BFncomms5114_MOESM898_ESM.txt
  for (rep in 1:trades){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
  }

  #### Define and Return Random Matrix ####
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

#' Randomize a binary matrix using the fastball algorithm
#'
#' `fastball` randomizes a binary matrix, preserving the row and column sums
#'
#' @param M matrix: a binary matrix (or a list; see details)
#' @param R integer: number of rows in `M`
#' @param C integer: number of columns in `M`
#' @param trades integer: number of trades; the default is 5R trades (approx. mixing time)
#'
#' @return
#' matrix: A random binary matrix with same row sums and column sums as M.
#'
#' @details
#' `fastball` is an optimized C++ implementation of the curveball algorithm (Strona et al.. 2014)
#'    for generating random binary matrices with fixed row and column sums. For maximum efficiency,
#'    the input should be oriented wide (more columns than rows) rather than long (more rows than columns).
#'    It is also more efficient is the input is supplied not as a matrix but as an adjacency list. Given a
#'    binary matrix `M`, it can be converted into an adjacency list using `apply(M==1, 1, which, simplify = FALSE)`.
#'    When `M` is supplied as an indexed list, `R` and `C` must be specified.
#'
#' The fastball algorithm is the fastest known method for randomly sampling binary matrices with given row and
#'    column sums, and is used by [fdsm()] to extract the backbone from a bipartite projection using the fixed
#'    degree sequence model.
#'
#' @references {Godard, Karl and Neal, Zachary P. 2021. fastball: A fast algorithm to sample binary matrices with fixed marginals. \href{https://arxiv.org/abs/2112.04017}{*arXiv:2112.04017*}}
#'
#' @export
#' @examples
#' M <- matrix(rbinom(200,1,0.5),10,20)  #A random 10x20 binary matrix
#' fastball(M)  #Fast generation of a random matrix
#' Mlist <- apply(M==1, 1, which, simplify = FALSE)  #Converting M to a list
#' fastball(Mlist, R = 10, C = 20)  #Even faster generation of a random matrix
fastball <- function(M, R = nrow(M), C = ncol(M), trades = 5*R) {
  force(R)  #Evaluate R and C now, setting defaults before matrix gets indexed
  force(C)
  if (methods::is(M, "matrix")) {M <- apply(M==1, 1, which, simplify = FALSE)}  #If a matrix is provided, convert to an indexed list
  return(fastball_cpp(M, c(R,C), trades))  #Run fastball C++
}

#' Estimate number of Monte Carlo trials needed for FDSM backbone
#'
#' `fdsm.trials` estimates the number of Monte Carlo trials needed to extract an FDSM backbone, correcting for
#'    familywise error rate, and given tolerance for Type-I and Type-II errors
#'
#' @param B graph: An unweighted bipartite graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param type1 numeric: Type-I error used in sample size calculation
#' @param type2 numeric: Type-II error used in sample size calculation
#' @param alpha numeric: Desired Type-I error for tests of edge significance
#' @param fwer boolean: If TRUE, `alpha` is interpreted as the desired familywise error rate.
#'    If FALSE, `alpha` is interpreted as the desired testwise error rate.
#' @param signed boolean: TRUE to estimate the number of trials needed to extract a signed backbone, FALSE to estimate
#'    the number of trials needed to extract a binary backbone
#' @param riskyp numeric: Expected riskiest edge p-value, as a proportion of `alpha` (see details)
#'
#' @details
#' This function uses sample size estimation equations 2.22 and 2.24 given by Fleiss et al. (2013).
#'
#' If `fwer = TRUE`, it assumes that a conservative Bonferroni correction will be used to maintain
#' the familywise error rate across the independent hypothesis tests required for every edge in
#' the bipartite projection of `B`.
#'
#' The required number of trials depends in part on the difference between an edge's estimated p-value and the
#' desired level of statistical significance. If an edge is deemed statistically significant when its p-value is
#' less than 0.05, then there is little risk in making a decision about an edge with an estimated p-value of 0,
#' and fewer trials are required. In contrast, if the edge's estimated p-value is 0.049, there is more risk of
#' making an error and more trials are required. The `riskyp` parameter specifies how close to `alpha` the riskiest
#' expected edge p-value, as a proportion of `alpha`. For example, if `alpha = 0.05` and `riskyp = 0.75`, then
#' the expected riskiest p-value is 0.0375.
#'
#' @references
#' {Fleiss, J. L., Levin, B., & Paik, M. C. (2013). Statistical methods for rates and proportions. John Wiley & Sons.}
#'
#' {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @return integer: estimated minimum number of Monte Carlo trials
#'
#' @export
#'
#' @examples
#' B <- matrix(rbinom(100*1000,1,0.5),100,1000)
#' fdsm.trials(B, riskyp = .75)
fdsm.trials <- function(B, type1 = 0.05, type2 = 0.05, alpha = 0.05, fwer = TRUE, signed = FALSE, riskyp = 0) {
  B <- suppressMessages(tomatrix(B))
  if (B$summary$bipartite == TRUE & B$summary$weighted == FALSE) {B <- B$G} else {stop("B must be a binary bipartite network")}

  if (signed) {  #For a signed backbone...
    tests <- (nrow(B) * (nrow(B) - 1)) / 2  #...every dyad requires an independent test
    p0 <- alpha / 2  #...a two-tailed test is used
    }
  if (!signed) {  #For a binary backbone...
    P <- tcrossprod(B)
    diag(P) <- 0
    tests <- sum((P!=0)*1)/2  #...only dyads with non-zero edge weight are tested
    p0 <- alpha  #...a one-tailed test is used
  }

  if (fwer) {p0 <- p0 / tests} #Bonferroni corrected per-test significance threshold
  p1 <- p0 * riskyp
  n_prime <- ceiling((((stats::qnorm(type1) * sqrt(p0*(1-p0))) + (stats::qnorm(type2) * sqrt(p1*(1-p1))))/(p0 - p1))^2)  #Equation 2.22
  n <- n_prime + (1/(abs(p1 - p0)))  #Equation 2.24
  return(n)
}
