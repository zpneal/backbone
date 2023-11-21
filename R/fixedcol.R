#' Extract backbone using the Fixed Column Model
#'
#' `fixedcol` extracts the backbone of a bipartite projection using the Fixed Column Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.

#' @details
#' This `fixedcol` function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the *column* vertex degrees are fixed but the row vertex degrees are allowed to vary.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger?) is performed for each edge. The resulting backbone
#'    contains edges whose weights are significantly *stronger* than expected in the null model. When `signed = TRUE`, a
#'    two-tailed test (is the weight stronger or weaker?) is performed for each edge. The resulting backbone contains
#'    positive edges for those whose weights are significantly *stronger*, and negative edges for those whose weights are
#'    significantly *weaker*, than expected in the null model.
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
#' @references fixedcol: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
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
#' bb <- fixedcol(B, alpha = 0.05, narrative = TRUE, class = "igraph") #A fixedcol backbone...
#' plot(bb) #...is sparse with clear communities

fixedcol <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

  #### Argument Checks ####
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

  #### Parameters for computing p-values ####
  cs <- colSums(B)
  m = dim(B)[1]
  pt <- ((cs*(cs-1))/(m*(m-1)))
  mu <- sum(pt)
  sigma <- sqrt(sum(pt*(1-pt)))
  gamma <- sum(pt*(1-pt)*(1-2*pt))

  #### Compute p-values using RNA-approximation of poisson-binomial ####
  Pupper <- (P-1+.5-mu)/sigma
  Pupper <- 1 - (stats::pnorm(Pupper)+gamma/(6*sigma^3)*(1-Pupper^2)*stats::dnorm(Pupper))

  if (signed) {
    Plower <- (P+.5-mu)/sigma
    Plower <- stats::pnorm(Plower)+gamma/(6*sigma^3)*(1-Plower^2)*stats::dnorm(Plower)
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing.as.zero) {
    Pupper[P == 0] <- NA
    if (signed) {Plower[P == 0] <- NA}
  }

  #### Assemble backbone object ####
  bb <- list(G = P,  #Preliminary backbone object
             Pupper = Pupper,
             model = "fixedcol",
             agents = nrow(B),
             artifacts = ncol(B),
             weighted = FALSE,
             bipartite = TRUE,
             symmetric = TRUE,
             class = class,
             trials = NULL)
  if (signed) {bb <- append(bb, list(Plower = Plower))}  #Add lower-tail values, if requested
  if (!is.null(attribs)) {bb <- append(bb, list(attribs = attribs))}  #Add node attributes, if present
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = class, narrative = narrative)
    return(backbone)
  }
}
