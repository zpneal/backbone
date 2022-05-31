#' Extract backbone using the Fixed Column Model
#'
#' `fixedcol` extracts the backbone of a bipartite projection using the Fixed Column Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param alpha real: significance level of hypothesis test(s)
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
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
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
#' bb <- fixedcol(B, alpha = 0.05, narrative = TRUE, class = "igraph") #A fixedcol backbone...
#' plot(bb) #...is sparse with clear communities

fixedcol <- function(B, alpha = 0.05, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

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
  cs <- colSums(B)
  m = dim(B)[1]

  #### Create Backbone Object ####
  pt <- ((cs*(cs-1))/(m*(m-1)))  #Parameters for computing Poissonbinomial
  Pupper <- matrix((pb(k = (P-1), p = pt, lower = FALSE)),nrow = m, ncol = m,byrow = TRUE)
  Plower <- matrix(pb(k = P, p = pt),nrow = m, ncol = m,byrow = TRUE)
  bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "fixedcol")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = "matrix")
    reduced_edges <- round(((sum(P!=0)-nrow(P)) - sum(backbone!=0)) / (sum(P!=0)-nrow(P)),3)*100  #Percent decrease in number of edges
    reduced_nodes <- round((max(sum(rowSums(P)!=0),sum(colSums(P)!=0)) - max(sum(rowSums(backbone)!=0),sum(colSums(backbone)!=0))) / max(sum(rowSums(P)!=0),sum(colSums(P)!=0)),3) * 100  #Percent decrease in number of connected nodes
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = FALSE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, mtc = mtc, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "fixedcol",
                                            reduced_edges = reduced_edges, reduced_nodes = reduced_nodes)}
    backbone <- frommatrix(backbone, attribs, convert = class)
    return(backbone)
  }
}
