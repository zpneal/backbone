#' Extract backbone using the Stochastic Degree Sequence Model
#'
#' `sdsm` extracts the backbone of a bipartite projection using the Stochastic Degree Sequence Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' @param ... optional arguments
#'
#' @details
#' The `sdsm` function compares an edge's observed weight in the projection \code{B*t(B)} to the distribution of weights
#'    expected in a projection obtained from a random bipartite network where both the row vertex degrees and column
#'    vertex degrees are *approximately* fixed at their values in `B`. It uses the Bipartite Configuration Model \link{bicm}
#'    to compute probabilities for the Poisson binomial distribution.
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
#' If `alpha` == NULL: An S3 backbone object containing (1) the weighted graph as a matrix, (2) upper-tail p-values as a
#'    matrix, (3, if `signed = TRUE`) lower-tail p-values as a matrix, (4, if present) node attributes as a dataframe, and
#'    (5) a string indicating the null model used to compute p-values, from which a backbone can subsequently be extracted
#'    using [backbone.extract()]. The `mtc`, `class`, and `narrative` parameters are ignored.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references sdsm: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @references sdsm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
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
#' bb <- sdsm(B, alpha = 0.05, narrative = TRUE, class = "igraph") #An SDSM backbone...
#' plot(bb) #...is sparse with clear communities

sdsm <- function(B, alpha = 0.05, signed = FALSE, mtc = "none", class = "original", narrative = FALSE, ...){

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

  #### Compute Probabilities for SDSM ####
  prob.mat <- bicm(B,...)

  #### Compute p-values (for unsigned backbone, ignore lower-tail p-values) ####
  if (!signed) {
    Pupper <- matrix(1, nrow(P), ncol(P))  #Set upper-tail p-value to 1 initially
    for (col in 1:ncol(P)) {  #Loop over lower triangle
      for (row in col:nrow(P)) {
        if (P[row,col] != 0) {  #Compute and update the upper-tail p-value only if the edge has non-zero weight
          pvalues <- pb(P[row,col], prob.mat[row,]*prob.mat[col,], lowertail = FALSE)
          Pupper[row,col] <- pvalues[2]
        }
      }
    }
    Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Add upper triangle
  }

  #### Compute p-values (for signed backbone) ####
  if (signed) {
    Pupper <- matrix(0, nrow(P), ncol(P))
    Plower <- matrix(0, nrow(P), ncol(P))
    for (col in 1:ncol(P)) {  #Loop over lower triangle
      for (row in col:nrow(P)) {
        pvalues <- pb(P[row,col], prob.mat[row,]*prob.mat[col,])
        Plower[row,col] <- pvalues[1]
        Pupper[row,col] <- pvalues[2]
      }
    }
    Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Add upper triangles
    Plower[upper.tri(Plower)] <- t(Plower)[upper.tri(Plower)]
  }

  #### Assemble backbone object ####
  bb <- list(G = P, Pupper = Pupper, model = "sdsm")  #Preliminary backbone object
  if (signed) {bb <- append(bb, list(Plower = Plower))}  #Add lower-tail values, if requested
  if (!is.null(attribs)) {bb <- append(bb, list(attribs = attribs))}  #Add node attributes, if present
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, mtc = mtc, class = "matrix")
    reduced_edges <- round(((sum(P!=0)-nrow(P)) - sum(backbone!=0)) / (sum(P!=0)-nrow(P)),3)*100  #Percent decrease in number of edges
    reduced_nodes <- round((max(sum(rowSums(P)!=0),sum(colSums(P)!=0)) - max(sum(rowSums(backbone)!=0),sum(colSums(backbone)!=0))) / max(sum(rowSums(P)!=0),sum(colSums(P)!=0)),3) * 100  #Percent decrease in number of connected nodes
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = FALSE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, mtc = mtc, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "sdsm",
                                            reduced_edges = reduced_edges, reduced_nodes = reduced_nodes)}
    backbone <- frommatrix(backbone, attribs, convert = class)
    return(backbone)
  }
}
