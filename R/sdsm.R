#' Extract backbone using the Stochastic Degree Sequence Model
#'
#' `sdsm` extracts the backbone of a bipartite projection using the Stochastic Degree Sequence Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.
#' @param alpha real: significance level of hypothesis test(s)
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
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
#' When `signed = FALSE`, a one-tailed test (is the weight stronger?) is performed for each edge. The resulting backbone
#'    contains edges whose weights are significantly *stronger* than expected in the null model. When `signed = TRUE`, a
#'    two-tailed test (is the weight stronger or weaker?) is performed for each edge. The resulting backbone contains
#'    positive edges for those whose weights are significantly *stronger*, and negative edges for those whose weights are
#'    significantly *weaker*, than expected in the null model.
#'
#' The bipartite network `B` may contain some edges that are *required* in the null model (i.e., structural 1s); these edges should
#'    have a weight of 11 (i.e., B_ik = 11). This network may also contain some edges that are *prohibited* in the null model
#'    (i.e., structural 0s); these edges should have a weight of 10 (i.e., B_ik = 10). When `B` contains required or prohibited edges,
#'    cellwise probabilities are computed using \link{logit} following Neal et al. (2024). Otherwise, cellwise probabilities are
#'    computed using the faster and more accurate Bipartite Configuration Model with \link{bicm} (Neal et al. 2021).
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
#' @references sdsm: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @references bicm: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#' @references logit: {Neal, Z. P. and Neal, J. W. (2024). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. *Proceedings of the 12th International Conference on Complex Networks and their Applications*. Springer.}
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

sdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", class = "original", narrative = FALSE, ...){

  #### Argument Checks ####
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Class Conversion ####
  convert <- tomatrix(B)
  if (class == "original") {class <- convert$summary$class}
  attribs <- convert$attribs
  B <- convert$G
  if (convert$summary$weighted==TRUE) {if (any(!(B%in%c(0,1,10,11)))) {stop("Graph must be unweighted.")}}  #If G is not binary, check whether it only contains structural 0s or 1s
  if (convert$summary$bipartite==FALSE){
    warning("This object is being treated as a bipartite network.")
    convert$summary$bipartite <- TRUE
    }

  #### Bipartite Projection ####
  B_unweight <- B
  B_unweight[B_unweight==10] <- 0  #Make structural 0s ordinary 0
  B_unweight[B_unweight==11] <- 1  #Make structural 1s ordinary 1
  P <- tcrossprod(B_unweight)  #Projection, not considering any structural 0s or 1s

  #### Compute Probabilities for SDSM ####
  if (any(B%in%c(10,11))) {probs <- logit(B)} else {probs <- bicm(B, ...)}
  probs <- lapply(seq_len(nrow(probs)), function(i) probs[i,])  #Store probabilities as list

  #### Compute p-values (for unsigned backbone, ignore lower-tail p-values) ####
  if (!signed) {
    Pupper <- matrix(NA, nrow(P), ncol(P))  #Set upper-tail p-value to NA initially
    for (col in 1:ncol(P)) {  #Loop over lower triangle
      for (row in col:nrow(P)) {

        if (missing.as.zero) {  #If missing edges should be treated as zero, test each one
          pvalues <- pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)
          Pupper[row,col] <- pvalues[2]
        }

        if (!missing.as.zero & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
          pvalues <- pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)
          Pupper[row,col] <- pvalues[2]
        }

      }
    }
    Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Add upper triangle
  }

  #### Compute p-values (for signed backbone) ####
  if (signed) {
    Pupper <- matrix(NA, nrow(P), ncol(P))
    Plower <- matrix(NA, nrow(P), ncol(P))
    for (col in 1:ncol(P)) {  #Loop over lower triangle
      for (row in col:nrow(P)) {

        if (missing.as.zero) {  #If missing edges should be treated as zero, test each one
          pvalues <- pb(P[row,col], unlist(Map('*',probs[row],probs[col])))
          Plower[row,col] <- pvalues[1]
          Pupper[row,col] <- pvalues[2]
        }

        if (!missing.as.zero  & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
          pvalues <- pb(P[row,col], unlist(Map('*',probs[row],probs[col])))
          Plower[row,col] <- pvalues[1]
          Pupper[row,col] <- pvalues[2]
        }

      }
    }
    Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Add upper triangles
    Plower[upper.tri(Plower)] <- t(Plower)[upper.tri(Plower)]
  }

  #### Assemble backbone object ####
  bb <- list(G = P,  #Preliminary backbone object
             Pupper = Pupper,
             model = "sdsm",
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
