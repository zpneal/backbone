#' Extract backbone using the Marginal Likelihood Filter
#'
#' `mlf` extracts the backbone of a weighted network using the Marginal Likelihood Filter
#'
#' @param W An integer-weighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `W`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `mlf` function applies the marginal likelihood filter (MLF; Dianati, 2016), which compares an edge's weight to
#'    its expected weight in a graph that preserves the total weight and preserves the degree sequence *on average*.
#'    The graph may be directed or undirected, however the edge weights must be positive integers.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'
#' If `W` is an unweighted bipartite graph, then the MLF is applied to its weighted bipartite projection.
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
#' @references mlf: {Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. *Physical Review E, 93*, 012304. \doi{10.1103/PhysRevE.93.012304}}
#' @export
#'
#' @examples
#' #A network with heterogeneous weights
#' net <- matrix(c(0,10,10,10,10,75,0,0,0,0,
#'                 10,0,1,1,1,0,0,0,0,0,
#'                 10,1,0,1,1,0,0,0,0,0,
#'                 10,1,1,0,1,0,0,0,0,0,
#'                 10,1,1,1,0,0,0,0,0,0,
#'                 75,0,0,0,0,0,100,100,100,100,
#'                 0,0,0,0,0,100,0,10,10,10,
#'                 0,0,0,0,0,100,10,0,10,10,
#'                 0,0,0,0,0,100,10,10,0,10,
#'                 0,0,0,0,0,100,10,10,10,0),10)
#'
#' net <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE)
#' plot(net, edge.width = sqrt(igraph::E(net)$weight)) #A stronger clique & a weaker clique
#'
#' strong <- igraph::delete.edges(net, which(igraph::E(net)$weight < mean(igraph::E(net)$weight)))
#' plot(strong) #A backbone of stronger-than-average edges ignores the weaker clique
#'
#' bb <- mlf(net, alpha = 0.05, narrative = TRUE) #An MLF backbone...
#' plot(bb) #...preserves edges at multiple scales
mlf <- function(W, alpha = 0.05, signed = FALSE, mtc = "none", class = "original", narrative = FALSE){

  #### Argument Checks ####
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Class Conversion ####
  convert <- tomatrix(W)
  G <- convert$G
  if (class == "original") {class <- convert$summary$class}
  attribs <- convert$attribs
  symmetric <- convert$summary$symmetric
  if (convert$summary$bipartite==TRUE){
    message("The input graph is bipartite; extraction is performed on its unipartite projection.")
    G <- tcrossprod(G)
    }

  # Check for possible bipartite projection
  if (all(G%%1==0) &                               #If all entries are integers, and
      any(!(diag(G)%in%c(0,1,NA))) &               #The diagonal is present, and not only 0s and 1s, and
      all((diag(G) == apply(G, 1, FUN=max)))) {    #The diagonal is the largest entry in each row
    message("This object looks like it could be a bipartite projection. If so, consider extracting the backbone using a model designed for bipartite projections: sdsm, fdsm, fixedfill, fixedrow, or fixedcol.")
  }

  #### Compute p-values ####
  if (symmetric) {
    Pupper <- matrix(1, nrow(G), ncol(G))  #Set p-values to 1 initially
    if (signed) {Plower <- matrix(1, nrow(G), ncol(G))}
    T <- sum(rowSums(G))/2
    p <- (rowSums(G) %*% t(rowSums(G))) / (2 * (T^2))
    for (col in 1:ncol(G)) {  #Loop over lower triangle
      for (row in col:nrow(G)) {
        if (G[row,col] != 0) {  #Compute and update p-values only if the edge has non-zero weight
          Pupper[row,col] <- stats::binom.test(G[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {Plower[row,col] <- stats::binom.test(G[row,col], T, p[row,col], alternative = "less")$p.value}
        }
      }
    }
    Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Add upper triangle
    if (signed) {Plower[upper.tri(Plower)] <- t(Plower)[upper.tri(Plower)]}
  }

  if (!symmetric) {
    Pupper <- matrix(1, nrow(G), ncol(G))  #Set p-values to 1 initially
    if (signed) {Plower <- matrix(1, nrow(G), ncol(G))}
    T <- sum(rowSums(G))
    p <- (rowSums(G) %*% t(colSums(G))) / (T^2)
    for (col in 1:ncol(G)) {  #Loop over full matrix
      for (row in 1:nrow(G)) {
        if (G[row,col] != 0) {  #Compute and update p-values only if the edge has non-zero weight
          Pupper[row,col] <- stats::binom.test(G[row,col], T, p[row,col], alternative = "greater")$p.value
          if (signed) {Plower[row,col] <- stats::binom.test(G[row,col], T, p[row,col], alternative = "less")$p.value}
        }
      }
    }
  }

  ### Create backbone object ###
  bb <- list(G = G,  #Preliminary backbone object
             Pupper = Pupper,
             model = "mlf",
             agents = nrow(G),
             artifacts = NULL,
             weighted = TRUE,
             bipartite = FALSE,
             symmetric = symmetric,
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
