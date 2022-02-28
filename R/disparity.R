#' Extract backbone using the Disparity Filter
#'
#' `disparity` extracts the backbone of a weighted network using the Disparity Filter.
#'
#' @param W A weighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a three-column dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param fwer string: type of familywise error rate correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `W`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `disparity` function applies the disparity filter (Serrano et al., 2009), which compares an edge's weight to
#'    its expected weight if a node's total degree was uniformly distributed across all its edges. The graph may be
#'    directed or undirected, however the edge weights must be positive.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*
#'
#' If `W` is an unweighted bipartite graph, any rows and columns that contain only zeros or only ones are removed, then
#'    the global threshold is applied to its weighted bipartite projection.
#'
#' @return
#' If `alpha` != NULL: Binary or signed backbone graph of class `class`.
#'
#' If `alpha` == NULL: An S3 backbone object containing three matrices (the weighted graph, edges' upper-tail p-values,
#'    edges' lower-tail p-values), and a string indicating the null model used to compute p-values, from which a backbone
#'    can subsequently be extracted using [backbone.extract()]. The `signed`, `fwer`, `class`, and `narrative` parameters
#'    are ignored.
#'
#' @references {Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. *Proceedings of the National Academy of Sciences, 106*, 6483-6488. \doi{10.1073/pnas.0808904106}}
#' @export
#'
#' @examples
#' #A network with heterogeneous (i.e. multiscale) weights
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
#' bb <- disparity(net, alpha = 0.05, narrative = TRUE) #A disparity backbone...
#' plot(bb) #...preserves edges at multiple scales
disparity <- function(W, alpha = 0.05, signed = FALSE, fwer = "none", class = "original", narrative = FALSE){

  #### Argument Checks ####
  if (!is.null(alpha)) {if (alpha < 0 | alpha > .5) {stop("alpha must be between 0 and 0.5")}}

  #### Class Conversion ####
  convert <- tomatrix(W)
  G <- convert$G
  if (class == "original") {class <- convert$summary$class}
  symmetric <- convert$summary$symmetric
  if (convert$summary$bipartite==TRUE){
    message("The input graph is bipartite; extraction is performed on its unipartite projection.")
    artifacts <- ncol(G)
    G <- tcrossprod(G)
    }

  # Check for possible bipartite projection
  if (all(G%%1==0) &                               #If all entries are integers, and
      any(!(diag(G)%in%c(0,1,NA))) &               #The diagonal is present, and not only 0s and 1s, and
      all((diag(G) == apply(G, 1, FUN=max)))) {    #The diagonal is the largest entry in each row
    message("This object looks like it could be a bipartite projection. If so, consider extracting the backbone using a model designed for bipartite projections: sdsm, fdsm, fixedfill, fixedrow, or fixedcol.")
  }

  #### Set Parameters and Compute p-values ####
  strength <- rowSums(G)
  binary <- (G>0)+0
  degree <- rowSums(binary)
  zeros <- G==0

  if (symmetric == TRUE){
    P <- G/strength
    pvalues <- (1-P)^(degree-1)
    Pupper <- as.matrix(pvalues)          #Asymmetric p-values, one from the perspective of each node
    Pupper <- pmin(Pupper,t(Pupper))  #From Serrano: "satisfy the above criterion for at least one of the two nodes"
    Plower <- 1-Pupper
  }

  if (symmetric == FALSE){
    ### Implies Directed ###
    outp <- G/strength
    outvalues <- (1-outp)^(degree-1)
    inp <- t(G)/(colSums(G))
    invalues <- t((1-inp)^(colSums(binary)-1))
    Pupper <- pmin(invalues,outvalues)
    Plower <- 1-Pupper
  }

  ### If edge weight was zero, set to 1 in positive and negative so edge is not in backbone ###
  Pupper[zeros] <- 1
  Plower[zeros] <- 1

  ### Create backbone object ###
  bb <- list(G = G, Pupper = Pupper, Plower = Plower, model = "disparity")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, fwer = fwer, class = "matrix")
    retained <- round((sum((backbone!=0)*1)) / (sum((G!=0)*1) - nrow(G)),3)*100
    if (narrative == TRUE) {write.narrative(agents = nrow(G), artifacts = NULL, weighted = TRUE, bipartite = FALSE, symmetric = symmetric,
                                            signed = signed, fwer = fwer, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "disparity", retained = retained)}
    backbone <- frommatrix(backbone, convert = class)
    return(backbone)
  }
}
