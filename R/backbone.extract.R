#' Extracts a backbone network from a backbone object
#'
#' `backbone.extract` returns a binary or signed adjacency matrix
#'      containing the backbone that retains only the significant edges.
#'
#' @param bb.object backbone: backbone S3 class object.
#' @param signed Boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param alpha Real: significance level of hypothesis test(s)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("matrix", "sparseMatrix", "igraph", "edgelist"), converted via \link{tomatrix}.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#' @return backbone graph: Binary or signed backbone graph of class given in parameter `class`.
#'
#' @details The "backbone" S3 class object is composed of (1) the weighted graph as a matrix, (2) upper-tail p-values as a
#'    matrix, (3, if `signed = TRUE`) lower-tail p-values as a matrix, (4, if present) node attributes as a dataframe, and
#'    (5) several properties of the original graph and backbone model
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*
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
#' backbone.object <- fixedrow(B, alpha = NULL)
#' bb <- backbone.extract(backbone.object, alpha = 0.05)
backbone.extract <- function(bb.object, signed = FALSE, alpha = 0.05, mtc = "none", class = bb.object$class, narrative = FALSE){

  #### Argument Checks ####
  if ((alpha >= 1) | (alpha <= 0)) {stop("alpha must be between 0 and 1")}
  if ((class != "matrix")
      & (class != "Matrix")
      & (class != "sparseMatrix")
      & (class != "igraph")
      & (class != "edgelist"))
  {stop("incorrect class type, must be one of c(matrix, Matrix, sparseMatrix, igraph, edgelist)")}
  if (signed == TRUE & is.null(bb.object$Plower)) {stop(paste0("This backbone object does not contain lower-tail p-values, so a signed backbone cannot be extracted.\n       To extract a signed backbone, please re-run ", bb.object$model, "() and specify `signed = TRUE`."))}

  #### Extract object components ####
  G <- bb.object$G
  Pupper <- bb.object$Pupper
  Plower <- bb.object$Plower
  attribs <- bb.object$attribs

  #### Extract signed backbone (two-tailed test; all dyads considered) ####
  if (signed) {
    alpha <- alpha / 2  #Use two-tailed test
    Psmaller <- pmin(Pupper,Plower)  #Find smaller p-value
    diag(Psmaller) <- NA
    Ptail <- (Pupper < Plower)  #Find tail of smaller p-value (TRUE if smaller p-value is in upper tail)
    diag(Ptail) <- NA

    if (mtc != "none") {  #Adjust p-values for familywise error, if requested
      if (isSymmetric(Psmaller)) {Psmaller[upper.tri(Psmaller)] <- NA}  #If undirected, ignore upper triangle
      p <- as.vector(Psmaller)  #Vector of p-values
      m <- sum((!is.na(p))*1)  #Number of p-values to evaluate, number of independent edges to test
      p <- stats::p.adjust(p, method = mtc, m)  #Adjust p-values
      Psmaller <- matrix(p, nrow = nrow(Psmaller), ncol = ncol(Psmaller))  #Put adjusted p-values in original p-value matrix
      if (all(is.na(Psmaller[upper.tri(Psmaller)]))) {Psmaller[upper.tri(Psmaller)] <- t(Psmaller)[upper.tri(Psmaller)]}  #If upper triangle is missing, put it back
    }

    backbone <- (Psmaller < alpha)*1  #Identify all significant edges
    backbone[which(Ptail==FALSE)] <- backbone[which(Ptail==FALSE)] * -1  #Make lower-tail significant edges negative
    backbone[which(is.na(backbone))] <- 0  #fill NAs with 0s
    rownames(backbone) <- rownames(G)
    colnames(backbone) <- rownames(G)
  }

  #### Extract binary backbone (one-tailed test; only positively-weighed edges considered) ####
  if (!signed) {
    Pupper[which(G==0)] <- NA  #Eliminate p-values for zero-weight edges; not relevant
    diag(Pupper) <- NA  #Eliminate p-values for loops; not relevant

    if (mtc != "none") {  #Adjust p-values for familywise error, if requested
      if (isSymmetric(Pupper)) {Pupper[upper.tri(Pupper)] <- NA}  #If undirected, ignore upper triangle
      p <- as.vector(Pupper)  #Vector of p-values
      m <- sum((!is.na(p))*1)  #Number of p-values to evaluate, number of independent edges to test
      p <- stats::p.adjust(p, method = mtc, m)  #Adjust p-values
      Pupper <- matrix(p, nrow = nrow(Pupper), ncol = ncol(Pupper))  #Put adjusted p-values in original p-value matrix
      if (all(is.na(Pupper[upper.tri(Pupper)]))) {Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]}  #If upper triangle is missing, put it back
    }

    backbone <- (Pupper < alpha)*1  #Identify all significant edges
    backbone[which(is.na(backbone))] <- 0  #fill NAs with 0s
    rownames(backbone) <- rownames(G)
    colnames(backbone) <- rownames(G)
  }

  #### Display narrative, if requested ####
  if (narrative) {
    reduced_edges <- round(((sum(G!=0)-nrow(G)) - sum(backbone!=0)) / (sum(G!=0)-nrow(G)),3)*100  #Percent decrease in number of edges
    reduced_nodes <- round((max(sum(rowSums(G)!=0),sum(colSums(G)!=0)) - max(sum(rowSums(backbone)!=0),sum(colSums(backbone)!=0))) / max(sum(rowSums(G)!=0),sum(colSums(G)!=0)),3) * 100  #Percent decrease in number of connected nodes
    write.narrative(agents = bb.object$agents,
                    artifacts = bb.object$artifacts,
                    weighted = bb.object$weighted,
                    bipartite = bb.object$bipartite,
                    symmetric = bb.object$symmetric,
                    signed = signed,
                    mtc = mtc,
                    alpha = alpha,
                    s = NULL,
                    ut = NULL,
                    lt = NULL,
                    trials = NULL,
                    model = bb.object$model,
                    reduced_edges = reduced_edges,
                    reduced_nodes = reduced_nodes)
  }

  #### Return result ####
  backbone <- frommatrix(backbone, attribs, convert = class)
  return(backbone)
}
