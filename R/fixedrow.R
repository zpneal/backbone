#' Extract backbone using the Fixed Row Model
#'
#' `fixedrow` extracts the backbone of a bipartite projection using the Fixed Row Model.
#'
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix, sparse \code{\link{Matrix}}, or dataframe; (2) an edgelist in the form of a two-column matrix, sparse \code{\link{Matrix}}, or dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param fwer string: type of familywise error rate correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative boolean: TRUE if suggested text & citations should be displayed.
#'
#' @details
#' The `fixedrow` function compares an edge's observed weight in the projection \eqn{B*t(B)} to the
#'     distribution of weights expected in a projection obtained from a random bipartite graph where
#'     the *row* vertex degrees are fixed but the column vertex degrees are allowed to vary.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected under the null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*

#' @return
#' If `alpha` != NULL: Binary or signed backbone graph of class `class`.
#'
#' If `alpha` == NULL: An S3 backbone object containing three matrices (the weighted graph, edges' upper-tail p-values,
#'    edges' lower-tail p-values), and a string indicating the null model used to compute p-values, from which a backbone
#'    can subsequently be extracted using [backbone.extract()]. The `signed`, `fwer`, `class`, and `narrative` parameters
#'    are ignored.
#'
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
#' bb <- fixedrow(B, alpha = 0.05, narrative = TRUE, class = "igraph") #A fixedrow backbone...
#' plot(bb) #...is sparse with clear communities

fixedrow <- function(B, alpha = NULL, signed = FALSE, fwer = "none", class = "original", narrative = FALSE){

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
  rs <- rowSums(B)

  #### Hypergeometric Distribution ####
  ### Set up df for values ###
  df <- data.frame(as.vector(P))
  names(df)[names(df)=="as.vector.P."] <- "projvalue"

  ### Compute row sums ###
  df$row_sum_i <- rep(rs, times = nrow(B))

  ### Match each row sum i with each row sum j and their Pij value ###
  df$row_sum_j <- rep(rs, each = nrow(B))

  ### Compute different in number of artifacts and row sum ###
  df$diff <- ncol(B)-df$row_sum_i

  ### Probability of Pij or less ###
  df$hgl <- stats::phyper(df$projvalue, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)

  ### Probability of Pij or more ###
  df$hgu <- stats::phyper(df$projvalue-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)

  #### Create Backbone Object ####
  Pupper <- matrix(as.numeric(df$hgu), nrow = nrow(B), ncol = nrow(B))
  Plower <- matrix(as.numeric(df$hgl), nrow = nrow(B), ncol = nrow(B))
  bb <- list(G = P, Pupper = Pupper, Plower = Plower, model = "fixedrow")
  class(bb) <- "backbone"

  #### Return result ####
  if (is.null(alpha)) {return(bb)}  #Return backbone object if `alpha` is not specified
  if (!is.null(alpha)) {            #Otherwise, return extracted backbone (and show narrative text if requested)
    backbone <- backbone.extract(bb, alpha = alpha, signed = signed, fwer = fwer, class = "matrix")
    retained <- round((sum((backbone!=0)*1)) / (sum((P!=0)*1) - nrow(P)),3)*100
    if (narrative == TRUE) {write.narrative(agents = nrow(B), artifacts = ncol(B), weighted = FALSE, bipartite = TRUE, symmetric = TRUE,
                                            signed = signed, fwer = fwer, alpha = alpha, s = NULL, ut = NULL, lt = NULL, trials = NULL, model = "fixedrow", retained = retained)}
    backbone <- frommatrix(backbone, convert = class)
    return(backbone)
  }
}

#' Wrapper for fixedrow()
#' @param B An unweighted bipartite graph, as: (1) an incidence matrix in the form of a matrix, sparse \code{\link{Matrix}}, or dataframe; (2) an edgelist in the form of a two-column matrix, sparse \code{\link{Matrix}}, or dataframe; (3) an \code{\link{igraph}} object; (4) a \code{\link{network}} object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#' @param alpha Real: significance level of hypothesis test(s)
#' @param signed Boolean: TRUE if signed backbone is to be returned, FALSE if binary backbone is to be returned
#' @param fwer string: type of familywise error rate correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("original", "matrix", "sparseMatrix", "igraph", "network", "edgelist").
#'     If "original", the backbone graph returned is of the same class as `B`.
#' @param narrative Boolean: TRUE if suggested text for a manuscript is to be returned.
#' @export
hyperg <- function(B, alpha = NULL, signed = FALSE, fwer = "none", class = "original", narrative = FALSE){
  warning("The hyperg() function is now called fixedrow(); please use fixedrow() instead.")
  fixedrow(B, alpha = alpha, signed = signed, fwer = fwer, class = class, narrative = narrative)
}
