#' Extract the backbone from a bipartite projection
#'
#' \code{backbone_from_bipartite} extracts the unweighted backbone from the weighted projection of a bipartite network.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix, a binary incidence \code{\link{Matrix}}, or a binary bipartite \code{\link{igraph}} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param model string: backbone model. This must be onf of: \code{"sdsm"}, \code{"fdsm"}, \code{"fixedrow"}, \code{"fixedcol"}, or \code{"fixedfill"}
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction; can be either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param missing_as_zero boolean: treat missing edges be treated as edges with zero weight and test them for significance
#' @param only_pvalues boolean: return only a matrix of edgewise p-values; all parameters except \code{model} are ignored
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' The \code{backbone_from_bipartite} extracts the backbone from the weighted projection of a bipartite network. The backbone is
#' is unweighted network that contains only edges whose weights in the projection are statistically significant. When \code{signed = FALSE},
#' the backbone contains edges that are statistically significantly strong under a one-tailed test. When \code{signed = TRUE},
#' the backbone contains positive edges that are statistically significantly strong, and negative edges that are statistically
#' significantly weak, under a two-tailed test.
#'
#' The \code{model} parameter controls the null model used to evaluate the statistical significance of edge weights. Each null model
#' imposes a unique set of constraints on \code{B} (Neal et al., 2021). In increasing order of constraints and computational complexity:
#' * \code{fixedfill} - Exactly constrain the total number of edges in (i.e., sum of) \code{B}
#' * \code{fixedrow} - Exactly constrain the agent degrees in (i.e., row sums of) \code{B}
#' * \code{fixedcol} - Exactly constrain the artifact degrees in (i.e., column sums of) \code{B}
#' * \code{sdsm} - Approximately constrain the agent and artifact degrees in \code{B} (the default)
#' * \code{fdsm} - Exactly constrain the agent and artifact degrees in \code{B}
#'
#' @return A backbone in the same class as \code{B} (or if \code{only_pvalues = TRUE}, a Matrix of edgewise p-values)
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references models: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @export
#'
#' @examples
#' #TBD
backbone_from_bipartite <- function(B,
                                    alpha = 0.05,
                                    model = "sdsm",
                                    signed = FALSE,
                                    mtc = "none",
                                    missing_as_zero = FALSE,
                                    only_pvalues = FALSE,
                                    narrative = FALSE) {

  #### Check parameters ####
  if (!is.numeric(alpha)) {stop("`alpha` must be a numeric value between 0 and 1")}
  if (alpha < 0 | alpha > 1) {stop("`alpha` must be a numeric value between 0 and 1")}
  if (!(model %in% c("sdsm", "fdsm", "fixedrow", "fixedcol", "fixedfill"))) {stop("`model` must be one of: \"sdsm\", \"fdsm\", \"fixedrow\", \"fixedcol\", or \"fixedfill\"")}
  if (!is.logical(signed)) {stop("`signed` must be either TRUE or FALSE")}
  if (!(mtc %in% c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {stop("`mtc` must be one of: \"none\", \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", or \"fdr\"")}
  if (!is.logical(missing_as_zero)) {stop("`missing_as_zero` must be either TRUE or FALSE")}
  if (!is.logical(only_pvalues)) {stop("`only_pvalues` must be either TRUE or FALSE")}
  if (!is.logical(narrative)) {stop("`narrative` must be either TRUE or FALSE")}

  #### Check and format input ####
  #Check that input is matrix, Matrix, or igraph (and if igraph, that it is bipartite)
  if (!methods::is(B,"matrix") & !methods::is(B,"Matrix") & !methods::is(B,"igraph")) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}
  if (methods::is(B,"igraph")) {if(!igraph::is_bipartite(B)) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}}

  #Convert input to sparse matrix
  if (methods::is(B,"matrix")) {mat <- Matrix::Matrix(B)}
  if (methods::is(B,"Matrix")) {mat <- B}
  if (methods::is(B,"igraph")) {mat <- igraph::as_biadjacency_matrix(B, names = FALSE, sparse = TRUE)}

  #Check if input may be a weighted projection
  if (!all(as.vector(mat) %in% c(0,1)) &    #The entries are not binary, and
      Matrix::isSymmetric(mat) &            #The matrix is symmetric, and
      all(as.vector(mat)%%1==0)) {          #The entries are all integers
      stop("`B` looks like it may be a weighted bipartite projection. The input to backbone_from_bipartite()
       must be the original bipartite network, not its weighted projection. If you only have the weighted
       bipartite projection, cautiously consider using backbone_from_weighted() instead.")}

  #Check that input is binary
  if (!all(as.vector(mat) %in% c(0,1))) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}

  #### Compute p-values ####


  #### Retain edges ####


  #### Display narrative ####


  #### Return backbone ####

}
