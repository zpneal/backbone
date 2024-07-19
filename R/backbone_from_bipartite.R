#' Extract the backbone from a bipartite projection
#'
#' \code{backbone_from_bipartite} extracts the unweighted backbone from the weighted projection of a bipartite network.
#'
#' @param B An unweighted bipartite network as a binary incidence matrix or a binary bipartite \code{\link{igraph}} object (see details)
#' @param alpha real: significance level of hypothesis test(s)
#' @param model string: backbone model. This must be one of: \code{"sdsm"}, \code{"sdsm-ec"} \code{"fdsm"}, \code{"fixedrow"}, \code{"fixedcol"}, or \code{"fixedfill"}
#' @param signed boolean: return a signed backbone
#' @param mtc string: type of Multiple Test Correction; can be either \code{"none"} or a method allowed by \code{\link{p.adjust}}.
#' @param missing_as_zero boolean: treat missing edges be treated as edges with zero weight and test them for significance
#' @param only_pvalues boolean: return only a matrix of edgewise p-values; all parameters except \code{model} are ignored
#' @param trials integer: number of Monte Carlo trials used for FDSM (ignored if \code{model != "fdsm"})
#' @param narrative boolean: display suggested text & citations
#'
#' @details
#' The \code{backbone_from_bipartite} extracts the backbone from the weighted projection of a bipartite network composed of *n* "agent"
#' nodes and *m* "artifact" nodes. The backbone is an unweighted unipartite network of agents that contains only edges whose weights
#' in the projection are statistically significant. When \code{signed = FALSE}, the backbone contains edges that are statistically
#' significantly strong under a one-tailed test. When \code{signed = TRUE}, the backbone contains positive edges that are statistically
#' significantly strong, and negative edges that are statistically significantly weak, under a two-tailed test.
#'
#' The \code{model} parameter controls the null model used to evaluate the statistical significance of edge weights, each of which
#' imposes a unique set of constraints on \code{B}:
#' * \code{fixedfill} - Use the "fixed fill" model (Neal et al., 2021), which exactly constrains the total number of edges (i.e., sum)
#' * \code{fixedrow} - Use the "fixed row" model (Neal et al., 2021), which exactly constrains the agent degrees (i.e., row sums)
#' * \code{fixedcol} - Use the "fixed column" model (Neal et al., 2021), which exactly constrains the artifact degrees (i.e., column sums)
#' * \code{sdsm} - Use the "Stochastic Degree Sequence Model" (SDSM; Neal et al., 2021), which pproximately constrains the agent and artifact degrees (the default)
#' * \code{sdsm-ec} - Use the "SDSM with Edge Constraints" (Neal & Neal, 2023), which approximately constrains the agent and artifact degrees, and exactly constrains edges that are prohibited (weight = 10) or required (weight = 11)
#' * \code{fdsm} - Use the "Fixed Degree Sequence Model" (Neal et al., 2021), which exactly constrain the agent and artifact degrees
#'
#' Although \cite{backbone_from_bipartite} extracts the backbone from a weighted bipartite projection, the input \code{B} must be the
#' bipartite network itself, and not the weighted projection. This is necessary because the backbone models use information in the bipartite
#' network that is missing from the projection. The "agent" nodes that appear in the projection must be represented by rows if \code{B}
#' is an incidence matrix, or \code{type = FALSE} nodes if \code{B} is a bipartite igraph object. In either case, the bipartite network
#' must be binary (i.e., unweighted), unless \code{model = "sdsm-ec"}, when prohibited" edges can be represented with weight = 10
#' and "required" edges can be represented with weight = 11.
#'
#' @return A backbone in the same class as \code{B} (or if \code{only_pvalues = TRUE}, a matrix of edgewise p-values)
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references sdsm-ec model: {Neal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. *International Conference on Complex Networks and Their Applications, 12*, 127-136. \doi{10.1007/978-3-031-53468-3_11}}
#' @references all other models: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
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
                                    trials = NULL,
                                    narrative = FALSE) {

  #### Check parameters ####
  if (!is.numeric(alpha)) {stop("`alpha` must be a numeric value between 0 and 1")}
  if (alpha < 0 | alpha > 1) {stop("`alpha` must be a numeric value between 0 and 1")}
  if (!(model %in% c("sdsm", "fdsm", "fixedrow", "fixedcol", "fixedfill"))) {stop("`model` must be one of: \"sdsm\", \"sdsm-ec\", \"fdsm\", \"fixedrow\", \"fixedcol\", or \"fixedfill\"")}
  if (!is.logical(signed)) {stop("`signed` must be either TRUE or FALSE")}
  if (!(mtc %in% c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {stop("`mtc` must be one of: \"none\", \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", or \"fdr\"")}
  if (!is.logical(missing_as_zero)) {stop("`missing_as_zero` must be either TRUE or FALSE")}
  if (!is.logical(only_pvalues)) {stop("`only_pvalues` must be either TRUE or FALSE")}
  if (!is.logical(narrative)) {stop("`narrative` must be either TRUE or FALSE")}

  #### Check and format input ####
  #Check that input is matrix or igraph (and if igraph, that it is bipartite)
  if (!methods::is(B,"matrix") & !methods::is(B,"igraph")) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}
  if (methods::is(B,"igraph")) {if(!igraph::is_bipartite(B)) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}}

  #Convert input to matrix
  if (methods::is(B,"matrix")) {I <- B}  #matrix --> matrix
  if (methods::is(B,"igraph")) {
    if ("weight" %in% igraph::edge_attr_names(B)) {I <- igraph::as_biadjacency_matrix(B, names = FALSE, sparse = TRUE, attr = "weight")}  #weighted igraph --> weighted incidence
    if (!("weight" %in% igraph::edge_attr_names(B))) {I <- igraph::as_biadjacency_matrix(B, names = FALSE, sparse = TRUE)}  #unweighted igraph --> binary incidence
  }

  #Check if input may be a weighted projection
  if (!all(I %in% c(0,1)) &    #The entries are not binary, and
      isSymmetric(I) &         #The matrix is symmetric, and
      all(I%%1==0)) {          #The entries are all integers
      stop("`B` looks like it may be a weighted bipartite projection. The input to backbone_from_bipartite()
       must be the original bipartite network, not its weighted projection. If you only have the weighted
       bipartite projection, cautiously consider using backbone_from_weighted() instead.")}

  #Check that input is binary, or contains structural values and model=SDSM
  if (model!="sdsm-ec" & !all(I %in% c(0,1))) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object")}

  if (model=="sdsm-ec" & !all(I %in% c(0,1,10,11))) {stop("`B` must be a binary incidence matrix or binary bipartite igraph object,
                                                          where required edges have weight 10 and prohibited edges have weight 11")}

  #### Compute p-values ####
  if (model == "sdsm") {p <- .sdsm(I, missing_as_zero, signed)}

  #### Retain edges ####
  backbone <- .retain(p, signed, alpha, mtc)

  #### Display narrative ####
  if (narrative) {
  # First sentence (descriptive)
  if (signed) {signed <- "signed"} else {signed <- "unweighted"}

  text <- paste0("We used the backbone package for R (v", utils::packageVersion("backbone"), "; Neal, 2022) to extract the ", signed, " backbone of the weighted projection of an unweighted bipartite network containing ", nrow(I), " agents and ", ncol(I), "artifacts.")

  # Second sentence (model)
  if (mtc == "bonferroni") {correction <- ", Bonferroni adjusted"}
  if (mtc == "holm") {correction <- ", Holm adjusted"}
  if (mtc == "hommel") {correction <- ", Hommel adjusted"}
  if (mtc == "hochberg") {correction <- ", Hochberg adjusted"}
  if (mtc == "BH" | mtc == "fdr") {correction <- ", Benjamini & Hochberg adjusted"}
  if (mtc == "BY") {correction <- ", Benjamini & Yekutieli adjusted"}

  if (model == "fixedfill") {desc <- "the fixed fill model (FFM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "fixedrow") {desc <- "the fixed row model (FRM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "fixedcol") {desc <- "the fixed column model (FCM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "sdsm") {desc <- "the stochastic degree sequence model (SDSM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "sdsm-ec") {desc <- "the stochastic degree sequence model with edge constraints (SDSM-EC; Neal & Neal, 2023)"}
  if (model == "fdsm") {desc <- paste0("the fixed degree sequence model (FDSM; Neal, Domagalski, and Sagan, 2021), where p-values were estimated from ", trials, " Monte Carlo trials")}

  text <- paste0(text, " An edge was retained in the backbone if its weight was statistically significant (alpha = ", alpha, correction, ") using ", desc, ".")

  # Third sentence (reduction)
  old <- sum(p$upper!=0, na.rm=TRUE)  #Number of edges in projection (i.e., number of edges tested, and that have an upper-tail p-value)
  new <- sum(backbone!=0)  #Number of edges in backbone
  reduced_edges <- round(((old - new) / old)*100,2)

  text <- paste0(text, " This reduced the number of edges by ", reduced_edges, "%.")

  # Display
  message("")
  message("=== Suggested text and citations ===")
  message(text)
  message("")
  message("Neal, Z. P. 2022. backbone: An R Package to Extract Network Backbones. PLOS ONE, 17, e0269137. https://doi.org/10.1371/journal.pone.0269137")
  message("")
  if (model %in% c("sdsm", "fdsm", "fixedrow", "fixedcol", "fixedfill")) {message("Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "sdsm-ec") {message("Neal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. International Conference on Complex Networks and Their Applications, 12, 127-136. https://doi.org/10.1007/978-3-031-53468-3_11")}
  }

  #### Return backbone ####
  if (methods::is(B,"matrix")) {
    rownames(backbone) <- rownames(B)
    colnames(backbone) <- rownames(B)
    return(backbone)
  }

  if (methods::is(B,"igraph")) {
    P <- igraph::bipartite_projection(B, which="false")  #Generate weighted projection, with any agent attributes
    igraph::E(P)$oldweight <- igraph::E(P)$weight  #Save old edge weights
    P <- igraph::delete_edge_attr(P, "weight")  #Delete weight attribute

    backbone <- igraph::graph_from_adjacency_matrix(backbone, mode = "undirected", weighted = TRUE)  #Generate igraph backbone
    backbone <- igraph::as_data_frame(backbone, what = "edges")  #Generate backbone edgelist

    #==> Need to set a "weight" attribute in P for each edge listed in backbone

    return(P)
  }

}
