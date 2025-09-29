#' Extract the backbone from a weighted bipartite or hypgergraph projection
#'
#' \code{backbone_from_projection()} extracts the unweighted backbone from the weighted projection of a bipartite network or hypergraph
#'
#' @param B An unweighted bipartite network or hypergraph as an incidence matrix or \link[Matrix]{Matrix}, or as a bipartite \link[igraph]{igraph} object
#' @param alpha real: significance level of hypothesis test(s)
#' @param model string: backbone model, one of: \code{"sdsm"}, \code{"fdsm"}, \code{"fixedrow"}, \code{"fixedcol"}, or \code{"fixedfill"}
#' @param signed logical: return a signed backbone
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param missing_as_zero logical: treat missing edges as edges with zero weight and test them for significance
#' @param narrative logical: display suggested text & citations
#' @param trials numeric: if \code{model = "fdsm"}, the number of graphs generated using fastball to approximate the edge weight distribution
#' @param backbone_only logical: return just the backbone (default), or a detailed backbone object
#'
#' @details
#' The \code{backbone_from_projection} function extracts the backbone from the weighted projection of a bipartite network or hypergraph.
#' The backbone is an unweighted unipartite network of agents that contains only edges whose weights in the projection are statistically
#' significant. When \code{signed = FALSE}, the backbone contains edges that are statistically significantly strong under a one-tailed test.
#' When \code{signed = TRUE}, the backbone contains positive edges that are statistically significantly strong, and negative edges that are
#' statistically significantly weak, under a two-tailed test.
#'
#' The \code{model} parameter controls the null model used to evaluate the statistical significance of edge weights. All available models
#' are *statistical models* that are controlled by \code{alpha}, and differ in the constraints they impose on \code{B}:
#' * \code{sdsm} (default) - The "Stochastic Degree Sequence Model" (SDSM; Neal et al., 2021) approximately constrains the agent and artifact degrees, and exactly constrains edges that are prohibited (weight = 10) or required (weight = 11; Neal & Neal, 2023)
#' * \code{fdsm} - The "Fixed Degree Sequence Model" (Neal et al., 2021) exactly constrains the agent and artifact degrees
#' * \code{fixedfill} - The "fixed fill" model (Neal et al., 2021) exactly constrains the total number of edges (i.e., sum)
#' * \code{fixedrow} - The "fixed row" model (Neal et al., 2021) exactly constrains the agent degrees (i.e., row sums)
#' * \code{fixedcol} - The "fixed column" model (Neal et al., 2021) exactly constrains the artifact degrees (i.e., column sums)
#'
#' Although \code{backbone_from_projection} extracts the backbone from the weighted projection of a bipartite network or hypergraph,
#' the input \code{B} *must be the bipartite network or hypergraph itself, and not the weighted projection*. This is necessary
#' because these backbone models use information in the bipartite network that is missing from the projection. The "agent" nodes that
#' appear in the projection must be represented by rows if \code{B} is an incidence matrix, or by \code{type = FALSE} nodes if \code{B}
#' is a bipartite igraph object. In either case, the source network must be binary (i.e., unweighted), unless \code{model = "sdsm"},
#' when "prohibited" edges can be represented with weight = 10 and "required" edges can be represented with weight = 11.
#'
#' @return A backbone in the same class as \code{B}, or if \code{backbone_only = FALSE}, then a backbone object.
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references sdsm-ec model: {Neal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. *International Conference on Complex Networks and Their Applications, 12*, 127-136. \doi{10.1007/978-3-031-53468-3_11}}
#' @references all other models: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @export
#'
#' @examples
#' #A binary bipartite network of 30 agents & 75 artifacts
#' #The agents form three communities
#' B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
#'                  matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.2),10)),
#'            cbind(matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.8),10),
#'                  matrix(rbinom(250,1,.2),10)),
#'            cbind(matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.2),10),
#'                  matrix(rbinom(250,1,.8),10)))
#' B <- igraph::graph_from_biadjacency_matrix(B)
#'
#' P <- igraph::bipartite_projection(B, which = "true")  #An ordinary weighted projection...
#' plot(P)                                               #...is a dense hairball
#'
#' bb <- backbone_from_projection(B)  #A backbone...
#' plot(bb)                           #...is sparse with clear communities
backbone_from_projection <- function(B,
                                     alpha = 0.05,
                                     model = "sdsm",
                                     signed = FALSE,
                                     mtc = "none",
                                     missing_as_zero = FALSE,
                                     narrative = FALSE,
                                     trials = NULL,
                                     backbone_only = TRUE) {

  call <- match.call()

  #### Check parameters and input ####
  I <- .check_and_coerce(N = B, source = "projection", model = model, alpha = alpha, signed = signed, mtc = mtc, missing_as_zero = missing_as_zero, narrative = narrative, trials = trials, backbone_only = backbone_only)
  if (model == "sdsm" & any(I %in% c(10,11))) {model <- "sdsm_ec"}

  #### Compute p-values ####
  if (model == "sdsm") {p <- .sdsm(I, missing_as_zero, signed)}
  if (model == "sdsm_ec") {p <- .sdsm_ec(I, missing_as_zero, signed)}
  if (model == "fixedrow") {p <- .fixedrow(I, missing_as_zero, signed)}
  if (model == "fixedcol") {p <- .fixedcol(I, missing_as_zero, signed)}
  if (model == "fixedfill") {p <- .fixedfill(I, missing_as_zero, signed)}
  if (model == "fdsm") {p <- .fdsm(I, missing_as_zero, signed, alpha, mtc, trials)}

  #### Retain edges ####
  backbone <- .retain(p, alpha, mtc)

  #### Construct narrative ####
  # First sentence (descriptive)
  if (signed) {type <- "signed"} else {type <- "unweighted"}

  text <- paste0("The backbone package for R (v", utils::packageVersion("backbone"), "; Neal, 2025) was used to extract the ", type, " backbone of the weighted projection of a bipartite network containing ", nrow(I), " agents and ", ncol(I), " artifacts.")

  # Second sentence (model and outcome)
  if (mtc == "none") {correction <- ""}
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
  if (model == "sdsm_ec") {desc <- "the stochastic degree sequence model with edge constraints (SDSM-EC; Neal & Neal, 2023)"}
  if (model == "fdsm") {desc <- paste0("the fixed degree sequence model (FDSM; Neal, Domagalski, and Sagan, 2021)")}

  old <- sum(p$upper!=0, na.rm=TRUE)  #Number of edges in projection (i.e., number of edges tested, and that have an upper-tail p-value)
  new <- sum(backbone!=0)  #Number of edges in backbone
  reduced_edges <- round(((old - new) / old)*100,2)

  text <- paste0(text, " An edge was retained in the backbone if its weight was statistically significant (alpha = ", alpha, correction, ") using ", desc, ", which reduced the number of edges by ", reduced_edges, "%.")

  #References
  text <- paste0(text, "\n\nNeal, Z. P. 2025. backbone: An R Package to Extract Network Backbones. CRAN. https://doi.org/10.32614/CRAN.package.backbone")
  if (model %in% c("sdsm", "fdsm", "fixedrow", "fixedcol", "fixedfill")) {text <- paste0(text, "\n\nNeal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "sdsm_ec") {text <- paste0(text, "\n\nNeal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. International Conference on Complex Networks and Their Applications, 12, 127-136. https://doi.org/10.1007/978-3-031-53468-3_11")}

  #### Display narrative ####
  if (narrative) {message(text)}

  #### Prepare backbone ####
  if (methods::is(B,"matrix")) {
    rownames(backbone) <- rownames(B)
    colnames(backbone) <- rownames(B)
    P <- tcrossprod(I)
  }

  if (methods::is(B,"Matrix")) {
    rownames(backbone) <- rownames(B)
    colnames(backbone) <- rownames(B)
    backbone <- Matrix::Matrix(backbone)
    P <- tcrossprod(I)
    P <- Matrix::Matrix(P)
  }

  if (methods::is(B,"igraph")) {
    tempB <- B  #Temporary Bipartite
    if (model=="sdsm_ec") {tempB <- igraph::delete_edges(tempB, which(igraph::E(tempB)$weight==10))}  #If there are prohibited edges in an igraph object, remove them
    P <- igraph::bipartite_projection(tempB, which="false")  #Generate weighted projection, with any agent attributes
    tempP <- P  #Placeholder for backbone
    igraph::E(tempP)$oldweight <- igraph::E(tempP)$weight  #Save old edge weights
    tempP <- igraph::delete_edge_attr(tempP, "weight")  #Delete weight attribute
    for (attr in igraph::vertex_attr_names(tempP)) {if (all(is.na(igraph::vertex_attr(tempP, attr)))) {tempP <- igraph::delete_vertex_attr(tempP, attr)}}  #Delete attributes of artifact nodes
    tempP <- igraph::set_edge_attr(tempP, "sign", value = backbone[igraph::as_edgelist(tempP, names = FALSE)])  #Insert edge retention marker as attribute
    tempP <- igraph::delete_edges(tempP, which(igraph::E(tempP)$sign==0))  #Delete any edges that should not be retained
    if (!signed) {tempP <- igraph::delete_edge_attr(tempP, "sign")}  #If backbone is not signed, remove edge retention marker
    backbone <- tempP
    if (!is.null(backbone$name)) {backbone$name <- paste0(model, " backbone of ", backbone$name)}
    if (is.null(backbone$name)) {backbone$name <- paste0(model, " backbone")}
    backbone$call <- call
    backbone$narrative <- text
    }

  #### Return ####
  if (backbone_only) {return(backbone)}
  if (!backbone_only) {return(structure(list(bipartite = B, projection = P, backbone = backbone, pvalues = p, narrative = text, model = model, alpha = alpha, call = call), class = "backbone"))}
}
