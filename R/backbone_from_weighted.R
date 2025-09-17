#' Extract the backbone from a weighted network
#'
#' \code{backbone_from_weighted()} extracts the unweighted backbone from a weighted network
#'
#' @param W A weighted network as a valued adjacency matrix or \link[Matrix]{Matrix}, or a weighted unipartite \link[igraph]{igraph} object
#' @param model string: backbone model, one of: \code{"disparity"}, \code{"lans"}, \code{"mlf"}, or \code{"global"}
#' @param alpha real: significance level of hypothesis test(s) in statistical models
#' @param signed boolean: return a signed backbone from a statistical model
#' @param mtc string: type of Multiple Test Correction, either \code{"none"} or a method allowed by [p.adjust()].
#' @param parameter real: parameter used to control structural backbone models
#' @param missing_as_zero boolean: treat missing edges as edges with zero weight and consider them for inclusion/exclusion in backbone
#' @param narrative boolean: display suggested text & citations
#' @param return string: return either only the \code{"backbone"} or \code{"everything"}
#'
#' @details
#' The \code{backbone_from_weighted} function extracts the backbone from a weighted unipartite network. The backbone is an unweighted
#' unipartite network that contains only edges whose weights are statistically significant (based on \code{alpha} for statistical models),
#' or which exhibit certain structural properties (based on \code{parameter} for structural models). For statistical models, when
#' \code{signed = FALSE}, the backbone contains edges that are statistically significantly strong under a one-tailed test. When
#' \code{signed = TRUE}, the backbone contains positive edges that are statistically significantly strong, and negative edges that are
#' statistically significantly weak, under a two-tailed test.
#'
#' The \code{model} parameter controls the model used to evaluate the edge weights. The available models include:
#'
#' *Statistical Models* (controlled by \code{alpha}, \code{signed}, and \code{mtc})
#' * \code{disparity} (default) - The disparity filter (Serrano et al., 2009)
#' * \code{lans} - Locally adaptive network sparsification (Foti et al., 2011)
#' * \code{mlf} - Marginal likelihood filter (Dianati, 2016)
#'
#' *Structural Models* (controlled by \code{parameter})
#' * \code{global} - \code{parameter} is a numeric vector of length 1 or 2. If \code{length(parameter)==1}, then edges with weights
#'   above \code{parameter} are preserved. If \code{length(parameter)==2}, then edges with weights above \code{max(parameter)} are
#'   preserved as positive, and edges with weights above \code{min(parameter)} are preserved as negative.
#'
#' The models implemented in \code{backbone_from_weighted()} can be applied to a weighted network that was obtained by projecting a
#' bipartite network or hypergraph. However, if the original bipartite network or hypergraph is available, it is better to use [backbone_from_projection()].
#'
#' @return If \code{return = "backbone"}, a backbone in the same class as \code{W}. If \code{return = "everything"}, then the backbone
#' is returned as an element in a list that also includes the original weighted network, (for statistical backbone models) the edgewise
#' p-values, a narrative description, and the original function call.
#'
#' @references package: {Neal, Z. P. (2025). backbone: An R Package to Extract Network Backbones. CRAN. \doi{10.32614/CRAN.package.backbone}}
#' @references disparity: {Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. *Proceedings of the National Academy of Sciences, 106*, 6483-6488. \doi{10.1073/pnas.0808904106}}
#' @references lans: {Foti, N. J., Hughes, J. M., & Rockmore, D. N. (2011). Nonparametric sparsification of complex multiscale networks. *PLOS One, 6*, e16431. \doi{10.1371/journal.pone.0016431}}
#' @references mlf: {Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. *Physical Review E, 93*, 012304. \doi{10.1103/PhysRevE.93.012304}}
#'
#' @export
#'
#' @examples
#' #A weighted network with heterogeneous (i.e. multiscale) weights
#' W <- matrix(c(0,10,10,10,10,75,0,0,0,0,
#'               10,0,1,1,1,0,0,0,0,0,
#'               10,1,0,1,1,0,0,0,0,0,
#'               10,1,1,0,1,0,0,0,0,0,
#'               10,1,1,1,0,0,0,0,0,0,
#'               75,0,0,0,0,0,100,100,100,100,
#'               0,0,0,0,0,100,0,10,10,10,
#'               0,0,0,0,0,100,10,0,10,10,
#'               0,0,0,0,0,100,10,10,0,10,
#'               0,0,0,0,0,100,10,10,10,0),10)
#'
#' W <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
#' plot(W, edge.width = sqrt(igraph::E(W)$weight)) #A stronger clique & a weaker clique
#'
#' mean_weight <- mean(igraph::E(W)$weight)  #Find average edge weight
#' bb <- backbone_from_weighted(W, model = "global", #A backbone with stronger-than-average edges...
#'       parameter = mean_weight)
#' plot(bb) #...ignores the weaker clique
#'
#' bb <- backbone_from_weighted(W, model = "disparity") #A disparity filter backbone...
#' plot(bb) #...preserves edges at multiple scales
backbone_from_weighted <- function(W,
                                   model = "disparity",
                                   alpha = 0.05,
                                   signed = FALSE,
                                   mtc = "none",
                                   parameter = 0,
                                   missing_as_zero = FALSE,
                                   narrative = FALSE,
                                   return = "backbone") {

  call <- match.call()

  #### Check parameters ####
  #All models
  if (!(model %in% c("disparity", "lans", "mlf", "global"))) {stop("`model` must be one of: \"disparity\", \"lans\", \"mlf\", or \"global\"")}
  if (!is.logical(missing_as_zero)) {stop("`missing_as_zero` must be either TRUE or FALSE")}
  if (!is.logical(narrative)) {stop("`narrative` must be either TRUE or FALSE")}
  if (!(return %in% c("backbone", "everything"))) {stop("`return` must be one of: \"backbone\", \"everything\"")}

  #Statistical models
  if (model %in% c("disparity", "lans", "mlf")) {
    if (!is.numeric(alpha)) {stop("`alpha` must be a numeric value between 0 and 1")}
    if (alpha < 0 | alpha > 1) {stop("`alpha` must be a numeric value between 0 and 1")}
    if (!is.logical(signed)) {stop("`signed` must be either TRUE or FALSE")}
    if (!(mtc %in% c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {stop("`mtc` must be one of: \"none\", \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", or \"fdr\"")}
  }

  #Structural models
  if (model %in% c("global")) {
    if (!is.numeric(parameter)) {stop("parameter must be a numeric vector of length 1 or 2")}
    if (length(parameter)<1 | length(parameter)>2) {stop("parameter must be a numeric vector of length 1 or 2")}
  }

  #### Check and format input ####
  #Check that input is a weighted adjacency matrix or weighted unipartite igraph
  if (!methods::is(W,"matrix") & !methods::is(W,"Matrix") & !methods::is(W,"igraph")) {stop("`W` must be an adjacency matrix or Matrix, or an igraph object")}

  if (methods::is(W,"matrix")) {
    if (dim(as.matrix(W))[1] != dim(as.matrix(W))[2]) {stop("`W` must be a square adjacency matrix")}
    if (all(as.matrix(W) %in% c(0,1))) {stop("The entries of `W` must represent edge weights")}
  }

  if (methods::is(W,"igraph")) {
    if (igraph::is_bipartite(W)) {stop("`W` must be a unipartite igraph object")}
    if (!"weight" %in% igraph::edge_attr_names(W)) {stop("`W` must contain an edge weight attribute")}
    }

  #Convert input to adjacency matrix
  if (methods::is(W,"matrix")) {A <- W}  #matrix --> matrix
  if (methods::is(W,"Matrix")) {A <- as.matrix(W)}  #Matrix --> matrix
  if (methods::is(W,"igraph")) {A <- igraph::as_adjacency_matrix(W, names = FALSE, sparse = FALSE, attr = "weight")}

  #### Statistical Models ####
  if (model == "disparity") {p <- .disparity(A, missing_as_zero, signed)}
  if (model == "lans") {p <- .lans(A, missing_as_zero, signed)}
  if (model == "mlf") {p <- .mlf(A, missing_as_zero, signed)}
  if (model == "disparity" | model == "lans" | model == "mlf") {backbone <- .retain(p, alpha, mtc)}

  #### Structural Models ####
  if (model == "global") {backbone <- .global(A, missing_as_zero, parameter)}

  #### Construct narrative ####
  # First sentence (descriptive)
  if (signed & (model == "disparity" | model == "lans" | model == "mlf")) {type <- "signed"} else {type <- "unweighted"}
  if (model == "global" & length(parameter)==2) {type <- "signed"} else {type <- "unweighted"}

  text <- paste0("We used the backbone package for R (v", utils::packageVersion("backbone"), "; Neal, 2025) to extract the ", type, " backbone of a weighted network containing ", nrow(A), " nodes.")

  # Second sentence (model and outcome)
  if (mtc == "none") {correction <- ""}
  if (mtc == "bonferroni") {correction <- ", Bonferroni adjusted"}
  if (mtc == "holm") {correction <- ", Holm adjusted"}
  if (mtc == "hommel") {correction <- ", Hommel adjusted"}
  if (mtc == "hochberg") {correction <- ", Hochberg adjusted"}
  if (mtc == "BH" | mtc == "fdr") {correction <- ", Benjamini & Hochberg adjusted"}
  if (mtc == "BY") {correction <- ", Benjamini & Yekutieli adjusted"}

  if (model == "disparity") {desc <- "the disparity filter (Serrano et al., 2009)"}
  if (model == "lans") {desc <- "locally adaptive network sparsification (LANS; Foti et al., 2011)"}
  if (model == "mlf") {desc <- "the marginal likelihood filter (MLF; Dianati, 2016)"}

  old <- sum(A!=0, na.rm=TRUE)  #Number of edges in weighted network
  new <- sum(backbone!=0)  #Number of edges in backbone
  reduced_edges <- round(((old - new) / old)*100,2)

  if (model == "disparity" | model == "lans" | model == "mlf") {text <- paste0(text, " An edge was retained in the backbone if its weight was statistically significant (alpha = ", alpha, correction, ") using ", desc, ", which reduced the number of edges by ", reduced_edges, "%.")}
  if (model == "global" & length(parameter)==1) {text <- paste0(text, " An edge was retained in the backbone if its weight was larger than ", max(parameter), ", which reduced the number of edges by ", reduced_edges, "%.")}
  if (model == "global" & length(parameter)==2) {text <- paste0(text, " An edge was retained in the backbone as positive if its weight was larger than ", max(parameter), " and as negative if its weight was smaller than ", min(parameter), " which reduced the number of edges by ", reduced_edges, "%.")}

  # References
  text <- paste0(text, "\n\nNeal, Z. P. 2025. backbone: An R Package to Extract Network Backbones. CRAN. https://doi.org/10.32614/CRAN.package.backbone")
  if (model == "disparity") {text <- paste0(text, "\n\nSerrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106, 6483-6488. https://doi.org/10.1073/pnas.0808904106")}
  if (model == "lans") {text <- paste0(text, "\n\nFoti, N. J., Hughes, J. M., & Rockmore, D. N. (2011). Nonparametric sparsification of complex multiscale networks. PLOS One, 6, e16431. https://doi.org/10.1371/journal.pone.0016431")}
  if (model == "mlf") {text <- paste0(text, "\n\nDianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. Physical Review E, 93, 012304. https://doi.org/10.1103/PhysRevE.93.012304")}

  #### Display narrative ####
  if (narrative) {message(text)}

  #### Prepare backbone ####
  if (methods::is(W,"matrix")) {
    rownames(backbone) <- rownames(W)
    colnames(backbone) <- rownames(W)
  }

  if (methods::is(W,"Matrix")) {
    rownames(backbone) <- rownames(W)
    colnames(backbone) <- rownames(W)
    backbone <- Matrix::Matrix(backbone)
  }

  if (methods::is(W,"igraph")) {
    temp <- W  #Placeholder for backbone
    igraph::E(temp)$oldweight <- igraph::E(temp)$weight  #Save old edge weights
    temp <- igraph::delete_edge_attr(temp, "weight")  #Delete weight attribute
    temp <- igraph::set_edge_attr(temp, "sign", value = backbone[igraph::as_edgelist(temp, names = FALSE)])  #Insert edge retention marker as attribute
    temp <- igraph::delete_edges(temp, which(igraph::E(temp)$sign==0))  #Delete any edges that should not be retained
    if (!signed & (model == "disparity" | model == "lans" | model == "mlf")) {temp <- igraph::delete_edge_attr(temp, "sign")}  #If backbone is not signed, remove edge retention marker
    if (length(parameter)!=2 & (model == "global")) {temp <- igraph::delete_edge_attr(temp, "sign")}  #If backbone is not signed, remove edge retention marker
    backbone <- temp
  }

  #### Return ####
  if (return == "backbone") {return(backbone)}
  if (return == "everything" & (model == "disparity" | model == "lans" | model == "mlf")) {return(list(weighted = W, backbone = backbone, pvalues = p, narrative = text, call = call))}
  if (return == "everything" & (model == "global")) {return(list(weighted = W, backbone = backbone, narrative = text, call = call))}
}
