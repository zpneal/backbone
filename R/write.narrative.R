#' Generates suggested manuscript text
#'
#' `write.narrative` displays suggested text and references that can be used in a manuscript to describe the extracted backbone
#'
#' @param agents integer: number of agents in a bipartite graph, or nodes in a unipartite graph
#' @param artifacts integer: number of artifacts in a bipartite graph
#' @param weighted boolean: TRUE if input graph was weighted
#' @param bipartite boolean: TRUE if input graph was bipartite
#' @param symmetric boolean: TRUE if input graph was symmetric
#' @param signed boolean: TRUE if a signed backbone was requested
#' @param mtc string: type of multiple test correction
#' @param alpha numeric: alpha significance threshold (used in statistical models)
#' @param s numeric: sparsification parameter (used in sparsification models)
#' @param ut numeric: upper threshold (used in global threshold)
#' @param lt numeric: lower threshold (used in global threshold)
#' @param trials integer: number of trials used to estimate FDSM or oSDSM p-values
#' @param escore string: Method for scoring edges' importance
#' @param normalize string: Method for normalizing edge scores
#' @param filter string: Type of filter to apply
#' @param umst boolean: TRUE if the backbone should include the union of minimum spanning trees, to ensure connectivity
#' @param model string: name of backbone model
#' @param reduced_edges numeric: percent reduction in number of edges
#' @param reduced_nodes numeric: percent reduction in number of connected nodes
#'
#' @return NULL; only displays text in the console
#' @keywords internal
write.narrative <- function(agents = 0, artifacts = 0, weighted = FALSE, bipartite = FALSE, symmetric = FALSE,
                            signed = FALSE, mtc = "none", alpha = NULL, s = NULL, ut = NULL, lt = NULL, trials = 0,
                            escore = NULL, normalize = NULL, filter = NULL, umst = NULL, #Recipe for custom sparsify
                            model = "", reduced_edges = NULL, reduced_nodes = NULL) {

  #### Prepare narrative components ####
  version <- utils::packageVersion("backbone")

  #Content of starting object
  if (bipartite) {contents <- paste0(agents, " agents and ", artifacts, " artifacts")} else {contents <- paste0(agents, " nodes")}
  if (weighted & bipartite) {type <- "the weighted projection of a weighted bipartite"}
  if (!weighted & bipartite) {type <- "the weighted projection of an unweighted bipartite"}
  if (weighted & symmetric & !bipartite) {type <- "a weighted and undirected unipartite"}
  if (weighted & !symmetric & !bipartite) {type <- "a weighted and directed unipartite"}
  if (!weighted & symmetric & !bipartite) {type <- "an unweighted and undirected unipartite"}
  if (!weighted & !symmetric & !bipartite) {type <- "an unweighted and directed unipartite"}
  if (signed) {signed <- "signed"} else {signed <- "unweighted"}
  correction <- ""

  #Multiple test correction
  if (mtc == "bonferroni") {correction <- ", Bonferroni adjusted"}
  if (mtc == "holm") {correction <- ", Holm adjusted"}
  if (mtc == "hommel") {correction <- ", Hommel adjusted"}
  if (mtc == "hochberg") {correction <- ", Hochberg adjusted"}
  if (mtc == "BH") {correction <- ", Benjamini & Hochberg adjusted"}
  if (mtc == "BY") {correction <- ", Benjamini & Yekutieli adjusted"}

  #Type of model
  if (model == "fixedfill") {desc <- "the fixed fill model (FFM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "fixedrow") {desc <- "the fixed row model (FRM; Neal, 2013)"}
  if (model == "fixedcol") {desc <- "the fixed column model (FCM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "sdsm") {desc <- "the stochastic degree sequence model (SDSM; Neal, 2014)"}
  if (model == "fdsm") {desc <- paste0("the fixed degree sequence model (FDSM; Neal, 2014), where p-values were estimated from ", trials, " Monte Carlo trials")}
  if (model == "osdsm") {desc <- "the ordinal stochastic degree sequence model (oSDSM; Neal, 2017)"}
  if (model == "disparity") {desc <- "the disparity filter (Serrano et al., 2009)"}
  if (model == "mlf") {desc <- "the marginal likelihood filter (Dianati, 2016)"}
  if (model == "lans") {desc <- "locally adaptive network sparsification (Foti et al., 2011)"}
  if (model == "skeleton") {desc <- "Karger's (1994) skeleton model"}
  if (model == "gspar") {desc <- "Satuluri et al.'s (2011) G-Spar model"}
  if (model == "lspar") {desc <- "Satuluri et al.'s (2011) L-Spar model"}
  if (model == "simmelian") {desc <- "Nick et al.'s (2013) Simmelian model"}
  if (model == "jaccard") {desc <- "Goldberg and Roth's (2003) Jaccard model"}
  if (model == "meetmin") {desc <- "Goldberg and Roth's (2003) MeetMin model"}
  if (model == "geometric") {desc <- "Goldberg and Roth's (2003) geometric model"}
  if (model == "hypergeometric") {desc <- "Goldberg and Roth's (2003) hypergeometric model"}
  if (model == "degree") {desc <- "Hamann et al.'s (2016) local degree model"}
  if (model == "quadrilateral") {desc <- "Nocaj et al.'s (2015) quadrilateral Simmelian model"}

  if (model == "sparsify") {
    desc <- paste0("In this model, edge weights are defined as ", escore, ", ")
    if (normalize == "none") {desc <- paste0(desc, "edges weights are not normalized, ")} else {desc <- paste0(desc, "edge weights are normalized using ", normalize, ", ")}
    desc <- paste0(desc, "and edges are filtered using ", filter, ". ")
    if (umst) {desc <- paste0(desc, "To ensure a connected backbone, edges in the union of minimum spanning trees were also included.")}
    }

  #### First sentence (descriptive) ####
  text <- paste0("We used the backbone package for R (v", version, "; Neal, 2022) to extract the ", signed, " backbone")
  text <- paste0(text, " of ", type, " network containing ", contents, ".")

  #### Second sentence (model) ####
  #Statistical models
  if (!is.null(alpha)) {
    text <- paste0(text, " An edge was retained in the backbone if its weight was statistically significant (alpha = ", alpha, correction, ") using ", desc, ".")
    }

  #Sparsify models
  if (!is.null(s)) {
    if (model != "sparsify") {text <- paste0(text, " Specifically, we used ", desc, " with a sparsification threshold of ", s, ".")}
    if (model == "sparsify") {text <- paste0(text, " Specifically, we used an unweighted graph sparsification model with a sparsification threshold parameter of ", s, ". ", desc)}
  }

  #Global threshold
  if (model == "global" & is.null(lt)) {  #Binary
    text <- paste0(text, " An edge was retained if its weight was greater than ", round(ut, 3), ".")
  }
  if (model == "global" & !is.null(lt)) {  #Signed
    text <- paste0(text, " An edge was retained as positive if its weight was greater than ", round(ut, 3), ", and as negative if its weight was less than ", round(lt, 3), ".")
  }

  #### Third sentence (reduction) ####
  text <- paste0(text, " This reduced the number of edges by ", reduced_edges, "%, and reduced the number of connected nodes by ", reduced_nodes, "%.")

  #### Display text ####
  message("")
  message("=== Suggested manuscript text and citations ===")
  message(text)
  message("")
  message("Neal, Z. P. 2022. backbone: An R Package to Extract Network Backbones. PLOS ONE, 17, e0269137. https://doi.org/10.1371/journal.pone.0269137")
  message("")
  if (model == "fixedfill") {message("Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "fixedrow") {message("Neal. Z. P. (2013). Identifying statistically significant edges in one-mode projections. Social Network Analysis and Mining, 3, 915-924. https://doi.org/10.1007/s13278-013-0107-y")}
  if (model == "fixedcol") {message("Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "sdsm") {message("Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. Social Networks, 39, 84-97. https://doi.org/10.1016/j.socnet.2014.06.001")}
  if (model == "fdsm") {message("Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. Social Networks, 39, 84-97. https://doi.org/10.1016/j.socnet.2014.06.001")}
  if (model == "osdsm") {message("Neal, Z. P. (2017). Well connected compared to what? Rethinking frames of reference in world city network research. Environment and Planning A, 49, 2859-2877. https://doi.org/10.1177/0308518X16631339")}
  if (model == "disparity") {message("Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488. https://doi.org/10.1073/pnas.0808904106")}
  if (model == "mlf") {message("Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. Physical Review E, 93, 012304. https://doi.org/10.1103/PhysRevE.93.012304")}
  if (model == "lans") {message("Foti, N. J., Hughes, J. M., and Rockmore, D. N. (2011). Nonparametric Sparsification of Complex Multiscale Networks. PLOS One, 6(2), e16431. https://doi.org/10.1371/journal.pone.0016431")}
  if (model == "skeleton") {message("Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. Mathematics of Operations Research, 24(2), 383-413. https://doi.org/10.1287/moor.24.2.383")}
  if (model == "gspar") {message("Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). https://doi.org/10.1145/1989323.1989399")}
  if (model == "lspar") {message("Satuluri, V., Parthasarathy, S., & Ruan, Y. (2011, June). Local graph sparsification for scalable clustering. In Proceedings of the 2011 ACM SIGMOD International Conference on Management of data (pp. 721-732). https://doi.org/10.1145/1989323.1989399")}
  if (model == "simmelian") {message("Nick, B., Lee, C., Cunningham, P., & Brandes, U. (2013, August). Simmelian backbones: Amplifying hidden homophily in facebook networks. In Proceedings of the 2013 IEEE/ACM international conference on advances in social networks analysis and mining (pp. 525-532). https://doi.org/10.1145/2492517.2492569")}
  if (model == "jaccard") {message("Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. Proceedings of the National Academy of Sciences, 100(8), 4372-4376. https://doi.org/10.1073/pnas.0735871100")}
  if (model == "meetmin") {message("Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. Proceedings of the National Academy of Sciences, 100(8), 4372-4376. https://doi.org/10.1073/pnas.0735871100")}
  if (model == "geometric") {message("Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. Proceedings of the National Academy of Sciences, 100(8), 4372-4376. https://doi.org/10.1073/pnas.0735871100")}
  if (model == "hypergeometric") {message("Goldberg, D. S., & Roth, F. P. (2003). Assessing experimentally derived interactions in a small world. Proceedings of the National Academy of Sciences, 100(8), 4372-4376. https://doi.org/10.1073/pnas.0735871100")}
  if (model == "degree") {message("Hamann, M., Lindner, G., Meyerhenke, H., Staudt, C. L., & Wagner, D. (2016). Structure-preserving sparsification methods for social networks. Social Network Analysis and Mining, 6(1), 22. https:://10.1007/s13278-016-0332-2")}
  if (model == "quadrilateral") {message("Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. Journal of Graph Algorithms and Applications: JGAA, 19(2), 595-618. https://10.7155/jgaa.00370")}
}
