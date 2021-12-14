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
#' @param fwer string: type of familywise error rate correction
#' @param alpha numeric: alpha significance threshold (used in statistical models)
#' @param s numeric: sparsification parameter (used in sparsification models)
#' @param ut numeric: upper threshold (used in global threshold)
#' @param lt numeric: lower threshold (used in global threshold)
#' @param trials: integer: number of trials used to estimate FDSM or oSDSM p-values
#' @param model string: name of backbone null model
#' @param retained numeric: percent of edges retained in the backbone
#'
#' @return NULL; only displays text in the console
#' @keywords internal
write.narrative <- function(agents, artifacts, weighted, bipartite, symmetric, signed, fwer, alpha, s, ut, lt,  trials, model, retained) {

  #### Prepare narrative components ####
  version <- utils::packageVersion("backbone")
  if (bipartite) {contents <- paste0(agents, " agents and ", artifacts, " artifacts")} else {contents <- paste0(agents, " nodes")}
  if (weighted & bipartite) {type <- "the weighted projection of a weighted bipartite"}
  if (!weighted & bipartite) {type <- "the weighted projection of an unweighted bipartite"}
  if (weighted & symmetric & !bipartite) {type <- "a weighted and undirected unipartite"}
  if (weighted & !symmetric & !bipartite) {type <- "a weighted and directed unipartite"}
  if (!weighted & symmetric & !bipartite) {type <- "an unweighted and undirected unipartite"}
  if (!weighted & !symmetric & !bipartite) {type <- "an unweighted and directed unipartite"}
  if (signed) {signed <- "signed"} else {signed <- "binary"}
  correction <- ""
  if (fwer == "bonferroni") {correction <- ", Bonferroni adjusted"}
  if (fwer == "holm") {correction <- ", Holm adjusted"}
  if (fwer == "hommel") {correction <- ", Hommel adjusted"}
  if (fwer == "hochberg") {correction <- ", Hochberg adjusted"}
  if (fwer == "BH") {correction <- ", Benjamini & Hochberg adjusted"}
  if (fwer == "BY") {correction <- ", Benjamini & Yekutieli adjusted"}

  if (model == "fixedfill") {desc <- "the fixed fill model (FFM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "fixedrow") {desc <- "the fixed row model (FRM; Neal, 2013)"}
  if (model == "fixedcol") {desc <- "the fixed column model (FCM; Neal, Domagalski, and Sagan, 2021)"}
  if (model == "sdsm") {desc <- "the stochastic degree sequence model (SDSM; Neal, 2014)"}
  if (model == "fdsm") {desc <- paste0("the fixed degree sequence model (FDSM; Neal, 2014), where p-values were estimated from ", trials, " Monte Carlo trials")}
  if (model == "osdsm") {desc <- "the ordinal stochastic degree sequence model (oSDSM; Neal, 2017)"}
  if (model == "disparity") {desc <- "the disparity filter (Serrano et al., 2009)"}
  if (model == "sparsify") {desc <- "a sparsification model"}
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

  #### First sentence (descriptive) ####
  text <- paste0("We used the backbone package for R (v", version, "; Domagalski, Neal, & Sagan, 2021) to extract the ", signed, " backbone")
  text <- paste0(text, " of ", type, " network containing ", contents, ".")

  #### Second sentence (model) ####
  #Statistical models
  if (!is.null(alpha)) {
    text <- paste0(text, " An edge was retained in the backbone if its weight was statistically significant")
    text <- paste0(text, " (alpha = ", alpha, correction, ") using ", desc,", which retained ", retained, "% of edges.")
    }

  #Sparsify models
  if (!is.null(s)) {
    text <- paste0(text, " Specifically, we used ", desc, " with a sparsification threshold of ", s, ", which retained ", retained, "% of edges.")
  }

  #Global threshold
  if (model == "global" & is.null(lt)) {  #Binary
    text <- paste0(text, " An edge was retained if its weight was greater than ", round(ut, 3), ", resulting in the retention of ", retained, "% of edges.")
  }
  if (model == "global" & !is.null(lt)) {  #Signed
    text <- paste0(text, " An edge was retained as positive if its weight was greater than ", round(ut, 3), ", and as negative if its weight was less than ", round(lt, 3), ", resulting in the retention of ", retained, "% of edges.")
  }

  #### Display text ####
  message("")
  message("=== Suggested manuscript text and citations ===")
  message(text)
  message("")
  message("Domagalski, R., Neal, Z. P., and Sagan, B. (2021). backbone: An R Package for Backbone Extraction of Weighted Graphs. PLoS ONE, 16, e0244363. https://doi.org/10.1371/journal.pone.0244363")
  message("")
  if (model == "fixedfill") {message("Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "fixedrow") {message("Neal. Z. P. (2013). Identifying statistically significant edges in one-mode projections. Social Network Analysis and Mining, 3, 915-924. https://doi.org/10.1007/s13278-013-0107-y")}
  if (model == "fixedcol") {message("Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. Scientific Reports, 11, 23929. https://doi.org/10.1038/s41598-021-03238-3")}
  if (model == "sdsm") {message("Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. Social Networks, 39, 84-97. https://doi.org/10.1016/j.socnet.2014.06.001")}
  if (model == "fdsm") {message("Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. Social Networks, 39, 84-97. https://doi.org/10.1016/j.socnet.2014.06.001")}
  if (model == "osdsm") {message("Neal, Z. P. (2017). Well connected compared to what? Rethinking frames of reference in world city network research. Environment and Planning A, 49, 2859-2877. https://doi.org/10.1177/0308518X16631339")}
  if (model == "disparity") {message("Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488. https://doi.org/10.1073/pnas.0808904106")}
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



