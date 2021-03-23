#' backbone: Extracts the Backbone from Weighted Graphs
#'
#' @description Provides methods for extracting from a weighted graph
#'     a binary or signed backbone that retains only the significant edges.
#'     The user may input a weighted graph, or a bipartite graph
#'     from which a weighted graph is first constructed via projection.
#'     Backbone extraction methods include:
#'     \itemize{
#'     \item the stochastic degree sequence model ({Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. Social Networks, 39, Elsevier: 84-97.}\doi{10.1016/j.socnet.2014.06.001}),
#'     \item hypergeometric model ({Neal, Zachary. 2013. “Identifying Statistically Significant Edges in One-Mode Projections.” Social Network Analysis and Mining 3 (4). Springer: 915–24.}\doi{10.1007/s13278-013-0107-y}),
#'     \item Poisson binomial model ({Neal, Sagan, Domagalski, upcoming})
#'     \item the fixed degree sequence model ({Zweig, Katharina Anna, and Michael Kaufmann. 2011. “A Systematic Approach to the One-Mode Projection of Bipartite Graphs.” Social Network Analysis and Mining 1 (3): 187–218.}\doi{10.1007/s13278-011-0021-0}),
#'     \item  as well as a universal threshold method.
#'     }
#'
#' @details To decide whether an edge of a bipartite projection \eqn{B*t(B)} is statistically significant, we compare the edge's observed weight to the distribution of weights expected in a projection obtained from a random bipartite graph under a null model.
#'     Given a null model that can generate a distribution of an edge's weights, we can then compute the probability that the observed weight is in the upper or lower tail of that distribution.
#'     The `backbone` package provides several different null models to use, all of which define different constraints on the row degrees and column degrees of a random bipartite graph in their corresponding distribution.
#' \itemize{
#' \item '\code{\link{sdsm}}': computes the probability of edge weights being above or below the observed edge weights in a bipartite projection using the stochastic degree sequence model.
#'     Here the expected row and column degrees of a random bipartite graph in the distribution match that of the original bipartite graph.
#' \item '\code{\link{fixedrow}}': computes the probability of edge weights being above or below the observed edge weights in a bipartite projection using the hypergeometric model.
#'     Here the row degrees of a random bipartite graph in the distribution exactly match that of the original bipartite graph. The column degrees are allowed to vary.
#' \item '\code{\link{fixedcol}}': computes the probability of edge weights being above or below the observed edge weights in a bipartite projection using the Poisson binomial model.
#'     Here the column degrees of a random bipartite graph in the distribution exactly match that of the original bipartite graph. The row degrees are allowed to vary.
#' \item '\code{\link{fdsm}}': computes the proportion of edge weights above or below the observed edge weights in a bipartite projection using the fixed degree sequence model.
#'     Here the row and column degrees of a random bipartite graph in the distribution exactly match that of the original bipartite graph.
#' \item '\code{\link{fixedfill}}':computes the probability of edge weights being above or below the observed edge weights in a bipartite projection where the number of edges in a random bipartite graph
#'     of the distribution exactly match the number of edges of the original bipartite graph.
#' }
#' Once one of these models is computed, use \code{\link{backbone.extract}} to return the backbone graph for a chosen alpha value.
#'
#' In addition to these methods, one can also use \code{\link{universal}} to return a backbone graph in which edge weights
#'     are set to 1 if above the given upper parameter threshold,
#'     and set to -1 if below the given lower parameter threshold, and are 0 otherwise.
#'
#' @details Additional functions that aid in the use of the above models are exported:
#' \itemize{
#' \item '\code{\link{bicm}}': finds a matrix that maximizes the entropy function, used in \code{\link{sdsm}}.
#' \item '\code{\link{curveball}}': generates a random 0/1 matrix with the same row and column sums as the input, used in \code{\link{fdsm}}.
#' }
#'
#' @details For additional documentation and background on the package functions, see \href{../doc/backbone.html}{\code{vignette("backbone", package = "backbone")}}.
#' @references {Domagalski, R., Neal, Z. P., and Sagan, B. (2021). backbone: An R Package for Backbone Extraction of Weighted Graphs. PLoS ONE. \doi{10.1371/journal.pone.0244363}}
#'
#' @docType package
#' @name backbone
NULL

