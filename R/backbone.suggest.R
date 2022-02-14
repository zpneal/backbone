#' Suggest a backbone model
#'
#' `backbone.suggest` suggests and optionally runs an appropriate backbone model for a graph object.
#'
#' @param G graph: A graph represented in an object of class matrix, sparse \code{\link{Matrix}}, dataframe, \code{\link{igraph}}, or \code{\link{network}}.
#' @param s numeric: If provided, a backbone is extracted using this value as the significance level or sparsification parameter.
#'
#' @return If `s` == NULL: NULL, but a message is displayed with a suggested model.
#'    If 0 <= `s` <= 1: A binary backbone graph in the same class as `G`, obtained by extracting the backbone
#'    at the `s` significance level (if a statistical model is suggested) or using sparisfication parameter `s`
#'    (if a sparsification model is suggested). The code used to perform the extraction, and suggested manuscript
#'    text are displayed.
#' @export
#'
#' @examples
#' M <- matrix(runif(100),10,10)  #A random weighted, directed graph
#' backbone <- backbone.suggest(M)
#' backbone <- backbone.suggest(M, s = 0.05)
backbone.suggest <- function(G, s = NULL) {

  #### Parameter check, Convert supplied object ####
  if (!is.null(s)) {if (!is.numeric(s) | s < 0 | s > 1) {stop("If supplied, s must be between 0 and 1.")}}
  if (is.null(s)) {G <- tomatrix(G)} else {G <- suppressMessages(tomatrix(G))}
  summary <- G$summary
  G <- G$G

  #### Unweighted bipartite ####
  if (summary$bipartite == TRUE & summary$weighted == FALSE) {
    if (is.null(s)) {message("The stochastic degree sequence model is suggested. Type \"?sdsm\" for more information.")}
    if (!is.null(s)) {
      message(paste0("Extracting backbone using: sdsm(B, alpha = ", s, ", signed = FALSE, fwer = \"none\", class = \"original\", narrative = TRUE)"))
      backbone <- sdsm(G, alpha = s, fwer = "none", class = summary$class, narrative = TRUE)
      return(backbone)
    }
  }

  #### Weighted bipartite ####
  if (summary$bipartite == TRUE & summary$weighted == TRUE) {

    if (any(G!=as.integer(G)) | any(G < 0)) {message("Backbone models for this type of network are not currently available.")}

    if (all(G==as.integer(G)) & all(G >= 0) & is.null(s)) {message("The ordinal stochastic degree sequence model is suggested. Type \"?osdsm\" for more information.")}

    if (all(G==as.integer(G)) & all(G >= 0) & !is.null(s)) {
      message(paste0("Extracting backbone using: osdsm(B, alpha = ", s, ", trials = 1000, signed = FALSE, fwer = \"none\", class = \"original\", narrative = TRUE)"))
      backbone <- osdsm(G, alpha = s, trials = 1000, signed = FALSE, fwer = "none", class = summary$class, narrative = TRUE)
      return(backbone)
    }
  }

  #### Unweighted unipartite ####
  if (summary$bipartite == FALSE & summary$weighted == FALSE) {
    if (is.null(s)) {message("The L-Spar sparsification model is suggested for revealing subgroups. Type \"?sparsify.with.lspar\" for more information.")}
    if (is.null(s)) {message("The Local Degree sparsification model is suggested for revealing hierarchy. Type \"?sparsify.with.localdegree\" for more information.")}
    if (!is.null(s)) {
      message(paste0("Extracting backbone using: sparsify.with.lspar(G, s = ", s, ", class = \"original\", narrative = TRUE)"))
      backbone <- sparsify.with.lspar(G, s = s, class = summary$class, narrative = TRUE)
      return(backbone)
    }
  }

  #### Weighted unipartite ####
  if (summary$bipartite == FALSE & summary$weighted == TRUE) {

    # Check for possible bipartite projection
    if (all(G%%1==0) &                               #If all entries are integers, and
        any(!(diag(G)%in%c(0,1,NA))) &               #The diagonal is present, and not only 0s and 1s, and
        all((diag(G) == apply(G, 1, FUN=max)))) {    #The diagonal is the largest entry in each row
      message("This object looks like it could be a bipartite projection.")
      message("If so, run backbone.suggest() on the original bipartite network, otherwise...")
    }

    if (is.null(s)) {message("The disparity filter is suggested. Type \"?disparity\" for more information.")}
    if (!is.null(s)) {
      message(paste0("Extracting backbone using: disparity(G, alpha = ", s, ", signed = FALSE, fwer = \"none\", class = \"original\", narrative = TRUE)"))
      backbone <- disparity(G, alpha = s, fwer = "none", class = summary$class, narrative = TRUE)
      return(backbone)
    }
  }

}
