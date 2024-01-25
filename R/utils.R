.onAttach <- function(lib,pkg) {
  local_version <- utils::packageVersion("backbone")
  packageStartupMessage(" ____   backbone v",local_version)
  packageStartupMessage("|  _ \\  Cite: Neal, Z. P., (2022). Backbone: An R package to extract network backbones.")
  packageStartupMessage("|#|_) |       PLOS ONE, 17, e0269137. https://doi.org/10.1371/journal.pone.0269137")
  packageStartupMessage("|# _ < ")
  packageStartupMessage("|#|_) | Help: type vignette(\"backbone\"); email zpneal@msu.edu; github zpneal/backbone")
  packageStartupMessage("|____/  Beta: type devtools::install_github(\"zpneal/backbone\", ref = \"devel\")")
}

#' Converts an input graph object to an adjacency/incidence matrix and identifies its characteristics
#'
#' @param graph A graph object of class "matrix", "Matrix", "dataframe", \link[igraph]{igraph}.
#'
#' @return a list(summary, G, attribs)
#'    `summary` is a dataframe containing characteristics of the supplied object
#'    `G` is an adjacency/incidence matrix
#'    `attribs` is a dataframe containing vertex attributes (if present in `graph`)
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(rbinom(5*5,1,.5),5,5)
#' test <- backbone:::tomatrix(M)
tomatrix <- function(graph){
  isbipartite <- FALSE
  if (!(methods::is(graph, "matrix")) & !(methods::is(graph, "Matrix")) & !(methods::is(graph, "igraph")) & !(methods::is(graph, "data.frame"))) {stop("input data must be a matrix, Matrix, dataframe, or igraph object.")}

  #### Convert from matrix or Matrix object ####
  if (((methods::is(graph, "matrix")) | (methods::is(graph, "Matrix")))) {
    if (methods::is(graph, "Matrix")) {class <- "Matrix"} else {class <- "matrix"} #Set class
    G <- as.matrix(graph)  #Coerce to matrix
    class(G) <- "numeric"  #Coerce to numeric
    if (any(is.na(G))) {stop("The object contains non-numeric entries")}

    ## Check if bipartite
    if (dim(G)[1]!=dim(G)[2]) {isbipartite <- TRUE}  #A rectangular matrix is treated as bipartite
    if (dim(G)[1]==dim(G)[2] & !is.null(rownames(G)) & !is.null(colnames(G))) { #A labeled square matrix is treated as bipartite if
      if (!identical(rownames(G),colnames(G)) &                                 #the row and column labels differ, and
          !isSymmetric(G)) {                                                    #it is not symmetric
        isbipartite <- TRUE
        }
      }
    }

  #### Convert from edge list ####
  if (methods::is(graph, "data.frame")) {
    if (ncol(graph)<2 | ncol(graph)>3) {stop("An edgelist must contain 2 or 3 columns")}
    class <- "edgelist"  #Update starting class as edgelist
    G <- graph
    colnames(G) <- LETTERS[1:dim(graph)[2]]  #Name columns A, B, and (if present) C
    isbipartite <- length(intersect(G[,1],G[,2])) == 0  #Treat as bipartite if there is no overlap in node lists A and B

      if (isbipartite == TRUE) { #If bipartite, create incidence matrix
        G <- igraph::graph_from_data_frame(G, directed = F)
        igraph::V(G)$type <- igraph::V(G)$name %in% graph[,2] #second column of edges is TRUE type
        if (dim(graph)[2] == 2) {G <- igraph::as_incidence_matrix(G, sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_incidence_matrix(G, attr = "C", sparse = FALSE)} #Weighted
      }

      if (isbipartite == FALSE) { #If unipartite, create adjacency matrix
        G <- igraph::graph_from_data_frame(G, directed = T)
        if (dim(graph)[2] == 2) {G <- igraph::as_adjacency_matrix(G, type = "both", sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_adjacency_matrix(G, type = "both", attr = "C", sparse = FALSE)} #Weighted
      }
    }

  #### Convert from igraph ####
  if (methods::is(graph, "igraph")) {
    class <- "igraph"

    if (("weight" %in% igraph::edge_attr_names(graph)) & igraph::any_multiple(graph)) {stop("A weighted igraph cannot contain multi-edges.")}
    if (!("weight" %in% igraph::edge_attr_names(graph)) & igraph::any_multiple(graph)) {igraph::E(graph)$weight <- 1}  #If unweighted with multi-edges, give each edge a weight of 1
    graph <- igraph::simplify(graph, edge.attr.comb = "sum") #Remove any loops, and sum multi-edges into a weight attribute

    ## For bipartite inputs
    if (igraph::is.bipartite(graph) == TRUE){
      isbipartite <- TRUE  #Set type of graph

      #Capture attributes of FALSE (i.e. row) nodes
      attribs <- as.data.frame(igraph::vertex_attr(graph, index = (igraph::V(graph)$type == FALSE)))  #Get attributes
      attribs <- attribs[,!colnames(attribs)%in%c("type","name"),drop=F]  #Do not need to keep type or name
      attribs <- Filter(function(x)!all(is.na(x)), attribs)  #Remove unused attributes
      if (ncol(attribs)==0) {rm(attribs)}  #Delete dataframe if empty

      #Convert to incidence matrix
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_incidence_matrix(graph, sparse = FALSE) #Unweighted bipartite
        } else {G <- igraph::as_incidence_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted bipartite
    }

    ## For unipartite inputs
    if (igraph::is.bipartite(graph) == FALSE){

      #Capture attributes of nodes
      attribs <- as.data.frame(igraph::vertex_attr(graph))  #Get attributes
      attribs <- attribs[,!colnames(attribs)%in%c("type","name"),drop=F]  #Do not need to keep type or name
      attribs <- Filter(function(x)!all(is.na(x)), attribs)  #Remove unused attributes
      if (ncol(attribs)==0) {rm(attribs)}  #Delete dataframe if empty

      #Convert to adjacency matrix
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_adjacency_matrix(graph, sparse = FALSE) #Unweighted unipartite
      } else {G <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted unipartite
    }
  }

  #### Summary dataframe ####
  if (isbipartite) {issymmetric <- FALSE} else {issymmetric <- isSymmetric(G)}  #Check if symmetric
  isweighted <- any(!(G%in%c(0,1)))  #Check if weighted
  summary <- data.frame(class = class, bipartite = isbipartite, symmetric = issymmetric, weighted = isweighted)  #Description

  #### Return result ####
  if (exists("attribs")) {return(list(summary = summary, G = G, attribs = attribs))} else {return(list(summary = summary, G = G))}
}

#' Converts a backbone adjacency matrix to a graph object of specified class
#'
#' @param mat an adjacency matrix
#' @param attribs dataframe: vertex attributes to be assigned in igraph object
#' @param convert class to convert to, one of "matrix", "Matrix", "igraph", or "edgelist"
#'
#' @return backbone graph: Binary or signed backbone graph of class `convert`.
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
#' test <- backbone:::frommatrix(M, "Matrix")
frommatrix <- function(mat, attribs = NA, convert = "matrix"){

  if (convert == "matrix"){G <- mat}

  if (convert == "Matrix"){G <- Matrix::Matrix(mat)}

  if (convert == "edgelist"){
    if (isSymmetric(mat)) {G <- igraph::graph.adjacency(mat, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(mat)) {G <- igraph::graph.adjacency(mat, mode = "directed", weighted = TRUE)}
    G <- igraph::as_data_frame(G, what = "edges")
  }

  if (convert == "igraph"){
    if (isSymmetric(mat)) {G <- igraph::graph.adjacency(mat, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(mat)) {G <- igraph::graph.adjacency(mat, mode = "directed", weighted = TRUE)}
    if (igraph::gsize(G)!=0) {igraph::E(G)$sign <- igraph::E(G)$weight}  #To facilitate use with library(signnet)

    if (!is.null(attribs)) {  #If attributes are supplied, add them to the igraph object
      attribs.to.add <- colnames(attribs)
      for (attrib in attribs.to.add) {
        G <- igraph::set_vertex_attr(graph = G, name = attrib, value = attribs[,attrib])
      }
    }

  }

  return(G)
}

#' Estimate number of monte carlo trials needed to estimate p-value
#'
#' `trials.needed` estimates the number of monte carlo trials needed to estimate edgewise p-values
#'
#' @param M matrix: An adjacency matrix representing a network from which a backbone is is being extracted
#' @param alpha real: significance level of hypothesis test(s)
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param how.close real: How close can the empirical p-value be to alpha and still be distinguishable, expressed as a proportion
#'
#' @return integer
#' @keywords internal
trials.needed <- function(M, alpha, signed, missing.as.zero, mtc, how.close = 0.95) {
  #In signed network, empirical p-value is compared to alpha/2
  if (signed == TRUE) {test.alpha <- alpha / 2} else {test.alpha <- alpha}

  #If a multiple test correction will be performed, conservatively adjust alpha using Bonferroni (i.e., based on number of tests performed)
  if (mtc != "none") {
    if (isSymmetric(M) & !missing.as.zero) {tests <- sum(lower.tri(M) & M!=0)}  #Non-zero entries in lower triangle
    if (isSymmetric(M) & missing.as.zero) {tests <- sum(lower.tri(M))}  #Entries in lower triangle
    if (!isSymmetric(M) & !missing.as.zero) {tests <- sum(lower.tri(M) & M!=0) + sum(upper.tri(M) & M!=0)}  #Non-zero entries in lower and upper triangles
    if (!isSymmetric(M) & missing.as.zero) {tests <- nrow(M)^2 - nrow(M)}  #Entries in upper and lower triangle
    test.alpha <- test.alpha / tests
  }

  #p1 = A close-to-alpha hypothetical empirical monte carlo p-value we want to evaluate
  #p2 = The alpha level against which we are evaluating p1, with any one-tailed or mtc adjustments
  #Because type-I errors (a false edge is included in the backbone) is as bad as type-II errors (a true edge is omitted from the backbone), therefore power = alpha
  trials <- ceiling((stats::power.prop.test(p1 = test.alpha * how.close, p2 = test.alpha, sig.level = alpha, power = (1-alpha), alternative = "one.sided")$n)/2)
  return(trials)
}
