.onAttach <- function(lib,pkg) {
  local_version <- utils::packageVersion("backbone")
  packageStartupMessage(" ____   backbone v",local_version)
  packageStartupMessage("|  _ \\  Cite: Domagalski, R., Neal, Z. P., & Sagan, B. (2021). Backbone: An")
  packageStartupMessage("|#|_) |       R package for extracting the backbone of bipartite projections. ")
  packageStartupMessage("|# _ <        PLoS ONE, 16(1), e0244363. https://doi.org/10.1371/journal.pone.0244363")
  packageStartupMessage("|#|_) | Help: type vignette(\"backbone\"); email zpneal@msu.edu; github zpneal/backbone")
  packageStartupMessage("|____/  Beta: type devtools::install_github(\"zpneal/backbone\", ref = \"devel\")")
}

#' Converts an input graph object to an adjacency/incidence matrix and identifies its characteristics
#'
#' @param graph A graph object of class "matrix", "sparseMatrix", "dataframe", \link[igraph]{igraph}, \link[network]{network}.
#'
#' @return a list(summary, G)
#'    `summary` is a dataframe containing characteristics of the supplied object
#'    `G` is an adjacency/incidence matrix
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(rbinom(5*5,1,.5),5,5)
#' test <- backbone:::tomatrix(M)
tomatrix <- function(graph){
  class <- class(graph)[1]
  isbipartite <- FALSE

  if (!(methods::is(graph, "matrix")) & !(methods::is(graph, "sparseMatrix")) & !(methods::is(graph, "Matrix")) & !(methods::is(graph, "igraph")) & !(methods::is(graph, "network")) & !(methods::is(graph, "data.frame"))) {stop("input bipartite data must be a matrix, edgelist dataframe, igraph, or network object.")}

  #### Convert from matrix object ####
  if (((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix")) | (methods::is(graph, "Matrix")))) {
    G <- as.matrix(graph)  #Coerce to matrix
    class(G) <- "numeric"  #Coerce to numeric
    if (any(is.na(G))) {stop("The object contains non-numeric entries")}

    if (dim(G)[1]!=dim(G)[2]) {isbipartite <- TRUE}  #A rectangular matrix is treated as bipartite
    if (dim(G)[1]==dim(G)[2] & !is.null(rownames(G)) & !is.null(colnames(G))) { #A labeled square matrix is treated as bipartite IFF
      if (!identical(rownames(G),colnames(G)) &                                 #the row and column labels differ, and
          !isSymmetric(G)) {                                                    #it is not symmetric
        isbipartite <- TRUE
        }
      }
    }

  #### Convert from edge list ####
  if (methods::is(graph, "data.frame")) {
    if (ncol(graph)==1 | ncol(graph)>3) {stop("An edgelist must contain 2 or 3 columns")}
    class <- "edgelist"  #Update starting class as edgelist
    G <- graph
    colnames(G) <- LETTERS[1:dim(graph)[2]]  #Name columns A, B, and (if present) C
    isbipartite <- length(intersect(G[,1],G[,2])) == 0  #Treat as bipartite if there is no overlap in node lists A and B

      if (isbipartite == TRUE) { #Bipartite
        G <- igraph::graph_from_data_frame(G, directed = F)
        igraph::V(G)$type <- igraph::V(G)$name %in% graph[,2] #second column of edges is TRUE type
        if (dim(graph)[2] == 2) {G <- igraph::as_incidence_matrix(G, sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_incidence_matrix(G, attr = "C", sparse = FALSE)} #Weighted
      }

      if (isbipartite == FALSE) { #Unipartite
        G <- igraph::graph_from_data_frame(G, directed = T)
        if (dim(graph)[2] == 2) {G <- igraph::as_adjacency_matrix(G, type = "both", sparse = FALSE)} #Unweighted
        if (dim(graph)[2] == 3) {G <- igraph::as_adjacency_matrix(G, type = "both", attr = "C", sparse = FALSE)} #Weighted
      }
    }

  #### Convert from igraph ####
  if (methods::is(graph, "igraph")) {

    if (igraph::is.bipartite(graph) == TRUE){ #Bipartite
      isbipartite <- TRUE
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_incidence_matrix(graph, sparse = FALSE) #Unweighted
      } else {G <- igraph::as_incidence_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted
    }

    if (igraph::is.bipartite(graph) == FALSE){ #Unipartite
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_adjacency_matrix(graph, sparse = FALSE) #Unweighted
      } else {G <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted
    }
  }

  #### Convert from statnet ####
  if (methods::is(graph, "network")) {

    if (network::is.bipartite(graph) == TRUE) { #Bipartite
      isbipartite <- TRUE
      if ("weight" %in% network::list.edge.attributes(graph)) {
        G <- network::as.matrix.network(graph, type = "incidence", attrname = "weight")} else { #Weighted
          G <- network::as.matrix.network(graph, type = "incidence")} #Unweighted
    }

    if (network::is.bipartite(graph) == FALSE) { #Unipartite
      if ("weight" %in% network::list.edge.attributes(graph)) {
        G <- network::as.matrix.network(graph, type = "adjacency", attrname = "weight")} else { #Weighted
          G <- network::as.matrix.network(graph, type = "adjacency")} #Unweighted
    }
  }

  #### If the graph is bipartite, remove rows/columns with zero sums ####
  if (isbipartite){
    R <- Matrix::rowSums(G)
    C <- Matrix::colSums(G)
    r <- which(R == 0)
    c <- which(C == 0)
    if (length(r)>0){G <- G[-r,]}
    if (length(c)>0){G <- G[,-c]}
  }

  #### Summary dataframe ####
  if (isbipartite) {issymmetric <- FALSE} else {issymmetric <- isSymmetric(G)}
  isweighted <- any(!(G%in%c(0,1)))
  summary <- data.frame(
    class = class,
    bipartite = isbipartite,
    symmetric = issymmetric,
    weighted = isweighted)

  #### Report input type and modifications ####
  if (issymmetric) {dir <- "undirected"} else {dir <- "directed"}
  if (isweighted) {weigh <- "a weighted"} else {weigh <- "an unweighted"}
  if (isbipartite) {
    message(paste0("This ", class, " object is being treated as ", weigh, " bipartite network of ", nrow(G), " agents and ", ncol(G), " artifacts."))
    if (length(r)>0) {message("These empty (i.e. all 0s) rows have been removed from the data: ", paste0(r, " "))}
    if (length(c)>0) {message("These empty (i.e. all 0s) columns have been removed from the data: ", paste0(c, " "))}
  }
  if (!isbipartite) {message(paste0("This ", class, " object is being treated as ", weigh, " ", dir, " network containing ", nrow(G), " nodes."))}

  return(list(summary = summary, G = G))
}

#' Converts a backbone adjacency matrix to an object of specified class
#'
#' @param graph a matrix
#' @param convert class to convert to, one of "matrix", "sparseMatrix", "igraph", "edgelist", or "network"
#'
#' @return backbone graph: Binary or signed backbone graph of class given in parameter `convert`.
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
#' test <- backbone:::frommatrix(M, "Matrix")
frommatrix <- function(graph, convert = "matrix"){

  if (convert == "matrix"){G <- graph}

  if (convert == "Matrix"){G <- Matrix::Matrix(graph)}

  if (convert == "sparseMatrix"){G <- methods::as(graph, "sparseMatrix")}

  if (convert == "network"){
    if (isSymmetric(graph)) {G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = FALSE)}
    if (!isSymmetric(graph)) {G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = TRUE)}
  }

  if (convert == "edgelist"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
    G <- igraph::as_data_frame(G, what = "edges")
  }

  if (convert == "igraph"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
    if (igraph::gsize(G)!=0) {igraph::E(G)$sign <- igraph::E(G)$weight}  #To facilitate use with library(signnet)
  }

  return(G)
}

#' Extracts a backbone network from a backbone object
#'
#' `backbone.extract` returns a binary or signed adjacency matrix
#'      containing the backbone that retains only the significant edges.
#'
#' @param bb.object backbone: backbone S3 class object.
#' @param signed Boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param alpha Real: significance level of hypothesis test(s)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("matrix", "sparseMatrix", "igraph", "network", "edgelist"), converted via \link{tomatrix}.
#' @return backbone graph: Binary or signed backbone graph of class given in parameter `class`.
#'
#' @details The "backbone" S3 class object is composed of three matrices (the weighted graph, edges' upper-tail p-values,
#'    edges' lower-tail p-values), and a string indicating the null model used to compute p-values.
#'
#' When `signed = FALSE`, a one-tailed test (is the weight stronger) is performed for each edge with a non-zero weight. It
#'    yields a backbone that perserves edges whose weights are significantly *stronger* than expected in the chosen null
#'    model. When `signed = TRUE`, a two-tailed test (is the weight stronger or weaker) is performed for each every pair of nodes.
#'    It yields a backbone that contains positive edges for edges whose weights are significantly *stronger*, and
#'    negative edges for edges whose weights are significantly *weaker*, than expected in the chosen null model.
#'    *NOTE: Before v2.0.0, all significance tests were two-tailed and zero-weight edges were evaluated.*
#'
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
#' backbone.object <- fixedrow(B, alpha = NULL)
#' bb <- backbone.extract(backbone.object, alpha = 0.05)
backbone.extract <- function(bb.object, signed = FALSE, alpha = 0.05, mtc = "none", class = "matrix"){

  #### Argument Checks ####
  if ((alpha >= 1) | (alpha <= 0)) {stop("alpha must be between 0 and 1")}
  if ((class != "matrix")
      & (class != "Matrix")
      & (class != "sparseMatrix")
      & (class != "igraph")
      & (class != "network")
      & (class != "edgelist"))
  {stop("incorrect class type, must be one of c(matrix, Matrix, sparseMatrix, igraph, network, edgelist)")}

  #### Extract object components ####
  G <- bb.object$G
  Pupper <- bb.object$Pupper
  Plower <- bb.object$Plower
  model <- bb.object$model

  #### Extract signed backbone (two-tailed test; all dyads considered) ####
  if (signed) {
    alpha <- alpha / 2  #Use two-tailed test
    Psmaller <- pmin(Pupper,Plower)  #Find smaller p-value
    diag(Psmaller) <- NA
    Ptail <- (Pupper < Plower)  #Find tail of smaller p-value (TRUE if smaller p-value is in upper tail)
    diag(Ptail) <- NA

    if (mtc != "none") {  #Adjust p-values for familywise error, if requested
      if (isSymmetric(Psmaller)) {Psmaller[upper.tri(Psmaller)] <- NA}  #If undirected, ignore upper triangle
      p <- as.vector(Psmaller)  #Vector of p-values
      m <- sum((!is.na(p))*1)  #Number of p-values to evaluate, number of independent edges to test
      p <- stats::p.adjust(p, method = mtc, m)  #Adjust p-values
      Psmaller <- matrix(p, nrow = nrow(Psmaller), ncol = ncol(Psmaller))  #Put adjusted p-values in original p-value matrix
      if (all(is.na(Psmaller[upper.tri(Psmaller)]))) {Psmaller[upper.tri(Psmaller)] <- t(Psmaller)[upper.tri(Psmaller)]}  #If upper triangle is missing, put it back
    }

    backbone <- (Psmaller < alpha)*1  #Identify all significant edges
    backbone[which(Ptail==FALSE)] <- backbone[which(Ptail==FALSE)] * -1  #Make lower-tail significant edges negative
    backbone[which(is.na(backbone))] <- 0  #fill NAs with 0s
    rownames(backbone) <- rownames(G)
    colnames(backbone) <- rownames(G)
  }

  #### Extract binary backbone (one-tailed test; only positively-weighed edges considered) ####
  if (!signed) {
    Pupper[which(G==0)] <- NA  #Eliminate p-values for zero-weight edges; not relevant
    diag(Pupper) <- NA  #Eliminate p-values for loops; not relevant

    if (mtc != "none") {  #Adjust p-values for familywise error, if requested
      if (isSymmetric(Pupper)) {Pupper[upper.tri(Pupper)] <- NA}  #If undirected, ignore upper triangle
      p <- as.vector(Pupper)  #Vector of p-values
      m <- sum((!is.na(p))*1)  #Number of p-values to evaluate, number of independent edges to test
      p <- stats::p.adjust(p, method = mtc, m)  #Adjust p-values
      Pupper <- matrix(p, nrow = nrow(Pupper), ncol = ncol(Pupper))  #Put adjusted p-values in original p-value matrix
      if (all(is.na(Pupper[upper.tri(Pupper)]))) {Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]}  #If upper triangle is missing, put it back
    }

    backbone <- (Pupper < alpha)*1  #Identify all significant edges
    backbone[which(is.na(backbone))] <- 0  #fill NAs with 0s
    rownames(backbone) <- rownames(G)
    colnames(backbone) <- rownames(G)
  }

  #### Return result ####
  backbone <- frommatrix(backbone, class)
  return(backbone)
}

#' Randomize a binary matrix using the fastball algorithm
#'
#' `fastball` randomizes a binary matrix, preserving the row and column sums
#'
#' @param M matrix: a binary matrix (see details)
#' @param trades integer: number of trades; the default is 5R trades (approx. mixing time)
#'
#' @return
#' matrix: A random binary matrix with same row sums and column sums as M.
#'
#' @details
#' Given a matrix `M`, `fastball` randomly samples a new matrix from the space of all matrices with the same row and column sums as `M`.
#'
#' @references {Godard, Karl and Neal, Zachary P. 2022. fastball: A fast algorithm to sample bipartite graphs with fixed degree sequences. \href{https://arxiv.org/abs/2112.04017}{*arXiv:2112.04017*}}
#'
#' @export
#' @examples
#' M <- matrix(rbinom(200,1,0.5),10,20)  #A random 10x20 binary matrix
#' Mrand <- fastball(M)  #Random matrix with same row and column sums
fastball <- function(M, trades = 5 * nrow(M)) {
  if (methods::is(M, "matrix")) {
    L <- apply(M==1, 1, which, simplify = FALSE)
    Lrand <- fastball_cpp(L, trades)
    Mrand <- matrix(0,nrow(M),ncol(M))
    for (row in 1:nrow(Mrand)) {Mrand[row,Lrand[[row]]] <- 1L}
    return(Mrand)
  } else {return(fastball_cpp(M, 5 * length(M)))}
}
