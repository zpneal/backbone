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
#' @param graph A graph object of class "matrix", "sparseMatrix", "dataframe", \link[igraph]{igraph}.
#'
#' @return a list(summary, G, attribs)
#'    `summary` is a dataframe containing characteristics of the supplied object
#'    `G` is an adjacency/incidence matrix
#'    `attribs` is a dataframe containing node characteristics if present in `graph`
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(rbinom(5*5,1,.5),5,5)
#' test <- backbone:::tomatrix(M)
tomatrix <- function(graph){
  class <- class(graph)[1]  #Identify class of supplied object
  isbipartite <- FALSE

  if (!(methods::is(graph, "matrix")) & !(methods::is(graph, "Matrix")) & !(methods::is(graph, "igraph")) & !(methods::is(graph, "data.frame"))) {stop("input bipartite data must be a matrix, Matrix, dataframe, or igraph object.")}

  #### Convert from matrix or Matrix object ####
  if (((methods::is(graph, "matrix")) | (methods::is(graph, "Matrix")))) {
    G <- as.matrix(graph)  #Coerce to matrix
    class(G) <- "numeric"  #Coerce to numeric
    if (any(is.na(G))) {stop("The object contains non-numeric entries")}

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

    graph <- igraph::simplify(graph) #Remove any multi-edges and loops

    if (igraph::is.bipartite(graph) == TRUE){ #Bipartite
      isbipartite <- TRUE  #Set type of graph

      #Remove isolates and maximally connected nodes; similar as code below to remove empty/full rows/columns
      isolates <- which(igraph::degree(graph)==0)  #Isolate nodes
      fullrows <- which(igraph::degree(graph, v = igraph::V(graph)$type==FALSE) == sum(igraph::V(graph)$type==TRUE))  #FALSE nodes connected to all TRUE nodes
      fullcols <- which(igraph::degree(graph, v = igraph::V(graph)$type==TRUE) == sum(igraph::V(graph)$type==FALSE))  #TRUE nodes connected to all FALSE nodes
      todrop <- c(isolates, fullrows, fullcols)
      while (length(todrop)>0) {
        graph <- igraph::delete.vertices(graph, todrop)
        isolates <- which(igraph::degree(graph)==0)  #Isolate nodes
        fullrows <- which(igraph::degree(graph, v = igraph::V(graph)$type==FALSE) == sum(igraph::V(graph)$type==TRUE))  #FALSE nodes connected to all TRUE nodes
        fullcols <- which(igraph::degree(graph, v = igraph::V(graph)$type==TRUE) == sum(igraph::V(graph)$type==FALSE))  #TRUE nodes connected to all FALSE nodes
        todrop <- c(isolates, fullrows, fullcols)
      }

      #Capture attributes of FALSE (i.e. row) nodes
      attribs <- as.data.frame(igraph::vertex_attr(graph, index = (igraph::V(graph)$type == FALSE)))  #Get attributes
      attribs <- attribs[,!colnames(attribs)%in%c("type","name"),drop=F]  #Do not need to keep type or name
      attribs <- Filter(function(x)!all(is.na(x)), attribs)  #Remove unused attributes
      if (ncol(attribs)==0) {rm(attribs)}  #Delete dataframe if empty

      #Convert to incidence matrix
      if (is.null(igraph::E(graph)$weight)) {G <- igraph::as_incidence_matrix(graph, sparse = FALSE) #Unweighted bipartite
        } else {G <- igraph::as_incidence_matrix(graph, attr = "weight", sparse = FALSE)} #Weighted bipartite
    }

    if (igraph::is.bipartite(graph) == FALSE){ #Unipartite

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

  #### If the graph is bipartite, remove rows/columns with zero sums (already done above if igraph ####
  if (isbipartite){
    while (any(rowSums(G)==0) | any(rowSums(G)==ncol(G)) | any(colSums(G)==0) | any(colSums(G)==nrow(G))) {  #If any rows/columns are empty/full
      G <- G[which(rowSums(G)!=0 & rowSums(G)!=ncol(G)),which(colSums(G)!=0 & colSums(G)!=nrow(G))]  #Remove them, check again
    }
  }

  #### Summary dataframe ####
  if (isbipartite) {issymmetric <- FALSE} else {issymmetric <- isSymmetric(G)}
  isweighted <- any(!(G%in%c(0,1)))
  summary <- data.frame(
    class = class,
    bipartite = isbipartite,
    symmetric = issymmetric,
    weighted = isweighted)

  if (exists("attribs")) {return(list(summary = summary, G = G, attribs = attribs))} else {return(list(summary = summary, G = G))}
}

#' Converts a backbone adjacency matrix to an object of specified class
#'
#' @param graph a matrix
#' @param attribs dataframe: vertex attributes to be assigned in igraph object
#' @param convert class to convert to, one of "matrix", "Matrix", "igraph", or "edgelist"
#'
#' @return backbone graph: Binary or signed backbone graph of class given in parameter `convert`.
#'
#' @keywords internal
#'
#' @examples
#' M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
#' test <- backbone:::frommatrix(M, "Matrix")
frommatrix <- function(graph, attribs = NA, convert = "matrix"){

  if (convert == "matrix"){G <- graph}

  if (convert == "Matrix"){G <- Matrix::Matrix(graph)}

  if (convert == "edgelist"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
    G <- igraph::as_data_frame(G, what = "edges")
  }

  if (convert == "igraph"){
    if (isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE)}
    if (!isSymmetric(graph)) {G <- igraph::graph.adjacency(graph, mode = "directed", weighted = TRUE)}
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

#' Extracts a backbone network from a backbone object
#'
#' `backbone.extract` returns a binary or signed adjacency matrix
#'      containing the backbone that retains only the significant edges.
#'
#' @param bb.object backbone: backbone S3 class object.
#' @param signed Boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
#' @param alpha Real: significance level of hypothesis test(s)
#' @param mtc string: type of Multiple Test Correction to be applied; can be any method allowed by \code{\link{p.adjust}}.
#' @param class string: the class of the returned backbone graph, one of c("matrix", "sparseMatrix", "igraph", "edgelist"), converted via \link{tomatrix}.
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
      & (class != "edgelist"))
  {stop("incorrect class type, must be one of c(matrix, Matrix, sparseMatrix, igraph, edgelist)")}

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
  backbone <- frommatrix(backbone, convert = class)
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
    #L <- apply(M==1, 1, which, simplify = FALSE)  #Slightly faster, but requires R > 4.1.0
    L <- lapply(asplit(M == 1, 1), which)  #Ensures result is returned as a list
    Lrand <- fastball_cpp(L, trades)
    Mrand <- matrix(0,nrow(M),ncol(M))
    for (row in 1:nrow(Mrand)) {Mrand[row,Lrand[[row]]] <- 1L}
    return(Mrand)
  } else {return(fastball_cpp(M, 5 * length(M)))}
}

#' Poisson binomial distribution function
#'
#' `pb` computes the poisson binomial distribution function using the refined normal approximation.
#'
#' @param k numeric: value where the pdf should be evaluated
#' @param p vector: vector of success probabilities
#' @param lower boolean: If TRUE return both upper & lower tail probabilities,
#'    if FALSE return only upper tail probability
#'
#' @details
#' The Refined Normal Approximation (RNA) offers a close approximation when `length(p)` is
#'    large (Hong, 2013). This function is based on `ppoibin()` from the `poibin` package.
#'
#' @return vector, length 2: The first value (if lower = TRUE) is the lower tail probability, the
#'    probability of observing `k` or fewer successes when each trial has probability `p` of success.
#'    The second value is the upper tail probability, the probability of observing `k` or more
#'    successes when each trial has probability `p` of success.
#'
#' @references
#' {Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. *Computational Statistics and Data Analysis, 59*, 41-51. \doi{10.1016/j.csda.2012.10.006}}
#'
#' @export
#'
#' @examples
#' pb(50,runif(100))
pb <-function(k, p, lower=TRUE) {
  
  #Compute parameters
  mu <- sum(p)
  sigma <- sqrt(sum(p*(1-p)))
  gamma <- sum(p*(1-p)*(1-2*p))
  
  #Lower tail p-value, if requested
  if (lower) {
    x <- (k+.5-mu)/sigma
    lower <- stats::pnorm(x)+gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x)
  } else {lower <- NA}
  
  #Upper tail p-value
  x <- (k-1+.5-mu)/sigma
  upper <- 1 - (stats::pnorm(x)+gamma/(6*sigma^3)*(1-x^2)*stats::dnorm(x))
  
  #Combine, truncate, return
  prob <- c(lower,upper)
  prob[prob<0] <- 0
  prob[prob>1] <- 1
  return(prob)
}
