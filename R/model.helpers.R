#' Convert graph object to adjacency matrix
#'
#' @param graph, matrix, sparse matrix, igraph, edgelist, or network object
#' @param convert, class to convert to, one of "matrix", "sparseMatrix", "igraph", "edgelist", or "network"
#'
#' @return list(class, adjacency), a list containing the class of parameter graph, and the adjacency matrix of the graph
#' @export
#' @examples
#' davis.sp <- as(davis, "sparseMatrix")
#' davis.graph <- igraph::graph.incidence(davis)
#' davis.nw <- network::network(davis, ignore.eval = FALSE, names.eval = "weight", loops = TRUE)
#' class.convert(davis, "matrix")
#' class.convert(davis.sp, "matrix")
#' class.convert(davis.graph, "matrix")
#' class.convert(davis.nw, "matrix")
#' bb.sdsm <- sdsm(davis)
#' bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)
#' class.convert(bb, "matrix")
#' class.convert(bb, "sparseMatrix")
#' class.convert(bb, "igraph")
#' class.convert(bb, "network")
class.convert <- function(graph, convert = "matrix"){
  class <- class(graph)
  if (convert == "matrix"){
    if ((methods::is(graph, "matrix")) | (methods::is(graph, "sparseMatrix"))) {
      if (dim(graph)[2] == 2){ # for edgelists, assumes bipartite!
        g <- igraph::graph.data.frame(graph, directed = F)
        igraph::V(g)$type <- igraph::V(g)$name %in% graph[,2] #second column of edges is TRUE type
        G <- igraph::get.incidence(g)
        class <- "edgelist"}
      else{G <- graph}
    }
    if (methods::is(graph, "igraph")) {
      if (igraph::is.bipartite(graph)){
        G <- igraph::get.incidence(graph)
      }
      else {
        if (length(igraph::edge.attributes(graph)) > 0){G <- igraph::get.adjacency(graph,attr = "weight")}
        else G <- igraph::get.adjacency(graph)
      }
    }
    if (methods::is(graph, "network")) {G <- as.matrix(graph, attr = "weight")}
  }
  if (convert == "sparseMatrix"){
    G <- methods::as(graph, "sparseMatrix")
  }
  if (convert == "igraph"){
    if (-1 %in% graph){G <- igraph::graph.adjacency(graph,mode = "undirected", weighted = TRUE)}
    else G <- igraph::graph.adjacency(graph, mode = "undirected")
  }
  if (convert == "network"){
    G <- network::network(graph, ignore.eval = FALSE, names.eval = "weight", directed = FALSE)
  }
  if (convert == "edgelist"){
    g <- igraph::graph.adjacency(graph, weighted = TRUE)
    G <- igraph::get.data.frame(g)
  }
  return(list(class, G))
}

#' Polytope Method for
#'
#' @param G matrix, an adjacency matrix representing a graph
#'
#' @return matrix of probabilities
#' @export
#'
#' @examples
#' polytope(davis)
polytope <- function(G){

  #Define Variable to solve for
  matrix <- CVXR::Variable(dim(G)[1], dim(G)[2])

  #For row & column sums
  mat1 <- matrix(1, dim(G)[2], 1)
  mat2 <- matrix(1, dim(G)[1], 1)

  #Define constraints
  constraint1 <- matrix >= 0
  constraint2 <- matrix <= 1
  if (methods::is(G, "sparseMatrix")) {constraint3 <- (matrix%*%mat1) == Matrix::rowSums(G)}
  else {constraint3 <- (matrix%*%mat1) == rowSums(G)}
  if (methods::is(G, "sparseMatrix")) {constraint4 <- t(matrix)%*%mat2 == Matrix::colSums(G)}
  else {constraint4 <- t(matrix)%*%mat2 == colSums(G)}
  constraints <- list(constraint1, constraint2, constraint3, constraint4)

  #Define Objective, the function to solve
  objective <- CVXR::Maximize(sum(CVXR::entr(matrix)+CVXR::entr(1-matrix)))

  #Define problem: objective with the constrants
  problem <- CVXR::Problem(objective, constraints)

  #Solve the problem
  result <- CVXR::psolve(problem, warm_start = TRUE)

  #Results
  new_matrix <- result$getValue(matrix)

  #Matrix of Probabilities
  return(new_matrix)
}


#' curveball algorithm
#'
#' @param m, matrix
#'
#' @return rm, random matrix with same row sums and column sums as m
#' @export
#'
#' @references Strona, G., Nappo, D., Boccacci, F., Fattorini, S., San-Miguel-Ayanz, J. (2014). A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals. Nature Communications, 5, 4114
#' @examples
#' curveball(davis)
curveball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

#' Poisson Binomial distribution computed with Refined Normal Approximation
#'
#' @param kk, values where the cdf is to be computed
#' @param pp, vector of success probabilities for indicators
#' @param wts, the weights for each probability
#'
#' @return cdf
#' @export
#'
#' @references Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. Computational Statistics & Data Analysis, Vol. 59, pp. 41-51.
#' @details These values are approximated using the Refined Normal Approximation (RNA method).
#'     These functions are originally described by \link[poibin]{ppoibin} and used here under GPL-2 License.
#'
#' @examples
#' probs <- polytope(davis)
#' P <- davis %*% t(davis)
#' prob.mat <- matrix(probs, nrow = nrow(davis), ncol = ncol(davis))
#' prob.imat <- sweep(prob.mat, MARGIN = 2, prob.mat[1,], `*`)
#' mapply(rna, kk= as.data.frame(t(P[1,])), pp = as.data.frame(t(prob.imat)))
rna <-function(kk,pp,wts=NULL)
{
  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }
  if(is.null(wts))
  {
    wts=rep(1,length(pp))
  }
  pp=rep(pp,wts)
  muk=sum(pp)
  sigmak=sqrt(sum(pp*(1-pp)))
  gammak=sum(pp*(1-pp)*(1-2*pp))
  ind=gammak/(6*sigmak^3)
  kk1=(kk+.5-muk)/sigmak
  vkk.r=stats::pnorm(kk1)+gammak/(6*sigmak^3)*(1-kk1^2)*stats::dnorm(kk1)
  vkk.r[vkk.r<0]=0
  vkk.r[vkk.r>1]=1
  res=vkk.r
  return(res)
}

#' Family-wise Error Rates: Holm-Bonferroni method
#'
#' @param backbone, backbone class object
#' @param alpha, numeric, an alpha value for significance testing
#' @param signed, Boolean, if a signed backbone should be returned instead of a binary backbone
#'
#' @return backbone, a binary or signed matrix
#' @export
#'
#' @examples
#' probs <- sdsm(davis)
#' fwer(probs, alpha = .4, signed = FALSE)
fwer <- function(backbone, alpha = 0.05, signed = TRUE){
  #Read in data
  matrix_positive <- backbone$positive
  matrix_negative <- backbone$negative

  #For the loop, not needed for a single run
  order <- dim(matrix_positive)[1]

  #Matrix of smallest pvalues
  if (signed == TRUE){
    pvalue_matrix <- pmin(matrix_positive, matrix_negative)
    #True if positive value wins, false if negative value wins
    sign <- matrix_positive < matrix_negative
  }
  else{
    #pvalue all positive p values
    pvalue_matrix <- matrix_positive
    #sign matrix is all trues
    sign <- matrix_positive
    sign[sign>-1]<-TRUE
  }

  #Create df to hold p values and keep track of row and col of edges
  df <- data.frame(as.vector(pvalue_matrix))
  df$row <- rep(1:nrow(pvalue_matrix), times=ncol(pvalue_matrix))
  df$col <- rep(1:ncol(pvalue_matrix), each=nrow(pvalue_matrix))

  #Separate sorted matrix, just to be extra careful
  #Sort based on p-value, smallest to largest
  df_sorted <- df[order(df$as.vector.pvalue_matrix.),]
  #Add sign
  df_sorted$sign <- as.vector(sign)
  #Add empty col for the edge weights after FWER
  df_sorted$newvalues <- NA

  #Compute the rank
  independent_analysis <- ((order[[1]])*(order[[1]]-1))/2
  #For iterations
  j = 0
  #Alpha value
  if (signed == TRUE){
    FWER <- (alpha/2)
  }
  else{
    FWER = alpha
  }

  #Run over all edges
  for (i in 1:dim(df_sorted)[1]){
    #If row index less than col index (so we only have to run upper triangle)
    if (df_sorted$row[i]<df_sorted$col[i]){
      #If value less than the fwer fraction, reject the null hyp
      if (df_sorted$as.vector.pvalue_matrix.[i] < ((FWER)/(independent_analysis-j))){
        #The new edge weight should be 1 if pos smaller, -1 if neg smaller
        df_sorted$newvalues[i] <- 2*sign[i]-1
        print(paste0(j, ": Edge (", df_sorted$row[i], ",", df_sorted$col[i], ") created with weight:", df_sorted$newvalues[i]," p-value was:", df_sorted$as.vector.pvalue_matrix.[i]))
        #Increase iteration
        j <- j+1
      } #end if pvalue < fwer
      else{
        #If value NOT less than fwer fraction, fail to reject everything further, save current values
        u <- df_sorted$row[i]
        v <- df_sorted$col[i]
        weight <- 2*sign[i]-1
        pval <- df_sorted$as.vector.pvalue_matrix.[i]
        val <- j-1
        break
      }
    } #end if row<col
  } #end for i in 1:dim[1]

  print(paste0("Added ", val, " edges to network"))
  print(paste0("Failed to reject H0 for pair (",u, ",", v, ") created with weight: ",pval, " and threshold ",((FWER/2)/(independent_analysis-i))))

  #Construct the backbone matrix
  new_mat_values <- df_sorted[order(df_sorted$row, df_sorted$col),]
  backbone_ut <- t(matrix(new_mat_values$newvalues, nrow = nrow(matrix_positive), ncol = ncol(matrix_positive)))+0
  backbone_lt <- matrix(new_mat_values$newvalues, nrow = nrow(matrix_positive), ncol = ncol(matrix_positive))
  backbone_ut[is.na(backbone_ut)]<-0
  backbone_lt[is.na(backbone_lt)]<-0
  backbone <- backbone_ut+backbone_lt
  row.names(backbone) <- colnames(matrix_positive)
  colnames(backbone) <- colnames(matrix_positive)
  return(backbone)
}
