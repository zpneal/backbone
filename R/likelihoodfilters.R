#' Marginal Likelihood Filter
#'
#' @param G graph: An graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#'
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value,
#'     `negative` is a matrix of probabilities of edge weights being equal to or below the observed value, and
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @references Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms for weighted complex networks. Physical Review E, 93(1), 012304. \doi{10.1103/PhysRevE.93.012304}

#' @export
#'
#' @examples
mlf <- function(G){
  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(G)
  class <- convert$summary[[1]]
  G <- convert$G
  rs <- rowSums(G)
  cs <- colSums(G)

  if (convert$summary[[3]] == TRUE){undirected <- TRUE} #checks for symmetry
  else{undirected <- FALSE}
  if (convert$summary[[4]]==FALSE){stop("Graph must be weighted.")}
  if (convert$summary[[2]]==TRUE){warning("This object is being treated as a unipartite network.")}


  #### Compute Probabilities ####
  ### it is assumed nothing on diagonal ###
  if (undirected){
    ## T = total edges = sum of G/2 ##
    t <- sum(G)/2
    ## m = G[i,j] edge weight ##
    ## p = k_i*k_j/(2*t^2), k_i deg of i, k_j deg of j
    p = outer(rs,rs)/(2*t**2)

    Negative = stats::pbinom(G, t, p, lower.tail = TRUE)
    Positive = stats::pbinom(G-1,t,p, lower.tail = FALSE)
  }
  else if (!undirected){
    ## directed case
    ## T = total edges = sum(G)
    t <- sum(G)
    # m = G[i,j]
    # p = k_i^out*k_j^in/T^2
    # here I'm assuming our directed adj mat is read as G[i,j] = directed i to j
    # k_i^out = rs[i]
    # k_j^in = cs[j]
    p = outer(rs,cs)/(t**2)

    Negative = stats::pbinom(G, t, p, lower.tail = TRUE)
    Positive = stats::pbinom(G-1,t,p, lower.tail = FALSE)
  }

  ### Add back in rownames ###
  rownames(Positive) <- rownames(G)
  colnames(Positive) <- colnames(G)
  rownames(Negative) <- rownames(G)
  colnames(Negative) <- colnames(G)

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- rs
  c <- cs
  a <- c("Input Class", "Model", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Maximum Likelihood Filter", dim(G)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(G)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = Positive, negative = Negative, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}

#' Global Likelihood Filter
#'
#' @param G
#'
#' @return
#' @export
#'
#' @examples
glf <- function(G){
  message("to be completed")
}

