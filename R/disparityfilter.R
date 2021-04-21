#' Disparity Filter Backbone Probabilities
#'
#' @param G graph: Graph object of class matrix, sparse matrix, igraph, edgelist, or network object.
#'     Any rows and columns of the associated bipartite matrix that contain only zeros are automatically removed before computations.
#'
#' @return backbone, a list(positive, negative, summary). Here
#'     `positive` is a matrix of probabilities of edge weights being equal to or above the observed value in the projection,
#'     `summary` is a data frame summary of the inputted matrix and the model used including: model name, number of rows, skew of row sums, number of columns, skew of column sums, and running time.
#' @references Serrano, M. Á., Boguñá, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483–6488. \doi{10.1073/pnas.0808904106}
#' @export
#'
#' @examples
#' disparityfilter(davis%*%t(davis))
disparityfilter <- function(G){
  ### Run Time ###
  run.time.start <- Sys.time()

  #### Class Conversion ####
  convert <- tomatrix(G)
  G <- convert$G
  class <- convert$summary[[1]]
  symmetric <- convert$summary[["symmetric"]]
  strength <- rowSums(G)
  binary <- (G>0)+0
  degree <- rowSums(binary) #equivalent to k

  #Find all weight 0 edges
  zeros <- G==0

  if (symmetric == TRUE){
    P <- G/strength
    pvalues <- (1-P)^(degree-1)
    positive <- as.matrix(pvalues)
    #negative <- 1-positive
  }

  if (symmetric == FALSE){
    outp <- G/strength
    outvalues <- (1-outp)^(degree-1)
    inp <- t(G)/(colSums(G))
    invalues <- t((1-inp)^(colSums(binary)-1))
    positive <- pmin(invalues,outvalues)
    #negative <- 1-positive
  }

  ### If edge weight was zero
  ##     set to 1 in positive so not in backbone
  positive[zeros] <- 1

  ### Run Time ###
  run.time.end <- Sys.time()
  total.time = (round(difftime(run.time.end, run.time.start, units = "secs"), 2))

  #### Compile Summary ####
  r <- rowSums(G)
  c <- colSums(G)

  a <- c("Input Class", "Method", "Number of Rows", "Mean of Row Sums", "SD of Row Sums", "Skew of Row Sums", "Number of Columns", "Mean of Column Sums", "SD of Column Sums", "Skew of Column Sums", "Running Time (secs)")
  b <- c(class[1], "Disparity Filter", dim(G)[1], round(mean(r),5), round(stats::sd(r),5), round((sum((r-mean(r))**3))/((length(r))*((stats::sd(r))**3)), 5), dim(G)[2], round(mean(c),5), round(stats::sd(c),5), round((sum((c-mean(c))**3))/((length(c))*((stats::sd(c))**3)), 5), as.numeric(total.time))
  model.summary <- data.frame(a,b, row.names = 1)
  colnames(model.summary)<-"Model Summary"

  #### Return Backbone Object ####
  bb <- list(positive = positive, summary = model.summary)
  class(bb) <- "backbone"
  return(bb)
}
