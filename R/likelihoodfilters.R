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
#' @param m integer, number of edges to keep
#' @param trials integer, number of matrices to try
#'
#' @return
#' @export
#'
#' @examples
glf <- function(G, m, trials){
  #### CAN DO WAY FASTER IF ONLY RECOMPUTE THE SINGLE EDGE
  #### AKA ONLY COMPUTE THE CHANGE IN LLH

  #### Define Parameters ####
  get_llh <- function(sigma){
    # sigma_ij = weight of ij edge = G

    rs <- rowSums(sigma)
    cs <- colSums(sigma)

    #Nbar = sum_{i<j} sigma_ij
    # sum_{i<j} = lower triangle minus diag
    Nbar = sum(sigma)
    #pij = (k_i *k_j)/(2T^2)
    ## p = k_i*k_j/(2*t^2), k_i deg of i, k_j deg of j
    t <- sum(sigma)/2
    p = outer(rs,cs)/(2*t**2)

    #llh <- log(factorial(Nbar))+sum_{i<j}[sigma_ij log p_ij - log(sigma_ij !)]
    #1st term
    #t1 <- sigma*(log(p))
    #2nd term
    #t2 <- log(factorial(sigma))
    sigma_lower <- sigma
    sigma_lower[upper.tri(sigma_lower)] <- NA
    sigma_lower[diag(sigma_lower)] <- NA

    p_lower <- p
    p_lower[upper.tri(p_lower)] <- NA
    p_lower[diag(p_lower)] <- NA

    inner_sum <- log((p_lower**sigma_lower)/(factorial(sigma_lower)))
    llh <- log(factorial(Nbar))+sum(na.omit(as.numeric(inner_sum)))
    print(llh)
    return(llh)
  } #end get_llh

  delta_llh <- function(old_llh, new_llh){
    prob <- (1-exp(-(new_llh-old_llh)))
  }

  #### Monte Carlo ####

  index <- matrix(1:length(G), nrow = nrow(G), ncol = ncol(G))
  index[upper.tri(index)] <- NA
  index[diag(index)] <- NA

  starting_index <- sample(na.omit(as.numeric(index)), m, replace = FALSE)
  starting_subgraph <- matrix(0, nrow = nrow(G), ncol = ncol(G))
  starting_subgraph[starting_index] <- G[starting_index]

  old_llh <- get_llh(sigma = starting_subgraph)
  current_index <- starting_index
  for (trial in 1:trials){
    print(trial)

    new_index <- index
    new_index[current_index] <- NA
    new_edge <- sample(na.omit(as.numeric(new_index)),1)
    old_edge_index <- sample(length(current_index),1)
    old_edge <- current_index[old_edge_index]
    current_index <- replace(current_index, old_edge_index, new_edge)
    new_subgraph <- matrix(0, nrow=nrow(G), ncol=ncol(G))
    new_subgraph[current_index] <- G[current_index]
    new_llh <- get_llh(sigma = new_subgraph)
    if (new_llh < old_llh){
      print("smaller")
      old_llh <- new_llh
    } else {
      prob <- delta_llh(old_llh,new_llh)
      # 0 back to original
      # 1 keep
      coin_flip <- sample(c(0,1), size = 1, prob = c(prob, 1-prob))
      if (coin_flip == 0){
        current_index <- replace(current_index, old_edge_index, old_edge)
      }#end coin_flip
    }#end else
  }#end trials
  bb <- matrix(0, nrow=nrow(G), ncol=ncol(G))
  bb[current_index] <- G[current_index]
  return(bb)
}#end glf
