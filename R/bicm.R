#' Computes the loglikelihood gradient for the \link{bicm} function
#'
#' @param x0 vector, probabilities given by current step in bicm function
#' @param args list, c(degree sequence of rows, degree sequence of cols, multiplicity of rows, multiplicity of columns)
#'
#' @return loglikelihood
#' @keywords internal
loglikelihood_prime_bicm <- function(x0, args){
  r_dseq_rows <- args[[1]]
  r_dseq_cols <- args[[2]]
  rows_multiplicity = args[[3]]
  cols_multiplicity = args[[4]]
  num_rows = length(r_dseq_rows)
  num_cols = length(r_dseq_cols)
  x = x0[1:num_rows]
  y = x0[(num_rows+1):length(x0)]
  rm <- rows_multiplicity*x
  cm <- cols_multiplicity*y
  denom <- outer(x,y)+1

  a <- -rowSums(1/sweep(denom, MARGIN = 2, FUN = "/", STATS = cm))
  b <- -colSums(rm/denom)
  a <- a + (r_dseq_rows/x)
  b <- b + (r_dseq_cols/y)
  c <- c(a,b)

  return(c)
}


#' Computes the loglikelihood hessian for the \link{bicm} function
#'
#' @param x0 vector, probabilities given by current step in bicm function
#' @param args list, c(degree sequence of rows, degree sequence of cols, multiplicity of rows, multiplicity of columns)
#'
#' @return hessian matrix
#' @keywords internal
loglikelihood_hessian_diag_bicm <- function(x0, args){
  r_dseq_rows <- args[[1]]
  r_dseq_cols <- args[[2]]
  rows_multiplicity = args[[3]]
  cols_multiplicity = args[[4]]
  num_rows = length(r_dseq_rows)
  num_cols = length(r_dseq_cols)
  x = x0[1:num_rows]
  y = x0[(num_rows+1):length(x0)]
  x2 = x**2
  y2 = y**2
  rm <- rows_multiplicity*x2
  cm <- cols_multiplicity*y2
  denom <- (outer(x,y)+1)**2

  a <- rowSums(1/sweep(denom, MARGIN = 2, FUN = "/", STATS = cm))
  b <- colSums(rm/denom)
  a <- a - (r_dseq_rows/x2)
  b <- b - (r_dseq_cols/y2)
  c <- matrix(c(a,b), 1, (num_rows+num_cols))

  return(c)
}

#' Computes the loglikelihood for the \link{bicm} function
#'
#' @param x0 vector, probabilities given by current step in bicm function
#' @param args list, c(degree sequence of rows, degree sequence of cols, multiplicity of rows, multiplicity of columns)
#'
#' @return loglikelihood, numeric
#' @keywords internal
loglikelihood_bicm <- function(x0, args){
  r_dseq_rows <- args[[1]] #originally 0, have to shift everything by 1
  r_dseq_cols <- args[[2]]
  rows_multiplicity = args[[3]]
  cols_multiplicity = args[[4]]
  num_rows = length(r_dseq_rows)
  num_cols = length(r_dseq_cols)
  x = x0[1:num_rows]
  y = x0[(num_rows+1):length(x0)]
  f = sum(rows_multiplicity*r_dseq_rows*log(x))+sum(cols_multiplicity*r_dseq_cols*log(y))
  f = f-sum((rows_multiplicity%o%cols_multiplicity)*log(x%o%y+1))
  return(f)
}

#' Bipartite Configuration Model
#'
#' `bicm` estimates cell probabilities under the bipartite configuration model
#'
#' @param M matrix: a binary matrix
#' @param fitness boolean: FALSE returns a matrix of probabilities, TRUE returns a list of row and column fitnesses only
#' @param tol numeric, tolerance of algorithm
#' @param max_steps numeric, number of times to run \link{loglikelihood_prime_bicm} algorithm
#' @param ... optional arguments
#'
#' @details
#' Given a binary matrix **M**, the Bipartite Configuration Model (BiCM; Saracco et. al. 2015) returns a valued matrix
#'    **B** in which Bij is the *approximate* probability that Mij = 1 in the space of all binary matrices with
#'    the same row and column marginals as **M**. The BiCM yields the closest approximations of the true probabilities
#'    compared to other estimation methods (Neal et al., 2021), and is used by [sdsm()] to extract the backbone of
#'    a bipartite projection using the stochastic degree sequence model.
#'    
#' Optionally (if `fitness = TRUE`), `bicm()` instead returns a list of row and column fitnesses, which is faster and
#'    requires less memory. Given the *i*th row's fitness Ri and the *j*th column's fitness Rj, the entry Bij in the
#'    matrix can be computed as Ri\*Rj/(1+(Ri\*Rj)).
#'    
#' **Note**: M cannot contain any rows or columns that contain all 0s or all 1s.
#'
#' @references
#' {Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. *Scientific Reports, 5*, 10595. \doi{10.1038/srep10595}}
#'
#' {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @return a matrix of probabilities, or a list of fitnesses
#' @export
#'
#' @examples
#' M <- matrix(c(0,0,1,0,1,0,1,0,1),3,3)  #A binary matrix
#' bicm(M)
bicm <- function(M, fitness = FALSE, tol = 1e-8, max_steps = 200, ...){

  #### initialize_graph ####
  n_rows <- dim(M)[1]
  n_cols <- dim(M)[2]
  rows_deg <- rowSums(M)
  cols_deg <- colSums(M)
  if (any(rows_deg==0) | any(rows_deg==n_cols) | any(cols_deg==0) | any(cols_deg==n_rows)) {stop("M cannot contain any rows/cols that contain all 0s or all 1s")}

  ### find the unique row and column degrees ###
  ### keep track of each cells row/col deg ###
  row_t <- as.data.frame(table(rows_deg))
  r_rows_deg <- as.numeric(as.vector(row_t$rows_deg))     #Unique row degree values
  r_invert_rows_deg <- match((rows_deg), r_rows_deg)      #Rows having each unique degree value
  rows_multiplicity <- as.numeric(as.vector(row_t$Freq))  #Number of rows having each unique degree value
  r_n_rows <- length(r_rows_deg)                          #Number of unique degree values
  col_t <- as.data.frame(table(cols_deg))
  r_cols_deg <- as.numeric(as.vector(col_t$cols_deg))
  r_invert_cols_deg <- match((cols_deg), r_cols_deg)
  cols_multiplicity <- as.numeric(as.vector(col_t$Freq))
  r_n_edges <- sum(rows_deg)

  #### apply the initial guess function ####
  r_x = r_rows_deg/(sqrt(r_n_edges)+1)
  r_y = r_cols_deg/(sqrt(r_n_edges)+1)
  x = c(r_x,r_y)

  #### initialize problem main ####
  args = list(r_rows_deg, r_cols_deg, rows_multiplicity, cols_multiplicity)

  n_steps <- 0
  f_x <- loglikelihood_prime_bicm(x,args)
  norm <- norm(as.matrix(f_x), type = "F")
  diff <- 1

  while ((norm > tol) & (diff > tol) & (n_steps < max_steps)){
    x_old <- x
    B <- as.array(loglikelihood_hessian_diag_bicm(x,args))
    dx <- -f_x/B
    alpha <- 1
    i <- 0
    while ((!(-loglikelihood_bicm(x,args)>-loglikelihood_bicm(x + alpha * dx,args)))&(i<50)){
      alpha <- alpha*.5
      i <- i+1
      }#end while sdc
    x <- x+alpha*dx
    f_x <- loglikelihood_prime_bicm(x,args)
    norm <- norm(as.matrix(f_x), type = "F")
    diff <- norm(as.matrix(x-x_old), type = "F")
    n_steps <- n_steps + 1
    }#end while

  #### set solved problem ####
  sx <- as.vector(rep(0,n_rows))
  sy <- as.vector(rep(0,n_cols))
  row_fitness <- x[1,1:r_n_rows][r_invert_rows_deg]
  col_fitness <- x[1,-(1:r_n_rows)][r_invert_cols_deg]

  #### Return result ####
  if (fitness) {
    fitnesses <- list(rowfit = row_fitness, colfit = col_fitness)
    return(fitnesses)
  } else {
    x <- row_fitness %*% t(col_fitness)
    probs <- x/(x+1)
    return(probs)
  }
  
}



