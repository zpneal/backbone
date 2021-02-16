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

#' bicm: Bipartite Configuration Model.
#'
#' @param graph matrix, a bipartite adjacency matrix of a graph
#' @param tol numeric, tolerance of algorithm
#' @param max_steps numeric, number of times to run \link{loglikelihood_prime_bicm} algorithm
#' @param progress Boolean: If \link[utils]{txtProgressBar} should be used to measure progress
#' @param ... optional arguments
#'
#' @details The Bipartite Configuration Model (Saracco et. al. 2015, 2017) produces a matrix of edge specific probabilities which are used in \link{sdsm} to find the p-values of the edges in the bipartite projection. This R code is adapted from the python BiCM package by Matteo Bruno under the MIT license.
#' @references python bicm: \href{https://github.com/mat701/BiCM}{Matteo Bruno, matteo.bruno<at>imtlucca.it, https://github.com/mat701/BiCM}
#' @references bicm: {Saracco, F., Straka, M. J., Clemente, R. D., Gabrielli, A., Caldarelli, G., & Squartini, T. (2017). Inferring monopartite projections of bipartite networks: An entropy-based approach. New Journal of Physics, 19(5), 053022. \doi{10.1088/1367-2630/aa6b38}}
#' @references bicm: {Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. Scientific Reports, 5(1), 10595. \doi{10.1038/srep10595}}
#'
#' @return matrix containing probabilities
#' @export
#'
#' @examples bicm(davis)
bicm <- function(graph,
                 tol = 1e-8,
                 max_steps = 200,
                 progress = FALSE,
                 ...){

  #### initialize_graph ####
  n_rows <- dim(graph)[1]
  n_cols <- dim(graph)[2]
  rows_deg <- Matrix::rowSums(graph)
  cols_deg <- Matrix::colSums(graph)

  #### initialize probability matrix ####
  probs <- matrix(0, nrow = n_rows, ncol = n_cols)
  r_bipart <- graph
  good_rows <- seq(1,n_rows)
  good_cols <- seq(1,n_cols)
  zero_rows <- which(rows_deg == 0)
  zero_cols <- which(cols_deg == 0)
  full_rows <- which(rows_deg == n_cols)
  full_cols <- which(cols_deg == n_rows)
  num_full_rows <- 0
  num_full_cols <- 0

  while ((length(zero_rows)
          + length(zero_cols)
          + length(full_rows)
          + length(full_cols)) > 0){
    if (length(zero_rows)>0){
      r_bipart <- r_bipart[-zero_rows,]
      good_rows <- good_rows[-zero_rows]
    }
    if (length(zero_cols)>0){
      r_bipart <- r_bipart[, -zero_cols]
      good_cols <- good_cols[-zero_cols]
    }
    full_rows = which(Matrix::rowSums(r_bipart) == dim(r_bipart)[2])
    full_cols = which(Matrix::colSums(r_bipart) == dim(r_bipart)[1])
    num_full_rows = num_full_rows + length(full_rows)
    num_full_cols = num_full_cols + length(full_cols)
    probs[good_rows[full_rows],good_cols] <- 1
    probs[good_rows,good_cols[full_cols]] <- 1
    if (length(full_rows)>0){
      good_rows <- good_rows[-full_rows]
      r_bipart <- r_bipart[-full_rows,]
    }
    if (length(full_cols)>0){
      good_cols <- good_cols[-full_cols]
      r_bipart <- r_bipart[,-full_cols]
    }
    zero_rows <- which(Matrix::rowSums(r_bipart) == 0)
    zero_cols <- which(Matrix::colSums(r_bipart) == 0)
  } #end while
  nonfixed_rows = good_rows
  fixed_rows = seq(1,n_rows)[-good_rows]
  nonfixed_cols = good_cols
  fixed_cols = seq(1,n_cols)[-good_cols]

  #### degree reduction ####
  if (length(nonfixed_rows)>0){
    rows_deg <- rows_deg[nonfixed_rows]
  }
  if (length(nonfixed_cols)>0){
    cols_deg <- cols_deg[nonfixed_cols]
  }

  ### find the unique row and column degrees ###
  ### keep track of each cells row/col deg ###
  row_t <- as.data.frame(table(rows_deg-num_full_cols))
  r_rows_deg <- as.numeric(as.vector(row_t$Var1))
  r_invert_rows_deg <- match((rows_deg-num_full_cols), r_rows_deg)
  rows_multiplicity <- as.numeric(as.vector(row_t$Freq))
  r_n_rows <- length(r_rows_deg)
  col_t <- as.data.frame(table(cols_deg - num_full_rows))
  r_cols_deg <- as.numeric(as.vector(col_t$Var1))
  r_invert_cols_deg <- match((cols_deg-num_full_rows), r_cols_deg)
  cols_multiplicity <- as.numeric(as.vector(col_t$Freq))
  r_n_edges <- sum(rows_deg-num_full_cols)

  #### apply the initial guess function ####
  ### an initial probability to start with ###
  r_x = as.vector(r_rows_deg)/(sqrt(r_n_edges)+1)
  r_y = r_cols_deg/(sqrt(r_n_edges)+1)
  ### save all together as a vector ###
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
    while ((!(-loglikelihood_bicm(x,args)>-loglikelihood_bicm(x + alpha * dx,args)))&
           (i<50)){
      alpha <- alpha*.5
      i <- i+1
    }#end while sdc

    x <- x+alpha*dx
    f_x <- loglikelihood_prime_bicm(x,args)
    norm <- norm(as.matrix(f_x), type = "F")
    diff <- norm(as.matrix(x-x_old), type = "F")

    if ((progress == TRUE)&(n_steps == 1)){
      mx <- max(diff, norm)
      pb <- utils::txtProgressBar(min = n_steps,
                                  max = max_steps,
                                  style = 3)
    }
    if ((progress==TRUE)&(n_steps>1)){
      utils::setTxtProgressBar(pb, n_steps)
    }
    n_steps <- n_steps + 1
  }#end while
  if (progress == "TRUE"){
    utils::setTxtProgressBar(pb,max_steps)
    close(pb)
  }

  #### set solved problem ####
  sx <- as.vector(rep(0,n_rows))
  sy <- as.vector(rep(0,n_cols))
  sr_xy <- x
  sr_x <- sr_xy[1:r_n_rows]
  sr_y <- sr_xy[-(1:r_n_rows)]
  sx[nonfixed_rows] <- sr_x[r_invert_rows_deg]
  sy[nonfixed_cols] <- sr_y[r_invert_cols_deg]

  #### plug everything into probability matrix ####
  f <- function(x,y){
    z <- x*y
    z/(1+z)
  }
  r_probs <- outer(sx[nonfixed_rows],sy[nonfixed_cols],f)

  probs[nonfixed_rows,nonfixed_cols] <- r_probs
  return(probs)
}



