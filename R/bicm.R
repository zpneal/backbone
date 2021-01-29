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

  f = integer(length(x0))
  flag = TRUE

  for (i in 1:num_rows){
    for (j in 1:num_cols){
      denom = 1+x[i]*y[j]
      f[i] <- f[i] - ((y[j]*cols_multiplicity[j])/denom)
      f[(j+num_rows)] = f[(j+num_rows)] - ((x[i]*rows_multiplicity[i])/denom)
      if (flag == TRUE){
        f[(j+num_rows)] = f[(j+num_rows)] + (r_dseq_cols[j]/y[j])
      }#end flag
    }#end for j
    f[i] = f[i]+(r_dseq_rows[i]/x[i])
    flag <- FALSE
  }#end for i
  return(f)
}


#' Computes the loglikelihood hessian for the \link{bicm} function
#'
#' @param x0 vector, probabilities given by current step in bicm function
#' @param args list, c(degree sequence of rows, degree sequence of cols, multiplicity of rows, multiplicity of columns)
#'
#' @return hessian matrix
#' @keywords internal
loglikelihood_hessian_diag_bicm <- function(x0, args){
  r_dseq_rows <- args[[1]] #originally 0, have to shift everything by 1
  r_dseq_cols <- args[[2]]
  rows_multiplicity = args[[3]]
  cols_multiplicity = args[[4]]
  num_rows = length(r_dseq_rows)
  num_cols = length(r_dseq_cols)
  x = x0[1:num_rows]
  y = x0[(num_rows+1):length(x0)]

  f = matrix(0, 1, (num_rows+num_cols))
  x2 = x**2
  y2 = y**2
  flag <- TRUE

  for (i in 1:num_rows){
    for (j in 1:num_cols){
      denom = (1+x[i]*y[j])**2
      f[i] = f[i] + (cols_multiplicity[j]*y2[j])/denom
      f[(j+num_rows)] = f[(j+num_rows)] + ((rows_multiplicity[i]*x2[i])/denom)
      if (flag == TRUE){
        f[j+num_rows] = f[j+num_rows] - (r_dseq_cols[j]/y2[j])
      }
    }
    f[i] = f[i] - (r_dseq_rows[i]/x2[i])
    flag <- FALSE
  }
  return(f)
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
  flag <- TRUE

  f = 0
  for (i in 1:num_rows){
    f = f+ (rows_multiplicity[i]*r_dseq_rows[i]*log(x[i]))
    for (j in 1:num_cols){
      if (flag == TRUE){
        f = f + (cols_multiplicity[j]*r_dseq_cols[j]*log(y[j]))
      }#end if flag
      f = f - (rows_multiplicity[i]*cols_multiplicity[j]*log((1+x[i]*y[j])))
    }#end for j
    flag <- FALSE
  }#end for i
  return(f)
}

#' bicm: Bipartite Configuration Model.
#'
#' @param graph matrix, a bipartite adjacency matrix of a graph
#' @param tol numeric, tolerance of algorithm
#' @param max_steps numeric, number of times to run \link{loglikelihood_prime_bicm} algorithm
#' @param progress Boolean: If \link[utils]{txtProgressBar} should be used to measure progress
#'
#' @details The Bipartite Configuration Model (Saracco et. al. 2015, 2017) produces a matrix of edge specific probabilities which are used in \link{sdsm} to find the p-values of the edges in the bipartite projection. This R code is adapted from the python BiCM package by Matteo Bruno under the MIT license.
#' @references python bicm: \href{https://github.com/mat701/BiCM}{Matteo Bruno, matteo.bruno<at>imtlucca.it, https://github.com/mat701/BiCM}
#' @references bicm: \href{https://doi.org/10.1088/1367-2630/aa6b38}{Saracco, F., Straka, M. J., Clemente, R. D., Gabrielli, A., Caldarelli, G., & Squartini, T. (2017). Inferring monopartite projections of bipartite networks: An entropy-based approach. New Journal of Physics, 19(5), 053022. DOI: 10.1088/1367-2630/aa6b38}
#' @references bicm: \href{https://doi.org/10.1038/srep10595}{Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. Scientific Reports, 5(1), 10595. DOI: 10.1038/srep10595}
#'
#' @return matrix containing probabilities
#' @export
#'
#' @examples bicm(davis)
bicm <- function(graph, tol = 1e-8, max_steps = 200, progress = FALSE){

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
      alpha <- alpha*.05
      i <- i+1
    }#end while sdc

    x <- x+alpha*dx
    f_x <- loglikelihood_prime_bicm(x,args)
    norm <- norm(as.matrix(f_x), type = "F")
    diff <- norm(as.matrix(x-x_old), type = "F")
    n_steps <- n_steps + 1
    if (n_steps == 1){
      if (progress == TRUE){
      pb <- utils::txtProgressBar(min = -log(abs(max(diff,norm)-tol)),
                                  max = tol,
                                  style = 3)
      }
    }
    if (progress==TRUE){
      utils::setTxtProgressBar(pb, -log(abs(max(diff,norm)-tol)))}
  }#end while
  if (progress == "TRUE"){close(pb)}

  #### set solved problem ####
  sx <- as.vector(rep(0,n_rows))
  sy <- as.vector(rep(0,n_cols))
  sr_xy <- x
  sr_x <- sr_xy[1:r_n_rows]
  sr_y <- sr_xy[-(1:r_n_rows)]
  sx[nonfixed_rows] <- sr_x[r_invert_rows_deg]
  sy[nonfixed_cols] <- sr_y[r_invert_cols_deg]

  #### plug everything into probability matrix ####
  r_probs <- matrix(0, length(sx[nonfixed_rows]), length(sy[nonfixed_cols]))
  for (i in 1:length(sx[nonfixed_rows])){
    for (j in 1:length(sy[nonfixed_cols])){
      xy <- sx[nonfixed_rows][i]*sy[nonfixed_cols][j]
      r_probs[i,j] <- xy/(1+xy)
    }#end j
  }#end i

  probs[nonfixed_rows,nonfixed_cols] <- r_probs
  return(probs)
}



