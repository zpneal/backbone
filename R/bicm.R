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
#' Matrix **M** is "conforming" if no rows and no columns contain only zeros or only ones. If **M** is conforming, then
#'   `bicm()` is faster. Additionally, if `fitness = TRUE`, then `bicm()` returns a list of row and column fitnesses,
#'    which requires less memory. Given the *i*th row's fitness Ri and the *j*th column's fitness Rj, the entry Bij in
#'    the probability matrix can be computed as Ri x Rj/(1+(Ri x Rj)).
#'
#' Matrix **M** is "non-conforming" if any rows or any columns contain only zeros or only ones. If **M** is non-conforming,
#'    then `bicm()` is slower and will only return a probability matrix.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references bicm: {Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015). Randomizing bipartite networks: The case of the World Trade Web. *Scientific Reports, 5*, 10595. \doi{10.1038/srep10595}}
#'
#' @return a matrix of probabilities or a list of fitnesses
#' @export
#'
#' @examples
#' M <- matrix(c(0,0,1,0,1,0,1,0,1),3,3)  #A binary matrix
#' bicm(M)
bicm <- function(M, fitness = FALSE, tol = 1e-8, max_steps = 200, ...){

  #### Compute parameters ####
  n_rows <- dim(M)[1]
  n_cols <- dim(M)[2]
  rows_deg <- rowSums(M)
  cols_deg <- colSums(M)
  nonconforming <- any(rows_deg==0) | any(rows_deg==n_cols) | any(cols_deg==0) | any(cols_deg==n_rows)

  #### BiCM for conforming matrices ####
  if (!nonconforming) {

    ## find the unique row and column degrees, and their multiplicity
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

    ## apply initial guess function
    r_x = r_rows_deg/(sqrt(r_n_edges)+1)
    r_y = r_cols_deg/(sqrt(r_n_edges)+1)
    x = c(r_x,r_y)

    ## initialize problem
    args = list(r_rows_deg, r_cols_deg, rows_multiplicity, cols_multiplicity)
    n_steps <- 0
    f_x <- loglikelihood_prime_bicm(x,args)
    norm <- norm(as.matrix(f_x), type = "F")
    diff <- 1

    ## solve problem
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

    ## format solution
    sx <- as.vector(rep(0,n_rows))
    sy <- as.vector(rep(0,n_cols))
    row_fitness <- x[1,1:r_n_rows][r_invert_rows_deg]
    col_fitness <- x[1,-(1:r_n_rows)][r_invert_cols_deg]

    ## return probability matrix or fitnesses
    if (fitness) {
      fitnesses <- list(rowfit = row_fitness, colfit = col_fitness)
      return(fitnesses)
    } else {
      x <- row_fitness %*% t(col_fitness)
      probs <- x/(x+1)
      return(probs)
    }
  }

  #### BiCM for conforming matrices ####
  if (nonconforming) {
    ## initialize probability matrix
    probs <- matrix(0, nrow = n_rows, ncol = n_cols)
    r_bipart <- M
    good_rows <- seq(1,n_rows)
    good_cols <- seq(1,n_cols)
    zero_rows <- which(rows_deg == 0)
    zero_cols <- which(cols_deg == 0)
    full_rows <- which(rows_deg == n_cols)
    full_cols <- which(cols_deg == n_rows)
    num_full_rows <- 0
    num_full_cols <- 0

    ## identify non-conforming rows and columns
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

    ## re-compute row and column degrees
    if (length(nonfixed_rows)>0){
      rows_deg <- rows_deg[nonfixed_rows]
    }
    if (length(nonfixed_cols)>0){
      cols_deg <- cols_deg[nonfixed_cols]
    }

    ## find the unique row and column degrees, and their multiplicity
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

    ## apply initial guess function
    r_x = r_rows_deg/(sqrt(r_n_edges)+1)
    r_y = r_cols_deg/(sqrt(r_n_edges)+1)
    x = c(r_x,r_y)

    ## initialize problem
    args = list(r_rows_deg, r_cols_deg, rows_multiplicity, cols_multiplicity)
    n_steps <- 0
    f_x <- loglikelihood_prime_bicm(x,args)
    norm <- norm(as.matrix(f_x), type = "F")
    diff <- 1

    ## solve problem
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
      n_steps <- n_steps + 1
    }#end while

    ## format solution
    sx <- as.vector(rep(0,n_rows))
    sy <- as.vector(rep(0,n_cols))
    sr_xy <- x
    sr_x <- sr_xy[1:r_n_rows]
    sr_y <- sr_xy[-(1:r_n_rows)]
    sx[nonfixed_rows] <- sr_x[r_invert_rows_deg]
    sy[nonfixed_cols] <- sr_y[r_invert_cols_deg]

    ## return probability matrix
    f <- function(x,y){
      z <- x*y
      z/(1+z)
    }
    r_probs <- outer(sx[nonfixed_rows],sy[nonfixed_cols],f)
    probs[nonfixed_rows,nonfixed_cols] <- r_probs
    return(probs)
  }

}



