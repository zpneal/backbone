#' Compute universal threshold backbone
#'
#' `universal` returns a unipartite backbone matrix in which
#'     values are set to 1 if above the given upper parameter threshold,
#'     and set to -1 if below the given lower parameter threshold, and are 0 otherwise.
#'
#' @param M Matrix: a weighted adjacency matrix or a bipartite adjacency matrix.
#' @param upper Real or FUN: upper threshold value or function to be applied to the edge weights. Default is 0.
#' @param lower Real or FUN: lower threshold value or function to be applied to the edge weights. Default is NULL.
#' @param bipartite Boolean: TRUE if bipartite matrix, FALSE if weighted matrix. Default is FALSE.
#'
#' @return backbone Matrix: Signed (or positive) adjacency matrix of backbone
#' @export
#'
#' @examples
#' test <- universal(davis%*%t(davis), upper = function(x)mean(x)+sd(x), lower=function(x)mean(x))
#' test2 <- universal(davis, upper = function(x)mean(x)+2*sd(x), lower = 2, bipartite = TRUE)
#' test3 <- universal(davis, upper = 4, lower = 2, bipartite = TRUE)

universal <- function(M,
                      upper = 0,
                      lower = NULL,
                      bipartite = FALSE){
  if ((class(upper)!="function") & (class(upper)!="numeric")) {stop("upper must be either function or numeric")}
  if ((class(lower)!="function") & (class(lower)!="numeric") & (class(lower)!="NULL")) {stop("lower must be either function or numeric")}

  if (bipartite == TRUE){
    P <- M%*%t(M)
  } else {
    P <- M
  }

  #Set threshold values
  if (class(upper) == "function"){
    ut <- upper(P)
  }
  else{ut <- upper}

  if (class(lower) == "function"){
    lt <- lower(P)
  }
  else{lt <- lower}

  #Create backbone matrix
  backbone <- matrix(0, nrow(P), ncol(P))
  negative <- (P<lt)+0
  positive <- (P>ut)+0

  if (length(lower) > 0){
    backbone <- backbone - negative
  }

  backbone <- backbone + positive

  diag(backbone) <- 0
  return(backbone)
}


