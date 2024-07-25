#' Compute edgewise p-values under the Fixed Row Model
#'
#' @param I a binary incidence matrix
#' @param missing_as_zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone
#'
#' @return
#' If `signed = FALSE` a list containing a matrix of upper-tail p-values.
#'
#' If `signed = TRUE` a list containing a matrix of lower-tail and upper-tail p-values
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references sdsm: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance, and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#'
#' @noRd
.fixedrow <- function(I, missing_as_zero, signed){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Prepare dyad list ####
  #Although P is symmetric, this is a complete list of dyads; AFTER TESTING, TRY TO RE-WRITE WITH ONLY NECESSARY EDGES
  df <- data.frame(as.vector(P))  #Dataframe of dyads
  names(df)[names(df)=="as.vector.P."] <- "weight"

  rs <- rowSums(I)  #Find row sums in bipartite (agent degrees)
  df$row_sum_i <- rep(rs, times = nrow(I))  #Add rowsums to dataframe
  df$row_sum_j <- rep(rs, each = nrow(I))
  df$diff <- ncol(I)-df$row_sum_i  #Difference in total number of artifacts and i's degree

  #### Compute p-values ####
  df$upper <- stats::phyper(df$weight-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)
  upper <- matrix(as.numeric(df$upper), nrow = nrow(I), ncol = nrow(I))

  if (signed) {
    df$lower <- stats::phyper(df$weight, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)
    lower <- matrix(as.numeric(df$lower), nrow = nrow(I), ncol = nrow(I))
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing_as_zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
