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
#' @references fixedrow: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @noRd
.fixedrow <- function(I, missing_as_zero, signed){

  P <- tcrossprod(I)  #Weighted bipartite projection

  #### Prepare dyad list ####
  df <- data.frame(row = row(P)[upper.tri(P)],            #Dataframe of dyads in upper triangle
                   col = col(P)[upper.tri(P)], 
                   weight = as.vector(P[upper.tri(P)]))  
  
  rs <- rowSums(I)  #Find row sums in bipartite (agent degrees)
  df$row_sum_i <- rs[df$row]
  df$row_sum_j <- rs[df$col]
  df$diff <- ncol(I)-df$row_sum_i  #Difference in total number of artifacts and i's degree

  if (!missing_as_zero) {df$weight[which(df$weight==0)] <- NA}  #If missing edges should not be tested, replace weight with NA
  
  #### Compute p-values ####
  df$upper <- stats::phyper(df$weight-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)
  upper <- matrix(NA, nrow = nrow(P), ncol = nrow(P))  #Start with empty matrix of upper-tail p-values
  upper[upper.tri(upper)] <- df$upper  #Insert upper-tail p-values
  upper[lower.tri(upper)] = t(upper)[lower.tri(upper)]  #Make symmetric

  if (signed) {
    df$lower <- stats::phyper(df$weight, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)
    lower <- matrix(NA, nrow = nrow(P), ncol = nrow(P))  #Start with empty matrix of upper-tail p-values
    lower[upper.tri(lower)] <- df$lower
    lower[lower.tri(lower)] = t(lower)[lower.tri(lower)]
  }

  #### If missing edges should *not* be treated as having zero weight, remove p-value and do not consider for backbone ####
  if (!missing.as.zero) {
    upper[P == 0] <- NA
    if (signed) {lower[P == 0] <- NA}
  }
  
  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
