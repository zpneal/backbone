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
  df <- data.frame(row = row(P)[upper.tri(P)],            #Dataframe of dyads in upper triangle
                   col = col(P)[upper.tri(P)], 
                   weight = as.vector(P[upper.tri(P)]))  

  if (missing_as_zero) {df <- df[which(df$weight!=0),]}  #If missing edges should not be tested, remove zero-weight dyads
  
  rs <- rowSums(I)  #Find row sums in bipartite (agent degrees)
  df$row_sum_i <- rs[df$row]
  df$row_sum_j <- rs[df$col]
  df$diff <- ncol(I)-df$row_sum_i  #Difference in total number of artifacts and i's degree

  #### Compute p-values ####
  df$upper <- stats::phyper(df$weight-1, df$row_sum_i, df$diff, df$row_sum_j, lower.tail=FALSE)
  upper <- matrix(NA, nrow = nrow(I), ncol = nrow(I))  #Start with empty matrix of upper-tail p-values
  for (i in 1:nrow(df)) {  #Insert p-values
    row <- df$row[i]
    col <- df$col[i]
    upper[row,col] <- df$upper[i]
    upper[col,row] <- df$upper[i]
  }

  if (signed) {
    df$lower <- stats::phyper(df$weight, df$row_sum_i, df$diff, df$row_sum_j, lower.tail = TRUE)
    lower <- matrix(NA, nrow = nrow(I), ncol = nrow(I))  #Start with empty matrix of upper-tail p-values
    for (i in 1:nrow(df)) {  #Insert p-values
      row <- df$row[i]
      col <- df$col[i]
      lower[row,col] <- df$lower[i]
      lower[col,row] <- df$lower[i]
    }
  }

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
