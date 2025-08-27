#' Normalize edge scores
#'
#' @param W a weighted adjacency matrix
#' @param normalize string: type of normalization
#'
#' @return A weighted adjacency matrix
#'
#' @noRd
.normalize <- function(W, normalize) {

  #### Neighborhood rank, from Satuluri et al. (2011) ####
  if (normalize == "rank" | normalize == "embeddedness") {
    for (i in 1:nrow(W)) {  #For each row (i.e., from the perspective of each node)
      x <- W[i,]  #Vector of values from this row
      old <- sort(unique(x))  #Find unique values
      new <- c((length(old)):1)  #Rank them 1 = highest, 2 = second highest, etc
      if (min(old)==0) {new[which(new==max(new))] <- 0}  #If zero was one of the values, rank them as 0
      x <- new[match(x, old)]  #Replace original values with corresponding ranks
      W[i,] <- x  #Put ranks into row
      }
    }

  #### Embeddedness, from Nick et al. (2013) ####
  if (normalize == "embeddedness") {  #Scores will already be transformed as neighborhood ranks
    scores <- matrix(0, nrow(W), ncol(W))  #Initialize matrix to hold embeddedness scores
    for (row1 in 1:(nrow(W)-1)) {
      for (row2 in (row1+1):nrow(W)) {  #Loop over each pair of rows
        list1 <- W[row1,-c(row1,row2)]  #Vector of ranked edges for row1, excluding row1 and row2
        list2 <- W[row2,-c(row1,row2)]  #Vector of ranked edges for row2, excluding row1 and row2

        #Find overlap between neighborhoods using non-parametric variant
        k <- max(list1,list2)
        if (k==0 | ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k))))==0) {  #If jaccard for max(k) is zero, stop
          scores[row1,row2] <- 0
        } else {  #Otherwise, compute jaccard for each k, use maximum
          j <- NULL
          for (k in 1:max(list1,list2)) {j <- c(j, ((sum((list1>0 & list1<=k) & (list2>0 & list2<=k))) / (sum((list1>0 & list1<=k) | (list2>0 & list2<=k)))))}
          scores[row1,row2] <- max(j)
        }
      }
    }
    W <- scores * ((W!=0)*1)
    W[lower.tri(W)] <- t(W)[lower.tri(W)]  #Fill in rest of matrix
  }

  return(W)
}
