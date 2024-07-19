#' Compute edgewise p-values under the Stochastic Degree Sequence Model
#'
#' @param I a binary incidence matrix
#' @param missing.as.zero boolean: should missing edges be treated as edges with zero weight and tested for significance
#' @param signed boolean: TRUE for a signed backbone, FALSE for a binary backbone (see details)
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
.sdsm <- function(I, missing.as.zero, signed){
  
  P <- tcrossprod(I)  #Weighted bipartite projection

  probs <- bicm(I)  #Bipartite edge probabilities under BiCM
  probs <- lapply(seq_len(nrow(probs)), function(i) probs[i,])  #Store probabilities as list
  
  #### Compute p-values ####
  Pupper <- matrix(NA, nrow(P), ncol(P), dimnames = list(rownames(P),colnames(P)))                #Set upper-tail p-value to NA initially
  if (signed) {Plower <- matrix(NA, nrow(P), ncol(P), dimnames = list(rownames(P),colnames(P)))}  #If signed, set lower-tail p-value to NA initially
    
  for (col in 1:ncol(P)) {  #Loop over lower triangle of projection
    for (row in col:nrow(P)) {

      if (missing.as.zero) {  #If missing edges should be treated as zero, test each one
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}
        
        if (signed) {Plower[row,col] <- pvalues[1]}
        Pupper[row,col] <- pvalues[2]
      }

      if (!missing.as.zero & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}
        
        if (signed) {Plower[row,col] <- pvalues[1]}
        Pupper[row,col] <- pvalues[2]
      }

    }
  }
  Pupper[upper.tri(Pupper)] <- t(Pupper)[upper.tri(Pupper)]  #Make upper-tail p-value matrix symmetric
  if (signed) {Plower[upper.tri(Plower)] <- t(Plower)[upper.tri(Plower)]}  #Make lower-tail p-value matrix symmetric

  if (signed) {return(list(Plower = Plower, Pupper = Pupper))}
  if (!signed) {return(list(Pupper = Pupper))}
  }