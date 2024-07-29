#' Compute edgewise p-values under the Stochastic Degree Sequence Model with Edge Constraints
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
#' @references sdsm-ec model: {Neal, Z. P. and Neal, J. W. (2023). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for Backbone Extraction. *International Conference on Complex Networks and Their Applications, 12*, 127-136. \doi{10.1007/978-3-031-53468-3_11}}
#'
#' @noRd
.sdsm_ec <- function(I, missing_as_zero, signed){

  #### Construct weighted projection ####
  I_unweighted <- I
  I_unweighted[I_unweighted==10] <- 0  #Make structural 0s ordinary 0
  I_unweighted[I_unweighted==11] <- 1  #Make structural 1s ordinary 1
  P <- tcrossprod(I_unweighted)  #Projection, not considering any structural 0s or 1s
  
  #### Compute probabilities with edge constraints using Logit ####
  # Prepare dyad list
  A <- data.frame(edge = as.vector(I),     #Data frame of bipartite dyads
                  row = as.vector(row(I)),
                  col = as.vector(col(I)))
  A$edge2 <- A$edge  #Set structural edges to NA so they're excluded from degrees and logit
  A$rowmarg <- stats::ave(A$edge2,A$row,FUN=sum)  #Compute rowsums (agent degree), excluding structural edges
  A$colmarg <- stats::ave(A$edge2,A$col,FUN=sum)  #Compute colsums (artifact degree), excluding structural edges

  #Compute probabilities on non-structural edges using logit
  model.estimates <- suppressWarnings(stats::glm(formula = edge2 ~ rowmarg + colmarg, family = stats::binomial(link="logit"), data=A))
  A$probs <- as.vector(suppressWarnings(stats::predict(model.estimates, newdata = A, type = "response")))
  
  #Insert structural probabilities
  A$probs[which(A$edge==10)] <- 0  #Structural zeros have probability = 0
  A$probs[which(A$edge==11)] <- 1  #Structural ones have probability = 1
  
  #Probability matrix
  probs <- matrix(A$probs, nrow = nrow(I), ncol = ncol(I))  #Probability matrix
  probs <- lapply(seq_len(nrow(probs)), function(i) probs[i,])  #Store probabilities as list

  #### Compute p-values ####
  upper <- matrix(NA, nrow(P), ncol(P))                #Set upper-tail p-value to NA initially, untested edges have p = NA
  if (signed) {lower <- matrix(NA, nrow(P), ncol(P))}  #If signed, set lower-tail p-value to NA initially

  for (col in 1:(ncol(P)-1)) {  #Loop over lower triangle of projection
    for (row in (col+1):nrow(P)) {

      if (missing_as_zero) {  #If missing edges should be treated as zero, test each one
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

      if (!missing_as_zero & P[row,col] != 0) {  #If missing edges should not be treated as zero, test only edges with non-zero weight
        if (!signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = FALSE)}
        if (signed) {pvalues <- .pb(P[row,col], unlist(Map('*',probs[row],probs[col])), lowertail = TRUE)}

        if (signed) {lower[row,col] <- pvalues[1]}
        upper[row,col] <- pvalues[2]
      }

    }
  }
  upper[upper.tri(upper)] <- t(upper)[upper.tri(upper)]  #Make upper-tail p-value matrix symmetric
  if (signed) {lower[upper.tri(lower)] <- t(lower)[upper.tri(lower)]}  #Make lower-tail p-value matrix symmetric

  if (signed) {return(list(lower = lower, upper = upper))}
  if (!signed) {return(list(upper = upper))}
  }
