#' Logit-based probabilities for SDSM
#'
#' `logit` estimates cell probabilities under the logit model
#'
#' @param M matrix
#'
#' @details
#' Given a matrix **M**, the logit model returns a valued matrix **B** in which Bij is the *approximate* probability
#'    that Mij = 1 in the space of all binary matrices with the same row and column marginals as **M**.
#'
#'    The Bipartite Configuration Model (BiCM), which is available using \link{bicm} is faster and yields slightly
#'    more accurate probabilities (Neal et al., 2021). Therefore, it is the default used in \link{sdsm}. However,
#'    the BiCM it requires the assumption that any cell in **M** can take a value of 0 or 1.
#'
#'    In contrast, the logit model allows constraints on specific cells. If **M** represents a bipartite graph, these
#'    constraints are equivalent to structural 0s (an edge that can never be present) and structural 1s (an edge that
#'    must always be present). To impose such constraints, **M** should be binary, except that structural 0s are
#'    represented with Mij = 10, and structural 1s are represented with Mij = 11.
#'
#' @references package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. *PLOS ONE, 17*, e0269137. \doi{10.1371/journal.pone.0269137}}
#' @references logit model: {Neal, Z. P. (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. *Social Networks, 39*, 84-97. \doi{10.1016/j.socnet.2014.06.001}}
#' @references probability comparison: {Neal, Z. P., Domagalski, R., and Sagan, B. (2021). Comparing alternatives to the fixed degree sequence model for extracting the backbone of bipartite projections. *Scientific Reports, 11*, 23929. \doi{10.1038/s41598-021-03238-3}}
#'
#' @return a matrix of probabilities
#' @export
#'
#' @examples
#' M <- matrix(c(0,0,1,0,1,0,1,0,1),3,3)  #A binary matrix
#' logit(M)
#' M <- matrix(c(0,10,1,0,1,0,1,0,11),3,3)  #A binary matrix with structural values
#' logit(M)
logit <- function(M) {

  #Vectorize the bipartite data
  A <- data.frame(as.vector(M))
  names(A)[names(A)=="as.vector.M."] <- "value"

  #Assign row and column IDs in the vectorized data
  A$row <- rep(1:nrow(M), times=ncol(M))
  A$col <- rep(1:ncol(M), each=nrow(M))

  #Set structural values to zero so they're not included in rowsum & columnsum
  A$value2 <- A$value
  A$value2[which(A$value>1)] <- 0

  #Compute and attach rowsums & columnsums, not counting structural values
  A$rowmarg <- stats::ave(A$value2,A$row,FUN=sum)
  A$colmarg <- stats::ave(A$value2,A$col,FUN=sum)

  #Set structural values to missing so they're not included in the logit
  A$value2[which(A$value>1)] <- NA

  #Compute probabilities on non-structural dyads using logit
  model.estimates <- suppressWarnings(stats::glm(formula = value2 ~ rowmarg + colmarg, family = stats::binomial(link="logit"), data=A))
  A$probs <- as.vector(suppressWarnings(stats::predict(model.estimates, newdata = A, type = "response")))

  #Insert structural probabilities
  A$probs[which(A$value==10)] <- 0  #Structural zeros have probability = 0
  A$probs[which(A$value==11)] <- 1  #Structural ones have probability = 1

  #Probability matrix
  prob.mat <- matrix(A$probs, nrow = nrow(M), ncol = ncol(M))  #Probability matrix
  return(prob.mat)
}
