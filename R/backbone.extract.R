#' Extracts the backbone of a weighted network using results from a null model
#'
#' `backbone.extract` returns a binary or signed adjacency matrix
#'      containing the backbone that retains only the significant edges.
#'
#' @param matrix list: list containing two matrices, positive and negative, and a summary list as returned by \link{sdsm}, \link{fdsm}, or \link{hyperg}.
#' @param signed Boolean: TRUE if signed backbone is to be returned, FALSE if binary backbone is to be returned
#' @param alpha Real: Precision of significance test (one-tailed if only the positive matrix supplied, two-tailed if positive and negative matrices supplied)
#'
#' @return backbone Matrix: Binary or signed adjacency matrix of backbone graph.
#' @export
#'
#' @examples
#' probs <- sdsm(davis)
#' bb <- backbone.extract(probs, alpha = .1)
backbone.extract <- function(matrix, signed = TRUE, alpha = 0.05){

  #Argument Checks
  if ((alpha >= 1) | (alpha <= 0)) {stop("alpha must be between 0 and 1")}
  #if (length(negative) != 0) {alpha <- alpha / 2}  #Use a two-tailed test for signed backbones

  positive <- matrix[[1]]
  negative <- matrix[[2]]
  summary <- matrix$summary
  class <- as.character(summary[1,1])

  #Convert values to matrix
  SignedPositive <- as.matrix((positive<=alpha)+0)
  SignedNegative <- as.matrix((negative<=alpha)+0)
  SignedNegative[SignedNegative==1] <- -1

  #Create backbone matrix
  if (signed == "FALSE") {backbone <- SignedPositive
  } else {backbone <- SignedPositive + SignedNegative}
  diag(backbone) <- 0
  rownames(backbone) <- rownames(positive)
  colnames(backbone) <- rownames(positive)

  if ((class == "dgCMatrix") | (class == "dgRMatrix")){
    class <- "sparseMatrix"
  }
  backbone <- class.convert(backbone, class)
  return(backbone[[2]])
}
