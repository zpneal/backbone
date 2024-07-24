#' Depricated function to extract SDSM backbone
#' See backbone v2.1.4 for original documentation
#' @export
#' @noRd
sdsm <- function(B, alpha = 0.05, missing.as.zero = FALSE, signed = FALSE, mtc = "none", narrative = FALSE, ...){
  message("The sdsm() function is depricated in backbone v3.")
  message("This command should mostly work for now, however use")
  message("backbone_from_bipartite(model = \"sdsm\") in the future.")

  return(
  backbone_from_bipartite(B,
                          model = "sdsm",
                          alpha = alpha,
                          signed = signed,
                          mtc = mtc,
                          missing_as_zero = missing.as.zero,
                          narrative = TRUE)
  )
}
