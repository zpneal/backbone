.onAttach <- function(lib,pkg) {
  local_version <- utils::packageVersion("backbone")
  packageStartupMessage(" ____  ")
  packageStartupMessage("|  _ \\   backbone v",local_version)
  packageStartupMessage("|#|_) |  Cite: Domagalski, R., Neal, Z. P., & Sagan, B. (2021). Backbone: An")
  packageStartupMessage("|# _ <         R package for extracting the backbone of bipartite projections. ")
  packageStartupMessage("|#|_) |        PLoS ONE. https://doi.org/10.1371/journal.pone.0244363")
  packageStartupMessage("|____/   For help: type vignette(\"backbone\"); email zpneal@msu.edu; github domagal9/backbone")
}
