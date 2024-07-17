.onAttach <- function(lib,pkg) {
  local_version <- utils::packageVersion("backbone")
  packageStartupMessage(" ____   backbone v",local_version)
  packageStartupMessage("|  _ \\  Cite: Neal, Z. P., (2022). Backbone: An R package to extract network backbones.")
  packageStartupMessage("|#|_) |       PLOS ONE, 17, e0269137. https://doi.org/10.1371/journal.pone.0269137")
  packageStartupMessage("|# _ < ")
  packageStartupMessage("|#|_) | Help: type vignette(\"backbone\"); email zpneal@msu.edu; github zpneal/backbone")
  packageStartupMessage("|____/  Beta: type devtools::install_github(\"zpneal/backbone\", ref = \"devel\")")
}
