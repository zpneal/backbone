#' Bill Sponsorship in the 108th U.S. Senate
#'
#' A bipartite network describing bill sponsorship in the 108th U.S. Senate
#'
#' @format ## `senate108`
#' An \link[igraph]{igraph} object containing a bipartite network with 100 senators and 3035 bills,
#'   where a senator is connected to a bill they sponsored or co-sponsored. Senator nodes include
#'   attributes for name, party affiliation, state, Govtrack ID, color (red for Republican, blue
#'   for Democrat, green for independent). Bill nodes include attributes for date introduced,
#'   title, legislative policy area, party affiliation of the sponsor, partisanship, and most
#'   recent status.

#' @source These data were generated using `incidentally::incidence.from.congress(session = 108, types = "s", format = "igraph")`
"senate108"
