% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparsify.R
\name{sparsify.with.skeleton}
\alias{sparsify.with.skeleton}
\title{Extract Karger's (1999) skeleton backbone}
\usage{
sparsify.with.skeleton(U, s, class = "original", narrative = FALSE)
}
\arguments{
\item{U}{An unweighted unipartite graph, as: (1) an adjacency matrix in the form of a matrix or sparse \code{\link{Matrix}}; (2) an edgelist in the form of a two-column dataframe; (3) an \code{\link{igraph}} object.}

\item{s}{numeric: Proportion of edges to retain, 0 < s < 1; smaller values yield sparser graphs}

\item{class}{string: the class of the returned backbone graph, one of c("original", "matrix", "Matrix", "igraph", "edgelist").
If "original", the backbone graph returned is of the same class as \code{U}.}

\item{narrative}{boolean: TRUE if suggested text & citations should be displayed.}
}
\value{
An unweighted, undirected, unipartite graph of class \code{class}.
}
\description{
\code{sparsify.with.skeleton} is a wrapper for \code{\link[=sparsify]{sparsify()}} that extracts the skeleton backbone described by Karger (1999),
which preserves a specified proportion of random edges. It is equivalent to \code{sparsify(escore = "random", normalize = "none", filter = "proportion", umst = FALSE)}.
}
\examples{
U <- igraph::erdos.renyi.game(60, .5)
plot(U) #A dense graph
sparse <- sparsify.with.skeleton(U, s = 0.25, narrative = TRUE)
plot(sparse) #A sparser graph
}
\references{
package: {Neal, Z. P. (2022). backbone: An R Package to Extract Network Backbones. \emph{PLOS ONE, 17}, e0269137. \doi{10.1371/journal.pone.0269137}}

model: {Karger, D. R. (1999). Random sampling in cut, flow, and network design problems. \emph{Mathematics of Operations Research, 24}, 383-413. \doi{10.1287/moor.24.2.383}}
}
