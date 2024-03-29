% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastball.R
\name{fastball}
\alias{fastball}
\title{Randomize a binary matrix using the fastball algorithm}
\usage{
fastball(M, trades = 5 * nrow(M))
}
\arguments{
\item{M}{matrix: a binary matrix (see details)}

\item{trades}{integer: number of trades; the default is 5R trades (approx. mixing time)}
}
\value{
matrix: A random binary matrix with same row sums and column sums as M.
}
\description{
\code{fastball} randomizes a binary matrix, preserving the row and column sums
}
\details{
Given a matrix \code{M}, \code{fastball} randomly samples a new matrix from the space of all matrices with the same row and column sums as \code{M}.
}
\examples{
M <- matrix(rbinom(200,1,0.5),10,20)  #A random 10x20 binary matrix
Mrand <- fastball(M)  #Random matrix with same row and column sums
}
\references{
fastball: {Godard, Karl and Neal, Zachary P. 2022. fastball: A fast algorithm to sample bipartite graphs with fixed degree sequences. \emph{Journal of Complex Networks} \doi{10.1093/comnet/cnac049}}
}
