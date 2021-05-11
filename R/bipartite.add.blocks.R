#' Adds a block structure to a bipartite network
#'
#' `bipartite.add.blocks` rewires a bipartite graph B to have a block structure such that edges are
#' located within-block with `density` probability, while preserving the degree distributions.
#'
#' @param B A bipartite network object of class "matrix", "sparseMatrix", \link{igraph}, matrix or dataframe edgelist, or \link[network]{network}
#' @param blocks numeric: number of blocks to add (between 2 and 26)
#' @param density numeric: desired within-block density
#'
#' @return
#' @export
#'
#' @examples
#' B <- bipartite.from.probability(R = 10, C = 10, P = .5)
#' B <- bipartite.add.blocks(B, blocks = 2, density = .7)
bipartite.add.blocks <- function(B,blocks=2,density=.5) {

  #Convert supplied object to matrix
  B <- suppressMessages(tomatrix(B))
  class <- B$summary$class
  B <- B$G

  # Parameter checks
  if (!is.numeric(blocks) | !is.numeric(density)) {stop("blocks and density must be numeric")}
  if (blocks<2 | blocks>26 | blocks%%1!=0) {stop("blocks must be a positive integer between 2 and 26")}
  if (density<0.5 | density>1) {stop("density must be between 0 and 1")}

  # Begin progress bar, assign nodes to blocks
  block.names <- LETTERS[seq(from = 1, to = blocks)]  #Generate list of group names
  rownames(B) <- paste0(sample(block.names,nrow(B),replace=TRUE),c(1:nrow(B)))  #Assign each agent to a group
  colnames(B) <- paste0(sample(block.names,ncol(B),replace=TRUE),c(1:ncol(B)))  #Assign each artifact to a group
  within.block <- sum((outer(substr(rownames(B),1,1), substr(colnames(B),1,1), `==`)*1)*B) / sum(B)  #Compute starting block density
  pb <- txtProgressBar(min = .49, max = density, style = 3)  #Initiate progress bar

  while (within.block < density) {
    # Pick agents
    agent1 <- sample(rownames(B),1)  #Pick a random agent
    agent2 <- sample(rownames(B)[which(substr(rownames(B),1,1)!=substr(agent1,1,1))],1)  #Pick a random agent from another group

    # Pick artifacts
    artifact1 <- sample(colnames(B)[which(substr(colnames(B),1,1)==substr(agent1,1,1))],1)  #Pick a random artifact from agent 1's group
    artifact2 <- sample(colnames(B)[which(substr(colnames(B),1,1)==substr(agent2,1,1))],1)  #Pick a random artifact from agent 2's group

    # If a checkboard swap would increase within-block density, make the swap and recompute
    if (all(matrix(c(0,1,1,0),nrow=2,ncol=2) == B[c(agent1,agent2),c(artifact1,artifact2)])) {
      B[c(agent1,agent2),c(artifact1,artifact2)] <- abs(B[c(agent1,agent2),c(artifact1,artifact2)] - 1)
      within.block <- sum((outer(substr(rownames(B),1,1), substr(colnames(B),1,1), `==`)*1)*B) / sum(B)
      setTxtProgressBar(pb, within.block)
    }
  }

  # Arrange by groups and close progress bar
  close(pb)
  B <- B[order(rownames(B)), order(colnames(B))]
  B <- frommatrix(B, class)
  return(B)
}
