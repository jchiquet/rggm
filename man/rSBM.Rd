% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_sbm.R
\name{rSBM}
\alias{rSBM}
\title{Sample of a binary Stochastic Block Model (SBM)}
\usage{
rSBM(nbNodes, connectParam, blockProp)
}
\arguments{
\item{nbNodes}{number of nodes in the graph}

\item{connectParam}{inter and intra block connection probabilities}

\item{blockProp}{proportion of node in each block of the SBM}
}
\value{
an igraph object with a vertex attribute called "memberships" for block belonging
}
\description{
This function samples a realization of a stochastic block model from a set a user-defined parameters in a igraĥ format.
}
\examples{
## graph parameters
nNodes  <- 90
blockProp <- c(.5, .25, .25)   # group proportions
nbBlock   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlock) + 0.05 # connectivity matrix: affiliation network

## Graph Sampling
mySBM <- rSBM(nNodes, connectParam, blockProp)

## Graph plotting
plot(mySBM, vertex.color = igraph::V(mySBM)$memberships)

}
