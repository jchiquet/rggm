#' Sample of a binary Stochastic Block Model (SBM)
#'
#' This function samples a realization of a stochastic block model from a set a user-defined parameters in a igraĥ format.
#'
#' @param nbNodes number of nodes in the graph
#' @param connectParam inter and intra block connection probabilities
#' @param blockProp proportion of node in each block of the SBM
#'
#' @return an igraph object with a vertex attribute called "memberships" for block belonging
#'
#' @examples
#' ## graph parameters
#' nNodes  <- 90
#' blockProp <- c(.5, .25, .25)   # group proportions
#' nbBlock   <- length(blockProp) # number of blocks
#' connectParam <- diag(.4, nbBlock) + 0.05 # connectivity matrix: affiliation network
#'
#' ## Graph Sampling
#' mySBM <- rSBM(nNodes, connectParam, blockProp)
#'
#' ## Graph plotting
#' plot(mySBM, vertex.color = igraph::V(mySBM)$memberships)
#'
#' @importFrom igraph sample_sbm set_vertex_attr plot.igraph V
#' @importFrom stats rmultinom
#' @export
rSBM <- function(nbNodes, connectParam, blockProp) {
  sizes <- as.numeric(stats::rmultinom(1, nbNodes, blockProp))
  mySBM <- igraph::sample_sbm(nbNodes, connectParam, sizes)
  mySBM <- igraph::set_vertex_attr(mySBM, "memberships", value = rep(1:length(sizes), sizes))
  mySBM
}

#' Sample random weights on edges of a binary Stochastic Block Model
#'
#' This function draws weights according to a user-defined distribution for the edges of an existing binary SBM.
#'
#' @param anSBM SBM sampled from the rSBM function with nbBlock
#' @param family character describing the distribution used for the weigths
#' @param theta list embedding parameters required for the distribution of the weights. See details.
#'
#' Elements in the \code{theta} should be named as follows, depending on the argument \code{family}, with
#' dimension matching the number of blocks in the original SBM:
#' \itemize{
#'   \item{Gaussian}{\code{theta$mu}, a (nbBlock x nbBlock) matrix of means; \code{theta$sigma}, a (nbBlock x nbBlock) matrix of standa deviation}
#'   \item{Poisson}{\code{theta$lambda}, , a (nbBlock x nbBlock) matrix of means.}
#'   \item{Others}{zoon available.}
#' }
#'
#' @return an SBM weight weigthed edges
#'
#' @examples
#' ## graph parameters
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25)   # group proportions
#' nbBlock   <- length(blockProp) # number of blocks
#' connectParam <- diag(.4, nbBlock) + 0.05 # connectivity matrix: affiliation network
#'
#' ## Graph Sampling
#' mySBM <- rSBM(nbNodes, connectParam, blockProp)
#'
#' ## Sampling Gaussian weights
#' mu_within   <- 4 ; sigma_within  <- .5
#' mu_between  <- 2 ; sigma_between <- .5
#' theta <- list()
#' theta$mu    <- matrix(mu_between   , nbBlock, nbBlock); diag(theta$mu)    <- mu_within    # means
#' theta$sigma <- matrix(sigma_between, nbBlock, nbBlock); diag(theta$sigma) <- sigma_within # sd
#'
#' mySBM <- rWeightSBM(mySBM, "gaussian", theta)
#'
#' hist(igraph::E(mySBM)$weight,  breaks = sqrt(igraph::gsize(mySBM)))
#'
#' @importFrom stats rnorm rpois
#' @importFrom igraph set_edge_attr get.data.frame graph_attr V E gsize
#' @export
rWeightSBM <- function(anSBM, family = c("gaussian", "poisson"), theta) {

  ## initializations and routines checks
  stopifnot(igraph::graph_attr(anSBM)$name == "Stochastic block-model")
  stopifnot(!is.null(igraph::V(anSBM)$memberships))
  family <- match.arg(family)
  stopifnot(family %in% c("gaussian", "poisson"))
  nbEdges  <- igraph::gsize(mySBM)
  nbBlocks <- length(unique(igraph::V(mySBM)$memberships))
  stopifnot(all(sapply(theta, dim) == nbBlocks))

  ## Extracting memberships per edges (pair of node)
  memberships <- igraph::V(mySBM)$memberships
  edges       <- igraph::get.data.frame(mySBM, what = "edges")
  memberships_pairs <- cbind(memberships[edges[, 1]], memberships[edges[, 2]])

  ## drawing Gaussian weigths conditional on node's membership
  weights <- switch(family,
                  "gaussian" = rnorm(nbEdges, theta$mu[memberships_pairs] , theta$sigma[memberships_pairs]),
                  "poisson"  = rpois(nbEdges, theta$lambda[memberships_pairs])
  )

  mySBM <- igraph::set_edge_attr(mySBM, "weight", value = weights)
  mySBM
}
