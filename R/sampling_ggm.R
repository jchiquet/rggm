#' Graph to Precision matrix
#'
#' From a (possibility weighted) igraph object, build a precision matrix from the Laplacian matrix of the graph,
#' which is set strictly positive-definite by iteratively adding small constant to its diagonal,
#' proportionally to the degree of each node.
#'
#' @param graph an igraph object
#' @param neg_prop double, the proportion of negative signs in the target precision matrix. Default is 0.5
#' @param cond_var a target vector of conditional variances (which equal the inverse of the diaognal in a precision matrix).
#' When NULL (the default), approximately equal to 1/degrees(graph).
#' @param epsilon double, the minimal eigen values to reach in the precision matrix. Default to 1e-2.
#' @param delta double, the quantity by which the diagonal is increased at each iteration. Default to 1e-1
#' @param maxIter integer for the maximal number of iteration to reach the target minimal eigen value. Default to 1e4
#'
#' @return a precision matrix with Matrix format
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
#' ## Precision matrix
#' Omega <- graph2prec(mySBM)
#'
#' @importFrom Matrix diag Diagonal
#' @importFrom igraph laplacian_matrix degree is_igraph
#' @export
graph2prec <- function(graph, neg_prop = .5, cond_var = NULL, epsilon = 1e-2, delta = 1e-1, maxIter = 1e4) {

  ## routine checks
  stopifnot(is_igraph(graph))

  ## initialization
  Omega   <- igraph::laplacian_matrix(graph)
  degrees <- igraph::degree(graph)
  nbNodes <- length(degrees)
  iter <- 0

  ## iterative to reach diagonal dominance
  while (min(eigen(Omega, only.values = TRUE)$values) < epsilon & iter < maxIter) {
    iter <- iter + 1
    diag(Omega) <- (1 + delta)^iter * degrees
  }

  if (!is.null(cond_var)) {
    stopifnot(length(cond_var) == nbNodes)
    D <- Diagonal(x = 1/(sqrt(diag(Omega) * cond_var)))
    Omega <- D %*% Omega %*% D
  }

  signs <- matrix(sample(c(1, -1), nbNodes*2, replace = TRUE, prob = c(1 - neg_prop, neg_prop)), nbNodes, nbNodes)
  signs <- signs * t(signs)
  diag(signs) <- 1
  Omega <- Omega * signs

  Omega
}

