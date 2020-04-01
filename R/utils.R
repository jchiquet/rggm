
dlaplace <- function(x, mu = 0, b = 1, log = FALSE) {
  res <- -abs(x-mu)/b - log(2*b)
  if (!log) res <- exp(res)
  res
}

#' @importFrom stats runif
rlaplace <- function(n, mu = 0, b = 1) {
  u <- runif(n, -0.5, 0.5)
  mu - b*sign(u)*log(1-2*abs(u))
}

