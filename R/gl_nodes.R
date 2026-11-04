#' Gauss-Legendre nodes and weights
#'
#' Computes the nodes and weights by the Golub-Welsch algorithm.
#'
#' @param n Number of nodes.
#'
#' @return A list with elements \code{x} (nodes on \code{[-1, 1]} in
#'   increasing order) and \code{w} (weights).
#' @keywords internal
#' @noRd
gl_nodes <- function(n) {
  if (n == 1L) {
    return(list(x = 0, w = 2))
  }
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i * i - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = 2 * e$vectors[1L, ord]^2)
}
