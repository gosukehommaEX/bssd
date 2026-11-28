#' Chebyshev coefficients of the antiderivative
#'
#' @param A Matrix of Chebyshev coefficients with one column per function
#'   (\code{n} rows, \code{n} at least 2).
#'
#' @return A list with elements \code{B}, the \code{n + 1} by \code{ncol(A)}
#'   matrix of coefficients of the antiderivatives, and \code{base}, their
#'   values at \code{x = -1}.
#' @keywords internal
#' @noRd
cheb_antideriv <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  m <- ncol(A)
  Ap <- rbind(A, matrix(0, 2L, m))
  B <- matrix(0, n + 1L, m)
  B[2L, ] <- Ap[1L, ] - Ap[3L, ] / 2
  # rows 3 to n + 1 at once (each element as in the recurrence
  # B[j + 1, ] = (Ap[j, ] - Ap[j + 2, ]) / (2 j))
  j <- 2L:n
  B[j + 1L, ] <- (Ap[j, , drop = FALSE] - Ap[j + 2L, , drop = FALSE]) /
    (2 * j)
  sgn <- (-1)^(0:n)
  list(B = B, base = colSums(B * sgn))
}
