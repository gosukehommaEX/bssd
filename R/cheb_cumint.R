#' Cumulative integral of Chebyshev interpolants
#'
#' @param anti Output of \code{cheb_antideriv}.
#' @param xq Upper limits on \code{[-1, 1]}, one per column of \code{anti$B};
#'   values outside the interval are moved to the nearest end.
#'
#' @return The integrals from \code{-1} to \code{xq}, one per column.
#' @keywords internal
#' @noRd
cheb_cumint <- function(anti, xq) {
  B <- anti$B
  deg <- 0:(nrow(B) - 1L)
  theta <- acos(pmin(pmax(xq, -1), 1))
  Tq <- cos(outer(deg, theta))
  colSums(B * Tq) - anti$base
}
