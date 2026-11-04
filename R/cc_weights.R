#' Clenshaw-Curtis weights
#'
#' @param n Number of points (at least 2).
#'
#' @return Weights for the points returned by \code{cheb_nodes(n)} that
#'   integrate the interpolating polynomial over \code{[-1, 1]}.
#' @keywords internal
#' @noRd
cc_weights <- function(n) {
  k <- 0:(n - 1L)
  moments <- numeric(n)
  even <- k %% 2L == 0L
  moments[even] <- 2 / (1 - k[even]^2)
  drop(crossprod(cheb_coef_matrix(n), moments))
}
