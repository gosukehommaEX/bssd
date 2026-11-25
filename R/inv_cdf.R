#' Quantiles from a distribution function by root finding
#'
#' @param cdf Nondecreasing distribution function.
#' @param p Probabilities.
#' @param lo,hi Interval for the root search.
#' @param extend Logical; if \code{TRUE}, the interval is extended until the
#'   root is bracketed (the function is increasing).
#'
#' @return Numeric vector of quantiles (\code{NA} where the search fails).
#' @keywords internal
#' @noRd
inv_cdf <- function(cdf, p, lo, hi, extend = FALSE) {
  vapply(p, function(q) {
    f <- function(x) cdf(x) - q
    tryCatch(
      if (extend) {
        uniroot(f, c(lo, hi), tol = 1e-12, extendInt = "upX")$root
      } else {
        uniroot(f, c(lo, hi), tol = 1e-12)$root
      },
      error = function(e) NA_real_)
  }, numeric(1))
}
