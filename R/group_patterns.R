#' Group subjects with the same data
#'
#' Subjects with the same event value and the same time contribute the same
#' factor to the likelihood, so they are counted once with a weight.
#'
#' @param y Event values (counts or indicators).
#' @param t Times.
#'
#' @return A list with the distinct pairs \code{y}, \code{t} and their
#'   numbers \code{w}.
#' @keywords internal
#' @noRd
group_patterns <- function(y, t) {
  o <- order(y, t)
  ys <- y[o]
  ts <- t[o]
  n <- length(ys)
  new <- c(TRUE, ys[-1L] != ys[-n] | ts[-1L] != ts[-n])
  list(y = as.numeric(ys[new]), t = as.numeric(ts[new]),
       w = as.numeric(tabulate(cumsum(new))))
}
