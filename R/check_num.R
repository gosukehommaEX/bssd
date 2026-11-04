#' Check a numeric argument
#'
#' Stops with an informative message when \code{x} is not a finite numeric
#' vector of the required length and range.
#'
#' @param x Value to check.
#' @param name Name of the argument used in the error message.
#' @param len Required length, or \code{NULL} for any positive length.
#' @param lower Lower limit.
#' @param upper Upper limit.
#' @param lower_open Logical; if \code{TRUE}, \code{x} must be strictly
#'   greater than \code{lower}.
#' @param upper_open Logical; if \code{TRUE}, \code{x} must be strictly
#'   less than \code{upper}.
#' @param integer Logical; if \code{TRUE}, \code{x} must be integer valued.
#'
#' @return \code{x}, invisibly.
#' @keywords internal
#' @noRd
check_num <- function(x, name, len = 1L, lower = -Inf, upper = Inf,
                      lower_open = FALSE, upper_open = FALSE,
                      integer = FALSE) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop("'", name, "' must be a finite numeric value.", call. = FALSE)
  }
  if (!is.null(len) && length(x) != len) {
    stop("'", name, "' must have length ", len, ".", call. = FALSE)
  }
  low_bad <- if (lower_open) any(x <= lower) else any(x < lower)
  up_bad <- if (upper_open) any(x >= upper) else any(x > upper)
  if (low_bad || up_bad) {
    stop("'", name, "' is out of range.", call. = FALSE)
  }
  if (integer && any(abs(x - round(x)) > 1e-8)) {
    stop("'", name, "' must be integer valued.", call. = FALSE)
  }
  invisible(x)
}
