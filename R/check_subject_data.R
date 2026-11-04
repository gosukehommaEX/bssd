#' Check subject-level event indicators and follow-up times
#'
#' @param y Event indicators.
#' @param t Follow-up times.
#'
#' @return A list with elements \code{y} (integer vector) and \code{t}
#'   (numeric vector).
#' @keywords internal
#' @noRd
check_subject_data <- function(y, t) {
  if (is.logical(y)) {
    y <- as.integer(y)
  }
  check_num(y, "y", len = NULL, lower = 0, upper = 1, integer = TRUE)
  check_num(t, "t", len = length(y), lower = 0, lower_open = TRUE)
  list(y = as.integer(round(y)), t = as.numeric(t))
}
