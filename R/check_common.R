#' Check the arguments shared by the posterior functions
#'
#' @param rr0 Critical value of the ratio.
#' @param cutoff Posterior probability threshold.
#' @param level Probability of the credible intervals.
#' @param details Logical flag.
#'
#' @return \code{NULL}, invisibly; an error is raised for invalid values.
#' @keywords internal
#' @noRd
check_common <- function(rr0, cutoff, level, details) {
  check_num(rr0, "rr0", lower = 0, lower_open = TRUE)
  check_num(cutoff, "cutoff", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  check_num(level, "level", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("'details' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}
