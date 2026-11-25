#' Convert a prior specification to a mixture of Beta distributions
#'
#' @param x Numeric vector \code{c(a, b)} of the shape parameters of one Beta
#'   distribution, or a data frame or matrix with columns \code{w},
#'   \code{shape1} and \code{shape2} (one row per component), such as the
#'   result of \code{\link{map_beta}}.
#' @param name Name of the argument for error messages.
#'
#' @return A data frame with columns \code{w} (weights summing to 1),
#'   \code{shape1} and \code{shape2}.
#' @keywords internal
#' @noRd
as_mixbeta <- function(x, name) {
  if (is.numeric(x) && is.null(dim(x))) {
    check_num(x, name, len = 2L, lower = 0, lower_open = TRUE)
    return(data.frame(w = 1, shape1 = x[1L], shape2 = x[2L]))
  }
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x) || !all(c("w", "shape1", "shape2") %in% names(x))) {
    stop("'", name, "' must be c(shape1, shape2) or a data frame with ",
         "columns w, shape1 and shape2.", call. = FALSE)
  }
  check_num(x$w, paste0(name, "$w"), len = NULL, lower = 0)
  check_num(x$shape1, paste0(name, "$shape1"), len = nrow(x), lower = 0,
            lower_open = TRUE)
  check_num(x$shape2, paste0(name, "$shape2"), len = nrow(x), lower = 0,
            lower_open = TRUE)
  if (sum(x$w) <= 0) {
    stop("The weights in '", name, "' must not all be zero.", call. = FALSE)
  }
  data.frame(w = x$w / sum(x$w), shape1 = x$shape1, shape2 = x$shape2)
}
