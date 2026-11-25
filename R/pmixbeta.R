#' Distribution function of a mixture of Beta distributions
#'
#' @param q Quantiles.
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#'
#' @return Probabilities.
#' @keywords internal
#' @noRd
pmixbeta <- function(q, mix) {
  vapply(q, function(x) sum(mix$w * pbeta(x, mix$shape1, mix$shape2)),
         numeric(1))
}
