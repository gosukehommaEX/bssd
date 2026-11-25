#' Posterior probabilities of the BDRIBS model on a grid of data
#'
#' Computes the posterior probability that the rate ratio exceeds
#' \code{rr0} under the BDRIBS model for every combination of numbers of
#' events, exposures and background rates. The results give the contour
#' plots and the sensitivity plots used by Waterhouse et al. (2022) to plan
#' and interpret blinded safety reviews.
#'
#' @param y Numbers of events (non-negative integers).
#' @param exposure Exposures (positive).
#' @param lambda_hp Background event rates of the control group.
#' @param k Randomization ratio of study drug to control (\code{k:1}).
#' @param pt_hp Total patient time of the historical control data, or
#'   \code{Inf} for known background rates (the default, as in the examples
#'   of Waterhouse et al., 2022).
#' @param ... Further arguments of \code{\link{bdribs_post}} (\code{df},
#'   \code{p_prior}, \code{rr0}, \code{cutoff}).
#'
#' @return A data frame with columns \code{y}, \code{exposure},
#'   \code{lambda_hp} and \code{prob}, one row per combination.
#'
#' @references
#' Waterhouse, B., Hartford, A., Mukhopadhyay, S., Ferguson, R. and
#' Hendrickson, B. A. (2022). Using the Bayesian detection of potential risk
#' using inference on blinded safety data (BDRIBS) method to support the
#' decision to refer an event for unblinded evaluation. \emph{Pharmaceutical
#' Statistics}, 21(2), 372--385. \doi{10.1002/pst.2175}
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{pooled_ratio}}
#'
#' @examples
#' # Example 1 of Waterhouse et al. (2022): 6, 7 or 8 events, background rate
#' # 3.1 per 100 patient-years, randomization 1:1
#' g <- bdribs_grid(y = 6:8, exposure = seq(60, 200, by = 20),
#'                  lambda_hp = 0.031, p_prior = c(0.5, 0.5))
#' head(g)
#'
#' @export
bdribs_grid <- function(y, exposure, lambda_hp, k = 1, pt_hp = Inf, ...) {
  check_num(y, "y", len = NULL, lower = 0, integer = TRUE)
  check_num(exposure, "exposure", len = NULL, lower = 0, lower_open = TRUE)
  check_num(lambda_hp, "lambda_hp", len = NULL, lower = 0, lower_open = TRUE)
  args <- list(...)
  if ("details" %in% names(args)) {
    stop("'details' cannot be given.", call. = FALSE)
  }
  grid <- expand.grid(y = y, exposure = exposure, lambda_hp = lambda_hp,
                      KEEP.OUT.ATTRS = FALSE)
  grid$prob <- vapply(seq_len(nrow(grid)), function(i) {
    do.call(bdribs_post,
            c(list(y = grid$y[i], exposure = grid$exposure[i], k = k,
                   lambda_hp = grid$lambda_hp[i], pt_hp = pt_hp,
                   details = FALSE), args))$prob
  }, numeric(1))
  grid
}
