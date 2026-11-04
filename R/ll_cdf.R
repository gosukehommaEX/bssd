#' Cumulative incidence of the log-logistic distribution
#'
#' Computes the cumulative incidence (distribution function) of the
#' log-logistic time-to-event distribution in the parameterization used by
#' Duncan (2026), which is the parameterization of the log-logistic model of
#' the SAS LIFEREG procedure:
#' \deqn{F(t) = 1 - \frac{1}{1 + \exp[(\log t - \mu) / \sigma]}.}{F(t) = 1 - 1 / (1 + exp((log(t) - mu) / sigma)).}
#'
#' @param t Numeric vector of times (positive).
#' @param mu Location (intercept) parameter on the log-time scale.
#' @param sigma Scale parameter (positive). The hazard is decreasing when
#'   \code{sigma >= 1} and first increases and then decreases when
#'   \code{sigma < 1}.
#'
#' @return Numeric vector of cumulative incidences. The arguments are
#'   recycled to a common length.
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{ll_mu}}
#'
#' @examples
#' # incidence at week 104 under the placebo assumption of scenario 6
#' ll_cdf(104, mu = 7.3, sigma = 0.8)
#'
#' @export
ll_cdf <- function(t, mu, sigma) {
  check_num(t, "t", len = NULL, lower = 0, lower_open = TRUE)
  check_num(mu, "mu", len = NULL)
  check_num(sigma, "sigma", len = NULL, lower = 0, lower_open = TRUE)
  plogis((log(t) - mu) / sigma)
}
