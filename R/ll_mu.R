#' Location parameter of the log-logistic distribution for a given incidence
#'
#' Returns the location parameter \eqn{\mu}{mu} for which the log-logistic
#' cumulative incidence (see \code{\link{ll_cdf}}) at time \code{t} equals
#' \code{ip}:
#' \deqn{\mu = \log t - \sigma \log[ip / (1 - ip)].}{mu = log(t) - sigma * log(ip / (1 - ip)).}
#'
#' @param ip Cumulative incidence at time \code{t}, in \code{(0, 1)}.
#' @param t Time (positive).
#' @param sigma Scale parameter (positive).
#'
#' @return Numeric vector of location parameters. The arguments are recycled
#'   to a common length.
#'
#' @seealso \code{\link{ll_cdf}}
#'
#' @examples
#' # 3.5% incidence at week 104 with sigma = 0.8 gives mu close to 7.3
#' ll_mu(0.035, t = 104, sigma = 0.8)
#'
#' @export
ll_mu <- function(ip, t, sigma) {
  check_num(ip, "ip", len = NULL, lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  check_num(t, "t", len = NULL, lower = 0, lower_open = TRUE)
  check_num(sigma, "sigma", len = NULL, lower = 0, lower_open = TRUE)
  log(t) - sigma * qlogis(ip)
}
