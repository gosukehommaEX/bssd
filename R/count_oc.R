#' Probability that the blinded event count reaches a critical value
#'
#' Computes the exact probability that the blinded total event count is at
#' least \code{y_crit} under true event probabilities or rates of study drug
#' and control: the operating characteristics of a rule based on
#' \code{\link{signal_count}} (Gould and Wang, 2017, Figure 3 and Table V).
#'
#' @details
#' With \code{likelihood = "binomial"} and \code{allocation = "fixed"}, the
#' count is the sum of \eqn{\mathrm{Binomial}(N_d, p_d)}{Binomial(N_d, p_d)}
#' and \eqn{\mathrm{Binomial}(N_c, p_c)}{Binomial(N_c, p_c)} counts with
#' \eqn{N_d = \mathrm{round}(k N / (k + 1))}{N_d = round(k * N / (k + 1))}
#' (the model of \code{\link{canary_post}}); with \code{allocation =
#' "random"}, it is \eqn{\mathrm{Binomial}(N, (p_c + k p_d) / (k + 1))}{Binomial(N, (p_c + k * p_d) / (k + 1))}
#' (the model of \code{\link{binbeta_post}}). With \code{likelihood =
#' "poisson"}, \code{n} is the total exposure and the count is
#' \eqn{\mathrm{Poisson}(n (\lambda_c + k \lambda_d) / (k + 1))}{Poisson(n * (lambda_c + k * lambda_d) / (k + 1))}.
#'
#' @param y_crit Critical counts (\code{NA} gives probability 0).
#' @param n Numbers of subjects, or total exposures for the Poisson
#'   likelihood.
#' @param rate_c,rate_d True event probabilities (binomial) or rates
#'   (Poisson) of control and study drug.
#' @param k Randomization ratio of study drug to control (\code{k:1}).
#' @param likelihood \code{"binomial"} or \code{"poisson"}.
#' @param allocation \code{"fixed"} or \code{"random"} arm sizes (binomial
#'   likelihood only).
#'
#' @return Numeric vector of probabilities; the arguments \code{y_crit},
#'   \code{n}, \code{rate_c} and \code{rate_d} are recycled to a common
#'   length.
#'
#' @references
#' Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
#' rate differences using data from blinded trials: the canary in the coal
#' mine. \emph{Statistics in Medicine}, 36(1), 92--104.
#' \doi{10.1002/sim.7129}
#'
#' @seealso \code{\link{signal_count}}
#'
#' @examples
#' # probability of 6 or more events among 100 subjects when the true
#' # probabilities are 0.02 (control) and 0.02, 0.04 or 0.08 (study drug)
#' count_oc(6, 100, rate_c = 0.02, rate_d = c(0.02, 0.04, 0.08))
#'
#' @export
count_oc <- function(y_crit, n, rate_c, rate_d, k = 1,
                     likelihood = c("binomial", "poisson"),
                     allocation = c("fixed", "random")) {
  likelihood <- match.arg(likelihood)
  allocation <- match.arg(allocation)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  len <- max(length(y_crit), length(n), length(rate_c), length(rate_d))
  y_crit <- rep_len(y_crit, len)
  n <- rep_len(n, len)
  rate_c <- rep_len(rate_c, len)
  rate_d <- rep_len(rate_d, len)
  check_num(n, "n", len = NULL, lower = 0, lower_open = TRUE)
  ok <- !is.na(y_crit)
  if (any(ok)) {
    check_num(y_crit[ok], "y_crit", len = NULL, lower = 0, integer = TRUE)
  }
  upper <- if (likelihood == "binomial") 1 else Inf
  check_num(rate_c, "rate_c", len = NULL, lower = 0, upper = upper)
  check_num(rate_d, "rate_d", len = NULL, lower = 0, upper = upper)
  wc <- 1 / (k + 1)
  vapply(seq_len(len), function(i) {
    if (is.na(y_crit[i])) {
      return(0)
    }
    yc <- round(y_crit[i])
    if (likelihood == "poisson") {
      return(ppois(yc - 1, n[i] * (wc * rate_c[i] + (1 - wc) * rate_d[i]),
                   lower.tail = FALSE))
    }
    nn <- round(n[i])
    if (allocation == "random") {
      return(pbinom(yc - 1, nn, wc * rate_c[i] + (1 - wc) * rate_d[i],
                    lower.tail = FALSE))
    }
    n1 <- round(k / (k + 1) * nn)
    x <- 0:n1
    sum(dbinom(x, n1, rate_d[i]) *
          pbinom(yc - x - 1, nn - n1, rate_c[i], lower.tail = FALSE))
  }, numeric(1))
}
