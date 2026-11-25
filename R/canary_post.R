#' Posterior of the risk elevation from the total event count
#'
#' Computes the posterior distribution of a metric of the difference in
#' event risk between study drug and control (rate ratio, risk difference or
#' odds ratio) from the blinded total event count, with the method of Gould
#' and Wang (2017), which uses external information about the control event
#' rate.
#'
#' @details
#' \strong{Binomial likelihood.} Of \eqn{N}{N} subjects,
#' \eqn{N_d = \mathrm{round}(k N / (k + 1))}{N_d = round(k * N / (k + 1))} are
#' on study drug and \eqn{N_c = N - N_d}{N_c = N - N_d} on control (the
#' expected arm sizes, as in Gould and Wang, 2017). The total count
#' \eqn{X}{X} is the sum of independent
#' \eqn{\mathrm{Binomial}(N_d, p_d)}{Binomial(N_d, p_d)} and
#' \eqn{\mathrm{Binomial}(N_c, p_c)}{Binomial(N_c, p_c)} counts, and its
#' probability is the convolution of the two. The control probability has the
#' prior \eqn{p_c \sim \mathrm{Beta}(a, b)}{p_c ~ Beta(a, b)}. The metric is
#' \itemize{
#'   \item \code{"ratio"}: \eqn{R = p_d / p_c}{R = p_d / p_c} with a
#'     lognormal prior truncated to \eqn{R \le 1 / p_c}{R <= 1 / p_c}, as in
#'     the code of Gould and Wang (2017);
#'   \item \code{"odds"}: the odds ratio \eqn{Q}{Q} with a lognormal prior;
#'   \item \code{"difference"}: \eqn{D = p_d - p_c}{D = p_d - p_c} with a
#'     uniform prior on \eqn{(-p_c, 1 - p_c)}{(-p_c, 1 - p_c)}, which is a
#'     uniform prior on \eqn{p_d}{p_d} independent of \eqn{p_c}{p_c}.
#' }
#' The two-dimensional posterior is integrated on
#' \eqn{\mathrm{logit}(p_c)}{logit(p_c)} (Gauss-Legendre quadrature) and the
#' logarithm of the metric, or \eqn{p_d}{p_d} for the difference (Chebyshev
#' interpolation). For each node of \eqn{\mathrm{logit}(p_c)}{logit(p_c)},
#' the interval of the second coordinate that holds the conditional
#' posterior mass is located separately; for the ratio it ends at the
#' truncation point \eqn{R = 1 / p_c}{R = 1 / p_c}. The Gauss-Legendre rule
#' is split where the posterior probability of the metric has a kink
#' (\eqn{p_c = 1 / crit}{p_c = 1 / crit} for the ratio,
#' \eqn{p_c = 1 - crit}{p_c = 1 - crit} or \eqn{p_c = -crit}{p_c = -crit}
#' for the difference), so that the results stay accurate when the event
#' probabilities are large.
#'
#' \strong{Poisson likelihood.} With the total exposure \eqn{E}{E} (\code{n})
#' and the control event rate \eqn{\lambda_c}{lambda_c} per unit of exposure,
#' \eqn{X \sim \mathrm{Poisson}(E \lambda_c (1 + k R) / (k + 1))}{X ~ Poisson(E * lambda_c * (1 + k * R) / (k + 1))}
#' with \eqn{\lambda_c \sim \mathrm{Gamma}(a, b)}{lambda_c ~ Gamma(a, b)}
#' (shape and rate) and a lognormal prior on \eqn{R}{R}. The control rate is
#' integrated out analytically and the posterior of \eqn{\log R}{log(R)} is
#' integrated numerically. Gould and Wang (2017) write the model with the
#' expected number of control events; their Gamma prior with rate
#' \eqn{b'}{b'} on that number corresponds to \eqn{b = b' E / (k + 1)}{b = b' * E / (k + 1)}.
#' The likelihood is that of the BDRIBS model (\code{\link{bdribs_post}});
#' only the prior of the rate ratio differs.
#'
#' The default \code{metric_prior} is the lognormal distribution with mean 1
#' and probability 0.95 between 0.1 and 10 (Gould and Wang, 2017, Section
#' 3.1); \code{elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)} gives its
#' parameters. A signal is declared when
#' \eqn{P(M > crit \mid X) > }{P(M > crit | X) > } \code{cutoff}.
#'
#' @param y Total number of subjects with the event (binomial) or of events
#'   (Poisson) in the blinded data.
#' @param n Number of subjects (binomial) or total exposure (Poisson).
#' @param k Randomization ratio of study drug to control (\code{k:1}).
#' @param prior_c Parameters of the prior of the control event probability
#'   (\code{c(a, b)} of a Beta distribution, binomial) or rate (shape and rate
#'   of a Gamma distribution, Poisson).
#' @param metric One of \code{"ratio"}, \code{"difference"} and
#'   \code{"odds"}; only \code{"ratio"} is available for the Poisson
#'   likelihood.
#' @param metric_prior \code{c(meanlog, sdlog)} of the lognormal prior of the
#'   ratio or odds ratio, or \code{NULL} for the default (see Details).
#'   Ignored for the difference.
#' @param likelihood \code{"binomial"} or \code{"poisson"}.
#' @param crit Critical value of the metric; by default 1 for the ratio and
#'   the odds ratio and 0 for the difference.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the equal-tailed credible interval of the
#'   metric.
#' @param details Logical; if \code{FALSE}, the quantiles are not computed
#'   (they are \code{NA}).
#' @param control List of settings of the numerical integration, as in
#'   \code{\link{exptime_post}} (\code{n_u} is not used for the Poisson
#'   likelihood).
#'
#' @return An object of class \code{bssd_post} with model \code{"CANARY"},
#'   \code{prob} (posterior probability that the metric exceeds \code{crit}),
#'   \code{signal}, \code{rr} (posterior mean, median and credible interval of
#'   the metric), \code{estimates} (posterior means of the control and study
#'   drug probabilities or rates and of the metric), labels and settings,
#'   \code{data}, \code{prior} and \code{integration}.
#'
#' @references
#' Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
#' rate differences using data from blinded trials: the canary in the coal
#' mine. \emph{Statistics in Medicine}, 36(1), 92--104.
#' \doi{10.1002/sim.7129}
#'
#' @seealso \code{\link{signal_count}} and \code{\link{count_oc}} for the
#'   critical event counts and their operating characteristics.
#'
#' @examples
#' # 12 subjects with the event among 100, control prior Beta(1, 49)
#' canary_post(y = 12, n = 100, prior_c = c(1, 49))
#' # the same data with the risk difference
#' canary_post(y = 12, n = 100, prior_c = c(1, 49), metric = "difference")
#'
#' @export
canary_post <- function(y, n, k = 1, prior_c,
                        metric = c("ratio", "difference", "odds"),
                        metric_prior = NULL,
                        likelihood = c("binomial", "poisson"), crit = NULL,
                        cutoff = 0.95, level = 0.95, details = TRUE,
                        control = list()) {
  metric <- match.arg(metric)
  likelihood <- match.arg(likelihood)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(prior_c, "prior_c", len = 2L, lower = 0, lower_open = TRUE)
  if (likelihood == "binomial") {
    check_num(n, "n", lower = 1, integer = TRUE)
    check_num(y, "y", lower = 0, upper = n, integer = TRUE)
  } else {
    if (metric != "ratio") {
      stop("Only metric = \"ratio\" is available for the Poisson likelihood.",
           call. = FALSE)
    }
    check_num(n, "n", lower = 0, lower_open = TRUE)
    check_num(y, "y", lower = 0, integer = TRUE)
  }
  if (is.null(crit)) {
    crit <- if (metric == "difference") 0 else 1
  }
  if (metric == "difference") {
    check_num(crit, "crit", lower = -1, upper = 1, lower_open = TRUE,
              upper_open = TRUE)
  } else {
    check_num(crit, "crit", lower = 0, lower_open = TRUE)
    if (is.null(metric_prior)) {
      metric_prior <- elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)
    }
    check_num(metric_prior, "metric_prior", len = 2L)
    if (metric_prior[2L] <= 0) {
      stop("The second element of 'metric_prior' (sdlog) must be positive.",
           call. = FALSE)
    }
  }
  check_common(1, cutoff, level, details)
  ctrl <- merge_control(control, list(n_coarse = 41L, n_u = 48L, n_v = 49L,
                                      check = TRUE, eps = 1e-8,
                                      max_iter = 12L))
  y <- round(y)
  if (likelihood == "binomial") {
    n <- round(n)
    res <- canary_binom(y, n, k, prior_c, metric, metric_prior, crit, level,
                        details, ctrl)
  } else {
    res <- canary_pois(y, n, k, prior_c, metric_prior, crit, level, details,
                       ctrl)
  }
  lab <- c(ratio = "R", difference = "D", odds = "Q")[[metric]]
  text <- c(ratio = "ratio R = p_d / p_c of the event probabilities",
            difference = "difference D = p_d - p_c of the event probabilities",
            odds = "odds ratio Q of the event probabilities")[[metric]]
  if (likelihood == "poisson") {
    text <- "rate ratio R = lambda_d / lambda_c"
  }
  out <- list(
    model = "CANARY",
    prob = res$prob,
    signal = res$prob > cutoff,
    rr = c(mean = res$mean_m, median = res$q[2L], lower = res$q[1L],
           upper = res$q[3L]),
    estimates = res$means,
    incidence = NULL,
    rr_label = lab,
    rr0 = crit,
    description = paste("Tested quantity:", text),
    cutoff = cutoff,
    level = level,
    t_star = NULL,
    data = list(n = n, n_event = y, k = k, likelihood = likelihood),
    prior = list(prior_c = prior_c, metric = metric,
                 metric_prior = metric_prior),
    integration = list(method = res$method, nodes = res$nodes,
                       box = res$box, error = res$error)
  )
  class(out) <- "bssd_post"
  out
}
