#' Posterior of the rate ratio from blinded event counts of several studies
#'
#' Computes the posterior distribution of a rate ratio common to several
#' blinded studies from the total numbers of subjects with the event in each
#' study, with the Bayesian model of Goren et al. (2020).
#'
#' @details
#' Study \eqn{s}{s} has \eqn{N_s}{N_s} enrolled subjects, \eqn{Y_s}{Y_s}
#' subjects with the event, and the probability
#' \eqn{r_s = k_s / (k_s + 1)}{r_s = k_s / (k_s + 1)} of study drug. The
#' model is
#' \deqn{Y_s \sim \mathrm{Binomial}(N_s, \theta_s (1 - r_s + r_s \delta)), \quad \theta_s \sim \mathrm{Beta}(a_s, b_s), \quad \log \delta \sim N(m, v^2),}{Y_s ~ Binomial(N_s, theta_s * (1 - r_s + r_s * delta)), theta_s ~ Beta(a_s, b_s), log(delta) ~ N(m, v^2),}
#' restricted to \eqn{\theta_s \le 1 / \delta}{theta_s <= 1 / delta}, so
#' that the study drug probability \eqn{\delta \theta_s}{delta * theta_s} is
#' at most 1. The control probabilities \eqn{\theta_s}{theta_s} may differ
#' between studies, and \eqn{\delta}{delta} is common. Given
#' \eqn{\delta}{delta}, each \eqn{\theta_s}{theta_s} is integrated out by
#' composite Gauss-Legendre quadrature on
#' \eqn{\mathrm{logit}(\theta_s)}{logit(theta_s)} over the region where the
#' integrand is not negligible, so the posterior of \eqn{\delta}{delta} is
#' one-dimensional; it is integrated at Chebyshev points of the prior
#' distribution function of \eqn{\log \delta}{log(delta)}. Goren et al.
#' (2020) sample the same posterior by Hamiltonian Monte Carlo.
#'
#' A signal is declared when
#' \eqn{P(\delta > rr_0 \mid \mathrm{data}) > }{P(delta > rr0 | data) > }
#' \code{cutoff}. With a single study, the likelihood is that of
#' \code{\link{binbeta_post}} with \eqn{\theta_d = \delta \theta_c}{theta_d = delta * theta_c}.
#'
#' The posterior mean of \eqn{\delta}{delta} is not reported, because with
#' the vague prior recommended by Goren et al. (2020) (\eqn{v = 5}{v = 5}) it
#' depends on the far right tail of the prior; the posterior median and the
#' credible interval are reported instead, and \code{estimates} holds the
#' posterior mean of \eqn{\log \delta}{log(delta)}.
#'
#' @param y Numbers of subjects with the event, one per study.
#' @param n Numbers of enrolled subjects, one per study.
#' @param k Randomization ratios of study drug to control, one per study or
#'   one for all.
#' @param prior_c Shape parameters of the Beta prior of the control
#'   probabilities: \code{c(a, b)} for all studies, or a matrix with one row
#'   \code{c(a_s, b_s)} per study.
#' @param lrr_prior Mean and standard deviation \code{c(m, v)} of the normal
#'   prior of \eqn{\log \delta}{log(delta)}; the default is that of Goren et
#'   al. (2020).
#' @param rr0 Critical value of the rate ratio.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the equal-tailed credible interval.
#' @param details Logical; if \code{FALSE}, the quantiles are not computed
#'   (they are \code{NA}).
#' @param control List of settings of the numerical integration:
#'   \code{n_coarse} (default 41), \code{n_v} (Chebyshev points, default 49),
#'   \code{n_theta} (Gauss-Legendre nodes per panel of the integral over each
#'   control probability, default 12), \code{check} (default \code{TRUE}),
#'   \code{eps} (default 1e-8) and \code{max_iter} (default 12).
#'
#' @return An object of class \code{bssd_post} with model
#'   \code{"MULTIBIN"}, \code{prob}, \code{signal}, \code{rr} (posterior
#'   median and credible interval of \eqn{\delta}{delta}; the mean is
#'   \code{NA}), \code{estimates} (posterior mean of
#'   \eqn{\log \delta}{log(delta)}), labels and settings, \code{data},
#'   \code{prior} and \code{integration}.
#'
#' @references
#' Goren, E., Lin, L.-A., Blaustein, R. O. and Ball, G. (2020). Bayesian
#' meta-analysis of safety outcomes using blinded clinical trial data.
#' \emph{Therapeutic Innovation & Regulatory Science}, 54, 1557--1565.
#' \doi{10.1007/s43441-020-00183-7}
#'
#' @examples
#' # four studies with randomization 1:1 and a control prior with mean 0.02
#' multibin_post(y = c(6, 10, 25, 47), n = c(200, 400, 800, 1600), k = 1,
#'               prior_c = c(2, 98))
#'
#' @export
multibin_post <- function(y, n, k = 1, prior_c, lrr_prior = c(0, 5),
                          rr0 = 1, cutoff = 0.95, level = 0.95,
                          details = TRUE, control = list()) {
  check_num(n, "n", len = NULL, lower = 1, integer = TRUE)
  n_s <- length(n)
  check_num(y, "y", len = n_s, lower = 0, integer = TRUE)
  if (any(y > n)) {
    stop("'y' must not exceed 'n'.", call. = FALSE)
  }
  if (length(k) == 1L) {
    k <- rep(k, n_s)
  }
  check_num(k, "k", len = n_s, lower = 0, lower_open = TRUE)
  if (is.matrix(prior_c)) {
    if (nrow(prior_c) != n_s || ncol(prior_c) != 2L) {
      stop("'prior_c' must be c(a, b) or a matrix with one row per study.",
           call. = FALSE)
    }
  } else {
    check_num(prior_c, "prior_c", len = 2L)
    prior_c <- matrix(prior_c, n_s, 2L, byrow = TRUE)
  }
  check_num(as.vector(prior_c), "prior_c", len = NULL, lower = 0,
            lower_open = TRUE)
  check_num(lrr_prior, "lrr_prior", len = 2L)
  if (lrr_prior[2L] <= 0) {
    stop("The standard deviation in 'lrr_prior' must be positive.",
         call. = FALSE)
  }
  check_common(rr0, cutoff, level, details)
  ctrl <- merge_control(control, list(n_coarse = 41L, n_v = 49L,
                                      n_theta = 12L, check = TRUE,
                                      eps = 1e-8, max_iter = 12L))
  y <- round(y)
  n <- round(n)
  r <- k / (k + 1)
  m <- lrr_prior[[1L]]
  sd <- lrr_prior[[2L]]
  w_of <- function(v) qnorm(pmin(pmax(v, 1e-15), 1 - 1e-15), m, sd)

  logpost <- function(nodes) {
    w <- w_of(nodes[[1L]])
    out <- 0
    for (j in seq_len(n_s)) {
      cc <- 1 - r[j] + r[j] * exp(w)
      out <- out + (if (y[j] > 0) y[j] * log(cc) else 0) +
        multibin_inner(w, y[j], n[j], r[j], prior_c[j, 1L], prior_c[j, 2L],
                       ctrl$n_theta)
    }
    out
  }
  box <- refine_box(logpost, 0, 1, ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
                    0, 1)
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)
  v0 <- pnorm(log(rr0), m, sd)

  fine <- function(n_v, quantiles) {
    g <- post1d(logpost, box, n_v)
    q <- rep(NA_real_, 3L)
    if (quantiles) {
      q <- exp(w_of(inv_cdf(g$cdf, probs, box$lo, box$hi)))
    }
    # the mean of log(delta) by adaptive quadrature of the Chebyshev series
    # of the density, because log(delta) is unbounded at the ends of the
    # prior distribution function
    lrr <- integrate(function(vv) w_of(vv) * g$dens(vv), box$lo, box$hi,
                     rel.tol = 1e-8, subdivisions = 1000L)$value
    list(prob = 1 - g$cdf(v0), q = q, lrr = lrr)
  }
  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_v, quantiles = FALSE)
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_v, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_v <- ctrl$n_v
    res <- fine(n_v, quantiles = details)
    err <- NA_real_
  }

  out <- list(
    model = "MULTIBIN",
    prob = res$prob,
    signal = res$prob > cutoff,
    rr = c(mean = NA_real_, median = res$q[2L], lower = res$q[1L],
           upper = res$q[3L]),
    estimates = c(log_delta = res$lrr),
    incidence = NULL,
    rr_label = "delta",
    rr0 = rr0,
    description = "Tested quantity: rate ratio delta common to the studies",
    cutoff = cutoff,
    level = level,
    t_star = NULL,
    data = list(n = n, n_event = y, k = k),
    prior = list(prior_c = prior_c, lrr_prior = lrr_prior),
    integration = list(method = "Chebyshev and Gauss-Legendre quadrature",
                       nodes = c(v = n_v, theta = ctrl$n_theta),
                       box = c(lower = box$lo, upper = box$hi), error = err)
  )
  class(out) <- "bssd_post"
  out
}
