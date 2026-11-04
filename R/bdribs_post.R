#' Posterior probability of a safety signal under the BDRIBS model
#'
#' Computes the posterior distribution of the relative risk under the
#' Poisson-based BDRIBS model of Mukhopadhyay et al. (2018) as specified in
#' Duncan (2026). With the total exposure time as \code{exposure}, the model
#' is BDRIBS(1); with the total time at risk (time after the first occurrence
#' of the event removed), it is BDRIBS(2). Both can be computed from
#' subject-level data with \code{\link{time_at_risk}}.
#'
#' @details
#' The model is
#' \deqn{Y \sim \mathrm{Poisson}[E \delta (r k + 1) / (k + 1)],}{Y ~ Poisson(E * delta * (r * k + 1) / (k + 1)),}
#' \deqn{\delta \sim \mathrm{Gamma}(\lambda_{HP} PT_{HP} / DF, PT_{HP} / DF),}{delta ~ Gamma(lambda_hp * pt_hp / df, pt_hp / df),}
#' \deqn{r = p / [k (1 - p)], \quad p \sim \mathrm{Beta}(a, b),}{r = p / (k * (1 - p)), p ~ Beta(a, b),}
#' where \eqn{Y}{Y} is the number of subjects with the event, \eqn{E}{E} is
#' the exposure, \eqn{k}{k} is the randomization ratio (study drug to placebo
#' \eqn{k:1}{k:1}), \eqn{\delta}{delta} is the event rate in the placebo
#' group, and \eqn{r}{r} is the rate ratio of study drug to placebo.
#' The Gamma distribution is parameterized by shape and rate.
#'
#' Because the Gamma prior is conjugate to the Poisson likelihood,
#' \eqn{\delta}{delta} is integrated out analytically, and the marginal
#' posterior density of \eqn{p}{p} is proportional to
#' \deqn{p^{a - 1} (1 - p)^{b + \alpha - 1} [E + \beta (k + 1) (1 - p)]^{-(\alpha + Y)}}{p^(a - 1) * (1 - p)^(b + alpha - 1) * (E + beta * (k + 1) * (1 - p))^(-(alpha + Y))}
#' with \eqn{\alpha = \lambda_{HP} PT_{HP} / DF}{alpha = lambda_hp * pt_hp / df}
#' and \eqn{\beta = PT_{HP} / DF}{beta = pt_hp / df}. All posterior
#' quantities are one-dimensional integrals over \eqn{\mathrm{logit}(p)}{logit(p)},
#' which are evaluated by adaptive quadrature (\code{\link[stats]{integrate}})
#' with a relative tolerance of \code{1e-10}.
#'
#' A signal is declared when \eqn{P(r > rr_0 \mid Y) > }{P(r > rr0 | Y) > }
#' \code{cutoff} (Duncan, 2026, Equation 3).
#'
#' When \code{t_star} is given, the posterior means of the cumulative
#' incidences at \code{t_star} are computed under the exponential
#' distribution implied by the constant event rates, that is,
#' \eqn{1 - \exp(-\delta t^*)}{1 - exp(-delta * t_star)} for placebo and
#' \eqn{1 - \exp(-r \delta t^*)}{1 - exp(-r * delta * t_star)} for study drug.
#'
#' With \code{pt_hp = Inf}, the prior of \eqn{\delta}{delta} is a point mass at
#' \code{lambda_hp}, so the inference is conditional on a known background
#' rate. This is the screening analysis of Mukhopadhyay et al. (2018) and
#' corresponds to \code{inf.type = 1} of the \pkg{bdribs} package; the marginal
#' posterior density of \eqn{p}{p} is then proportional to
#' \deqn{p^{a - 1} (1 - p)^{b - Y - 1} \exp[-E \lambda_{HP} / ((k + 1) (1 - p))].}{p^(a - 1) * (1 - p)^(b - Y - 1) * exp(-E * lambda_hp / ((k + 1) * (1 - p))).}
#'
#' The time unit of \code{exposure}, \code{lambda_hp}, \code{pt_hp} and
#' \code{t_star} must be the same.
#'
#' @param y Number of subjects with the event in the blinded data
#'   (non-negative integer).
#' @param exposure Total exposure time (BDRIBS(1)) or total time at risk
#'   (BDRIBS(2)) of the blinded population (positive).
#' @param k Randomization ratio of study drug to placebo (\code{k:1}).
#' @param lambda_hp Event rate in the historical placebo population.
#' @param pt_hp Total patient time at risk in the historical placebo
#'   population, or \code{Inf} for a known background rate
#'   \code{lambda_hp}.
#' @param df Discount factor for the historical data (\code{df >= 1} discounts
#'   the historical information).
#' @param p_prior Shape parameters \code{c(a, b)} of the Beta prior for
#'   \eqn{p}{p}. The default \code{c(0.5, 0.5)} is the prior used by Duncan
#'   (2026).
#' @param rr0 Critical value of the rate ratio (\eqn{\theta_c}{theta_c} in
#'   Duncan, 2026).
#' @param cutoff Posterior probability threshold for a signal.
#' @param t_star Optional time at which the cumulative incidences are
#'   summarized.
#' @param level Probability of the equal-tailed credible interval of the rate
#'   ratio.
#' @param details Logical; if \code{FALSE}, only \code{prob}, \code{signal}
#'   and \code{incidence} are computed, and \code{rr} and \code{estimates}
#'   are \code{NA}. This is faster and is used in \code{\link{sim_oc}}.
#'
#' @return An object of class \code{bssd_post}, a list with elements
#' \describe{
#'   \item{\code{model}}{\code{"BDRIBS"}.}
#'   \item{\code{prob}}{Posterior probability \eqn{P(r > rr_0 \mid Y)}{P(r > rr0 | Y)}.}
#'   \item{\code{signal}}{Logical; \code{prob > cutoff}.}
#'   \item{\code{rr}}{Posterior mean, median, and the limits of the
#'     equal-tailed credible interval of \eqn{r}{r}.}
#'   \item{\code{estimates}}{Posterior means of \eqn{p}{p}, \eqn{r}{r}, the
#'     placebo event rate \eqn{\delta}{delta} (\code{lambda_p}) and the study
#'     drug event rate \eqn{r \delta}{r * delta} (\code{lambda_d}).}
#'   \item{\code{incidence}}{Posterior means of the cumulative incidences at
#'     \code{t_star} for placebo and study drug and of their difference, or
#'     \code{NULL} when \code{t_star} is not given.}
#'   \item{\code{rr_label}, \code{rr0}, \code{cutoff}, \code{level},
#'     \code{t_star}}{Settings used.}
#'   \item{\code{data}, \code{prior}}{Data summary and prior parameters.}
#'   \item{\code{integration}}{Description of the numerical integration.}
#' }
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
#' detection of potential risk using inference on blinded safety data.
#' \emph{Pharmaceutical Statistics}, 17(6), 823--834. \doi{10.1002/pst.1898}
#'
#' @seealso \code{\link{time_at_risk}}, \code{\link{expbern_post}},
#'   \code{\link{llbern_post}}
#'
#' @examples
#' # 5 subjects with the event in 500 patient-years, randomization 2:1,
#' # historical placebo rate 18 / 4000 per patient-year
#' fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
#'                    pt_hp = 4000)
#' fit
#'
#' # known background rate of 0.011 per patient-year and uniform prior for p
#' # (hypothetical example 2 of Mukhopadhyay et al., 2018)
#' bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011, pt_hp = Inf,
#'             p_prior = c(1, 1))
#'
#' @export
bdribs_post <- function(y, exposure, k = 1, lambda_hp, pt_hp, df = 1,
                        p_prior = c(0.5, 0.5), rr0 = 1, cutoff = 0.95,
                        t_star = NULL, level = 0.95, details = TRUE) {
  check_num(y, "y", lower = 0, integer = TRUE)
  check_num(exposure, "exposure", lower = 0, lower_open = TRUE)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(lambda_hp, "lambda_hp", lower = 0, lower_open = TRUE)
  fixed <- is.numeric(pt_hp) && length(pt_hp) == 1L && isTRUE(pt_hp == Inf)
  if (!fixed) {
    check_num(pt_hp, "pt_hp", lower = 0, lower_open = TRUE)
  }
  check_num(df, "df", lower = 0, lower_open = TRUE)
  check_num(p_prior, "p_prior", len = 2L, lower = 0, lower_open = TRUE)
  check_num(rr0, "rr0", lower = 0, lower_open = TRUE)
  check_num(cutoff, "cutoff", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  check_num(level, "level", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  if (!is.null(t_star)) {
    check_num(t_star, "t_star", lower = 0, lower_open = TRUE)
  }
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("'details' must be TRUE or FALSE.", call. = FALSE)
  }
  y <- round(y)
  alpha <- lambda_hp * pt_hp / df
  beta <- pt_hp / df
  a <- p_prior[1L]
  b <- p_prior[2L]
  shape_post <- alpha + y

  # z = logit(p) = log(k * r); log marginal posterior density of z
  if (fixed) {
    rate_hp <- exposure * lambda_hp / (k + 1)
    lk <- function(z) {
      a * plogis(z, log.p = TRUE) + (b - y) * plogis(-z, log.p = TRUE) -
        rate_hp * (1 + exp(z))
    }
  } else {
    lk <- function(z) {
      a * plogis(z, log.p = TRUE) + (b + alpha) * plogis(-z, log.p = TRUE) -
        shape_post * log(exposure + beta * (k + 1) * plogis(-z))
    }
  }
  # posterior rate of delta given z is beta + exposure * (1 + exp(z)) / (k + 1)
  inv_rate <- function(z) {
    1 / (beta + exposure * (1 + exp(z)) / (k + 1))
  }
  r_inv_rate <- function(z) {
    (1 / k) / (beta * exp(-z) + exposure * (exp(-z) + 1) / (k + 1))
  }

  z0 <- log(k * rr0)
  l_up <- log_integrate_exp(lk, z0, Inf)
  l_lo <- log_integrate_exp(lk, -Inf, z0)
  prob <- 1 / (1 + exp(l_lo - l_up))
  log_z <- log_integrate_exp(lk)

  post_mean <- function(lfun) exp(log_integrate_exp(lfun) - log_z)

  rq <- rep(NA_real_, 3L)
  r_mean <- NA_real_
  est <- c(p = NA_real_, r = NA_real_, lambda_p = NA_real_,
           lambda_d = NA_real_)
  if (details) {
    cdf_z <- function(z) exp(log_integrate_exp(lk, -Inf, z) - log_z)
    q_z <- function(q) {
      tryCatch(uniroot(function(z) cdf_z(z) - q, c(-50, 50),
                       extendInt = "upX", tol = 1e-10)$root,
               error = function(e) NA_real_)
    }
    probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)
    rq <- exp(vapply(probs, q_z, numeric(1))) / k
    r_mean <- if (fixed || b + alpha > 1) {
      post_mean(function(z) lk(z) + z) / k
    } else {
      Inf
    }
    p_mean <- post_mean(function(z) lk(z) + plogis(z, log.p = TRUE))
    if (fixed) {
      est <- c(p = p_mean, r = r_mean, lambda_p = lambda_hp,
               lambda_d = r_mean * lambda_hp)
    } else {
      est <- c(p = p_mean, r = r_mean,
               lambda_p = post_mean(function(z) lk(z) + log(shape_post) +
                                      log(inv_rate(z))),
               lambda_d = post_mean(function(z) lk(z) + log(shape_post) +
                                      log(r_inv_rate(z))))
    }
  }

  inc <- NULL
  if (!is.null(t_star)) {
    if (fixed) {
      surv_p <- exp(-lambda_hp * t_star)
      surv_d <- post_mean(function(z) lk(z) - exp(z) / k * lambda_hp * t_star)
    } else {
      # Gamma moment generating function: E(exp(-s * delta)) given z
      surv_p <- post_mean(function(z) lk(z) - shape_post *
                            log1p(t_star * inv_rate(z)))
      surv_d <- post_mean(function(z) lk(z) - shape_post *
                            log1p(t_star * r_inv_rate(z)))
    }
    inc <- c(placebo = 1 - surv_p, drug = 1 - surv_d,
             difference = surv_p - surv_d)
  }

  out <- list(
    model = "BDRIBS",
    prob = prob,
    signal = prob > cutoff,
    rr = c(mean = unname(r_mean), median = rq[2L], lower = rq[1L],
           upper = rq[3L]),
    estimates = est,
    incidence = inc,
    rr_label = "r",
    rr0 = rr0,
    cutoff = cutoff,
    level = level,
    t_star = t_star,
    data = list(y = y, exposure = exposure, k = k),
    prior = list(lambda_hp = lambda_hp, pt_hp = pt_hp, df = df,
                 delta_shape = alpha, delta_rate = beta, p_prior = p_prior),
    integration = list(method = "adaptive quadrature (stats::integrate)",
                       error = NA_real_)
  )
  class(out) <- "bssd_post"
  out
}
