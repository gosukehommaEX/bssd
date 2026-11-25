#' Posterior of the event proportions under the blinded binomial model
#'
#' Computes the posterior distribution of the proportions of subjects with
#' the event on study drug and control under the Bayesian hierarchical model
#' of Ball and Schnell (2016), from the blinded number of enrolled subjects
#' and the number of subjects with the event.
#'
#' @details
#' With \eqn{N}{N} enrolled subjects of whom \eqn{Z}{Z} had the event, the
#' model is
#' \deqn{M_c \sim \mathrm{Binomial}(N, 1 / (k + 1)), \quad Y_c \mid M_c \sim \mathrm{Binomial}(M_c, \theta_c), \quad Y_d \mid M_c \sim \mathrm{Binomial}(N - M_c, \theta_d),}{M_c ~ Binomial(N, 1 / (k + 1)), Y_c | M_c ~ Binomial(M_c, theta_c), Y_d | M_c ~ Binomial(N - M_c, theta_d),}
#' \deqn{Z = Y_c + Y_d, \quad \theta_c \sim \mathrm{Beta}(a, b), \quad \theta_d \sim \mathrm{Beta}(c, d).}{Z = Y_c + Y_d, theta_c ~ Beta(a, b), theta_d ~ Beta(c, d).}
#' Under simple randomization each subject is on control with probability
#' \eqn{1 / (k + 1)}{1 / (k + 1)}, so
#' \eqn{Z \sim \mathrm{Binomial}(N, \pi)}{Z ~ Binomial(N, pi)} with
#' \eqn{\pi = (\theta_c + k \theta_d) / (k + 1)}{pi = (theta_c + k * theta_d) / (k + 1)},
#' and the posterior is two-dimensional. It is integrated numerically on the
#' coordinates \eqn{\mathrm{logit}(\theta_c)}{logit(theta_c)} and the prior
#' distribution function of \eqn{\theta_d}{theta_d} as in
#' \code{\link{exptime_post}}.
#'
#' When \code{p0} is given, a signal is declared when
#' \eqn{P(\theta_d > p_0 \mid Z) > }{P(theta_d > p0 | Z) > } \code{cutoff}
#' (the rule of Ball and Schnell, 2016); otherwise the ratio
#' \eqn{\theta_d / \theta_c}{theta_d / theta_c} is tested against
#' \code{rr0}. With a single study, the likelihood is that of the model of
#' Goren et al. (2020), see \code{\link{multibin_post}}.
#'
#' @param y Number of subjects with the event (non-negative integer).
#' @param n Number of enrolled subjects (positive integer, at least
#'   \code{y}).
#' @param k Randomization ratio of study drug to control (\code{k:1}).
#' @param prior_c Shape parameters \code{c(a, b)} of the Beta prior of the
#'   control proportion.
#' @param prior_d Shape parameters \code{c(c, d)} of the Beta prior of the
#'   study drug proportion.
#' @param p0 Optional critical value of the study drug proportion (the
#'   maximum acceptable proportion).
#' @param rr0 Critical value of the ratio of the proportions, used for the
#'   signal when \code{p0} is \code{NULL}.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the equal-tailed credible intervals.
#' @param details Logical; if \code{FALSE}, the quantiles are not computed
#'   (they are \code{NA}).
#' @param control List of settings of the numerical integration, as in
#'   \code{\link{exptime_post}}.
#'
#' @return An object of class \code{bssd_post} as described in
#'   \code{\link{exptime_post}}, with model \code{"BINBETA"}, the control and
#'   study drug proportions \code{theta_c} and \code{theta_d}, and their ratio
#'   \code{rr} in place of the rates and the rate ratio.
#'
#' @references
#' Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
#' the FDA IND reporting final rule. In J. Lin et al. (eds.),
#' \emph{Statistical Applications from Clinical Trials and Personalized
#' Medicine to Finance and Business Analytics}, ICSA Book Series in
#' Statistics, 201--209. Springer, Cham.
#' \doi{10.1007/978-3-319-42568-9_17}
#'
#' @seealso \code{\link{signal_count}} for the largest acceptable number of
#'   subjects with the event, \code{\link{canary_post}} for the model with
#'   fixed arm sizes of Gould and Wang (2017).
#'
#' @examples
#' # 40 of 300 subjects with the event, randomization 1:1, control prior from
#' # 30 of 300 historical control subjects, vague prior for study drug
#' fit <- binbeta_post(y = 40, n = 300, k = 1, prior_c = c(31, 271),
#'                     prior_d = c(1, 1), p0 = 0.15)
#' fit
#'
#' @export
binbeta_post <- function(y, n, k = 1, prior_c, prior_d, p0 = NULL, rr0 = 1,
                         cutoff = 0.95, level = 0.95, details = TRUE,
                         control = list()) {
  check_num(n, "n", lower = 1, integer = TRUE)
  check_num(y, "y", lower = 0, upper = n, integer = TRUE)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(prior_c, "prior_c", len = 2L, lower = 0, lower_open = TRUE)
  check_num(prior_d, "prior_d", len = 2L, lower = 0, lower_open = TRUE)
  if (!is.null(p0)) {
    check_num(p0, "p0", lower = 0, upper = 1, lower_open = TRUE,
              upper_open = TRUE)
  }
  check_common(rr0, cutoff, level, details)
  ctrl <- merge_control(control, list(n_coarse = 41L, n_u = 48L, n_v = 49L,
                                      check = TRUE, eps = 1e-8,
                                      max_iter = 12L))
  y <- round(y)
  n <- round(n)
  wc <- 1 / (k + 1)
  ac <- prior_c[[1L]]
  bc <- prior_c[[2L]]
  ad <- prior_d[[1L]]
  bd <- prior_d[[2L]]

  res <- mix2_post(
    loglik = function(tc, td) {
      pm <- outer(wc * tc, (1 - wc) * td, "+")
      ll <- 0
      if (y > 0) {
        ll <- ll + y * log(pm)
      }
      if (n - y > 0) {
        ll <- ll + (n - y) * log1p(-pm)
      }
      ll + matrix(0, length(tc), length(td))
    },
    c_of = plogis,
    lprior_u = function(u) {
      ac * plogis(u, log.p = TRUE) + bc * plogis(-u, log.p = TRUE) -
        lbeta(ac, bc)
    },
    d_of = function(v) qbeta(v, ad, bd),
    F_d = function(x) pbeta(x, ad, bd),
    lo_u = max(qlogis(qbeta(1e-10, ac, bc)), -40),
    hi_u = min(qlogis(qbeta(1 - 1e-10, ac, bc)), 40),
    d_max = 1, rate0 = p0, rr0 = rr0, level = level, details = details,
    ctrl = ctrl)

  mix2_output(
    res, model = "BINBETA", rate0 = p0, rr0 = rr0, cutoff = cutoff,
    level = level,
    labels = c(control = "theta_c", drug = "theta_d", ratio = "rr",
               drug_text = "proportion theta_d of study drug",
               ratio_text = "ratio rr = theta_d / theta_c"),
    data = list(n = n, n_event = y, k = k),
    prior = list(prior_c = prior_c, prior_d = prior_d))
}
