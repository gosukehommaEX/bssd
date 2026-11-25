#' Posterior probability of a higher event probability on study drug with a
#' meta-analytic predictive prior
#'
#' Computes the posterior probability that the event probability per unit
#' time on study drug exceeds that on control, from blinded subject-level
#' data, with the mixture likelihood and the meta-analytic predictive (MAP)
#' prior of Brock et al. (2023).
#'
#' @details
#' Let \eqn{\pi_c}{pi_c} and \eqn{\pi_d}{pi_d} be the probabilities of at
#' least one event in one unit of time on control and study drug. With a
#' constant hazard, the probability of an event of subject \eqn{i}{i} with
#' time at risk \eqn{t_i}{t_i} is \eqn{1 - (1 - \pi)^{t_i}}{1 - (1 - pi)^t_i},
#' and the likelihood of the blinded data is the mixture
#' \deqn{\prod_i [p_B(y_i \mid 1 - (1 - \pi_c)^{t_i}) / (k + 1) + k p_B(y_i \mid 1 - (1 - \pi_d)^{t_i}) / (k + 1)],}{prod_i (pB(y_i | 1 - (1 - pi_c)^t_i) / (k + 1) + k * pB(y_i | 1 - (1 - pi_d)^t_i) / (k + 1)),}
#' where \eqn{p_B}{pB} is the Bernoulli probability. The prior of
#' \eqn{\pi_c}{pi_c} is a mixture of Beta distributions, typically the MAP
#' prior from historical control data (\code{\link{map_beta}}); the prior of
#' \eqn{\pi_d}{pi_d} is a Beta distribution, by default the vague
#' Beta(1/3, 1/3) used by Brock et al. (2023). The quantity of interest is
#' \eqn{\gamma = P(\pi_d > rr_0 \pi_c \mid \mathrm{data})}{gamma = P(pi_d > rr0 * pi_c | data)};
#' Brock et al. (2023) use \eqn{rr_0 = 1}{rr0 = 1} and mention larger values.
#'
#' The likelihood is that of the exponential Bernoulli mixture model
#' (\code{\link{expbern_post}}) with hazards
#' \eqn{-\log(1 - \pi)}{-log(1 - pi)}; only the priors differ. The posterior
#' is integrated numerically on the coordinates
#' \eqn{\mathrm{logit}(\pi_c)}{logit(pi_c)} and the prior distribution
#' function of \eqn{\pi_d}{pi_d}, as in \code{\link{exptime_post}} and
#' \code{\link{binbeta_post}}.
#'
#' A prior derived with the \pkg{RBesT} package can be used: the rows
#' \code{w}, \code{a} and \code{b} of a \code{mixbeta} object are the
#' columns \code{w}, \code{shape1} and \code{shape2} of \code{prior_c}.
#'
#' @param y Event indicators of the subjects (1 if the subject had at least
#'   one event).
#' @param t Times at risk of the subjects (positive), in the unit of time of
#'   the probabilities.
#' @param k Randomization ratio of study drug to control (\code{k:1}); Brock
#'   et al. (2023) use \code{k = 3}.
#' @param prior_c Prior of the control probability: \code{c(shape1, shape2)}
#'   of a Beta distribution, or a data frame with columns \code{w},
#'   \code{shape1} and \code{shape2} of a mixture of Beta distributions, such
#'   as the result of \code{\link{map_beta}}.
#' @param prior_d Shape parameters of the Beta prior of the study drug
#'   probability.
#' @param rate0 Optional critical value of the study drug probability; if
#'   given, it defines the signal instead of \code{rr0}.
#' @param rr0 Critical value of the ratio \eqn{\pi_d / \pi_c}{pi_d / pi_c}.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the equal-tailed credible intervals.
#' @param details Logical; if \code{FALSE}, the quantiles are not computed
#'   (they are \code{NA}).
#' @param control List of settings of the numerical integration, as in
#'   \code{\link{exptime_post}}.
#'
#' @return An object of class \code{bssd_post} as described in
#'   \code{\link{exptime_post}}, with model \code{"MAPBERN"}, the
#'   probabilities \code{pi_c} and \code{pi_d}, and their ratio \code{rr}.
#'   \code{prob} is \eqn{\gamma}{gamma}.
#'
#' @references
#' Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
#' Penard, N. (2023). A Bayesian method for safety signal detection in
#' ongoing blinded randomised controlled trials. \emph{Pharmaceutical
#' Statistics}, 22(2), 378--395. \doi{10.1002/pst.2278}
#'
#' @seealso \code{\link{map_beta}}
#'
#' @examples
#' # MAP prior from three historical control arms of 100 subjects each
#' prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
#' set.seed(2)
#' n <- 200
#' drug <- rbinom(n, 1, 0.75)
#' y <- rbinom(n, 1, ifelse(drug == 1, 0.18, 0.095))
#' mapbern_post(y, t = rep(1, n), k = 3, prior_c = prior)
#'
#' @export
mapbern_post <- function(y, t, k = 1, prior_c, prior_d = c(1 / 3, 1 / 3),
                         rate0 = NULL, rr0 = 1, cutoff = 0.95, level = 0.95,
                         details = TRUE, control = list()) {
  dat <- check_subject_data(y, t)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  mix <- as_mixbeta(prior_c, "prior_c")
  check_num(prior_d, "prior_d", len = 2L, lower = 0, lower_open = TRUE)
  if (!is.null(rate0)) {
    check_num(rate0, "rate0", lower = 0, upper = 1, lower_open = TRUE,
              upper_open = TRUE)
  }
  check_common(rr0, cutoff, level, details)
  ctrl <- merge_control(control, list(n_coarse = 41L, n_u = 48L, n_v = 49L,
                                      check = TRUE, eps = 1e-8,
                                      max_iter = 12L))
  pat <- group_patterns(dat$y, dat$t)
  wc <- 1 / (k + 1)
  ad <- prior_d[[1L]]
  bd <- prior_d[[2L]]
  lw <- log(mix$w)
  lb <- lbeta(mix$shape1, mix$shape2)

  # log prior density of logit(pi_c) (mixture times the Jacobian)
  lprior_u <- function(u) {
    l1 <- plogis(u, log.p = TRUE)
    l0 <- plogis(-u, log.p = TRUE)
    L <- outer(l1, mix$shape1) + outer(l0, mix$shape2) +
      rep(lw - lb, each = length(u))
    m <- apply(L, 1L, max)
    m + log(rowSums(exp(L - m)))
  }
  lo_u <- max(min(qlogis(qbeta(1e-10, mix$shape1, mix$shape2))), -40)
  hi_u <- min(max(qlogis(qbeta(1 - 1e-10, mix$shape1, mix$shape2))), 40)

  res <- mix2_post(
    loglik = function(pc, PD) {
      # one row of study drug probabilities per control probability
      hc <- -log1p(-pc)
      HD <- -log1p(-PD)
      out <- matrix(0, nrow(PD), ncol(PD))
      for (i in seq_len(nrow(PD))) {
        out[i, ] <- bernmix_loglik_cpp(pat$t, pat$y, pat$w, wc, hc[i],
                                       HD[i, ])
      }
      out
    },
    c_of = plogis, u_of = qlogis, lprior_u = lprior_u,
    pfun = function(x, lower.tail = TRUE, log.p = FALSE) {
      pbeta(x, ad, bd, lower.tail = lower.tail, log.p = log.p)
    },
    qfun = function(p, lower.tail = TRUE, log.p = FALSE) {
      qbeta(p, ad, bd, lower.tail = lower.tail, log.p = log.p)
    },
    lo_u = lo_u, hi_u = hi_u, d_max = 1,
    u_split = if (rr0 > 1) qlogis(1 / rr0) else numeric(0), rate0 = rate0,
    rr0 = rr0, level = level, details = details, ctrl = ctrl)

  mix2_output(
    res, model = "MAPBERN", rate0 = rate0, rr0 = rr0, cutoff = cutoff,
    level = level,
    labels = c(control = "pi_c", drug = "pi_d", ratio = "rr",
               drug_text = "event probability pi_d of study drug",
               ratio_text = "ratio rr = pi_d / pi_c"),
    data = list(n = length(dat$y), n_event = sum(dat$y),
                exposure = sum(dat$t), k = k),
    prior = list(prior_c = mix, prior_d = prior_d))
}
