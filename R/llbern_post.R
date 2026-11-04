#' Posterior probability of a safety signal under the LL_BERN model
#'
#' Computes the posterior distribution of the relative risk under the
#' log-logistic Bernoulli (LL_BERN) mixture model of Duncan (2026), which
#' allows the hazard of the event to change over time.
#'
#' @details
#' For subject \eqn{i}{i} with follow-up time \eqn{t_i}{t_i}, the model is
#' \deqn{P(Y_i = 1) = w_p F(t_i; \mu_p, \sigma_p) + w_d F(t_i; \mu_d, \sigma_d),}{P(Y_i = 1) = w_p * F(t_i; mu_p, sigma_p) + w_d * F(t_i; mu_d, sigma_d),}
#' where \eqn{F}{F} is the log-logistic cumulative incidence
#' (\code{\link{ll_cdf}}), \eqn{w_p = 1 / (k + 1)}{w_p = 1 / (k + 1)} and
#' \eqn{w_d = k / (k + 1)}{w_d = k / (k + 1)}. The priors are
#' \deqn{\mu_p \sim \mathrm{Gamma}(PT_{HP} \mu / DF_p, PT_{HP} / DF_p), \quad \sigma_p \sim \mathrm{Gamma}(PT_{HP} \sigma / DF_p, PT_{HP} / DF_p),}{mu_p ~ Gamma(pt_hp * mu_hp / df_p, pt_hp / df_p), sigma_p ~ Gamma(pt_hp * sigma_hp / df_p, pt_hp / df_p),}
#' \deqn{\mu_d \sim \mathrm{Gamma}(\mu, 1), \quad \sigma_d \sim \mathrm{Gamma}(\sigma, 1),}{mu_d ~ Gamma(mu_hp, 1), sigma_d ~ Gamma(sigma_hp, 1),}
#' with shape and rate parameters, where \eqn{\mu}{mu} and
#' \eqn{\sigma}{sigma} are the historical placebo values (\code{mu_hp} and
#' \code{sigma_hp}). The priors for the study drug parameters can be changed
#' with \code{mu_d_prior} and \code{sigma_d_prior}.
#'
#' Duncan (2026) does not state which relative risk is used for the signal
#' criterion of this model. This function uses the ratio of the cumulative
#' incidences at the time \code{t_star},
#' \deqn{RR(t^*) = F(t^*; \mu_d, \sigma_d) / F(t^*; \mu_p, \sigma_p),}{RR(t_star) = F(t_star; mu_d, sigma_d) / F(t_star; mu_p, sigma_p),}
#' and declares a signal when
#' \eqn{P[RR(t^*) > rr_0 \mid \mathrm{data}] > }{P(RR(t_star) > rr0 | data) > }
#' \code{cutoff}.
#'
#' The four-dimensional posterior is integrated numerically on the
#' coordinates \eqn{\log \mu_p}{log(mu_p)}, \eqn{\log \sigma_p}{log(sigma_p)},
#' \eqn{\log \mu_d}{log(mu_d)} and the prior distribution function of
#' \eqn{\sigma_d}{sigma_d}; the last one keeps the heavy left tail of the
#' prior of \eqn{\sigma_d}{sigma_d} under control. A coarse grid first locates
#' the region that holds all but \code{control$eps} of the posterior mass in
#' each marginal. On that region, Gauss-Legendre quadrature is used for three
#' coordinates and Chebyshev interpolation for \eqn{\log \mu_d}{log(mu_d)}.
#' Since \eqn{RR(t^*) > rr_0}{RR(t_star) > rr0} holds exactly when
#' \eqn{\mu_d}{mu_d} is below a value determined by the other three
#' parameters, the posterior probability and the quantiles of
#' \eqn{RR(t^*)}{RR(t_star)} are obtained from the cumulative integral of the
#' interpolant without a discontinuous integrand. When \code{control$check}
#' is \code{TRUE}, the computation is repeated with about 1.25 times as many
#' Gauss-Legendre nodes and 1.5 times as many Chebyshev points, the finer
#' result is returned, and the absolute difference in \code{prob} is reported
#' as the integration error estimate.
#'
#' When there are few events, the prior of \eqn{\sigma_d}{sigma_d} used by
#' Duncan (2026) puts much mass near zero, where the likelihood changes
#' sharply in \eqn{\mu_d}{mu_d} near the logarithms of the follow-up times.
#' The default number of Chebyshev points is large for this reason; if the
#' reported error estimate is too large for the purpose, increase
#' \code{control$n_inner}.
#'
#' The posterior mean of \eqn{RR(t^*)}{RR(t_star)} can be unstable because
#' its posterior distribution may have a heavy right tail; the posterior
#' median is a more stable summary.
#'
#' @param y Vector of event indicators (1 if the subject had at least one
#'   occurrence of the event during follow-up, 0 otherwise).
#' @param t Vector of follow-up times (positive), in the time unit for which
#'   \code{mu_hp} and \code{sigma_hp} were estimated.
#' @param k Randomization ratio of study drug to placebo (\code{k:1}).
#' @param mu_hp Location parameter \eqn{\mu}{mu} of the log-logistic
#'   distribution estimated from historical placebo data (positive).
#' @param sigma_hp Scale parameter \eqn{\sigma}{sigma} of the log-logistic
#'   distribution estimated from historical placebo data (positive).
#' @param pt_hp Total patient time of the historical placebo population.
#' @param df_p Discount factor for the historical placebo data.
#' @param mu_d_prior Shape and rate of the Gamma prior for \eqn{\mu_d}{mu_d};
#'   the default is \code{c(mu_hp, 1)}.
#' @param sigma_d_prior Shape and rate of the Gamma prior for
#'   \eqn{\sigma_d}{sigma_d}; the default is \code{c(sigma_hp, 1)}.
#' @param t_star Time at which the relative risk and the cumulative
#'   incidences are evaluated.
#' @param rr0 Critical value of the relative risk.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the equal-tailed credible interval of the
#'   relative risk.
#' @param details Logical; if \code{FALSE}, the quantiles of
#'   \eqn{RR(t^*)}{RR(t_star)} are not computed (faster; used in
#'   \code{\link{sim_oc}}).
#' @param control List of settings for the numerical integration:
#'   \code{n_coarse} (cells per coordinate of the coarse grid, default 12),
#'   \code{n_outer} (Gauss-Legendre nodes per coordinate, default 16),
#'   \code{n_inner} (Chebyshev points, default 121), \code{check} (default
#'   \code{TRUE}), \code{eps} (default \code{1e-8}) and \code{max_iter}
#'   (default 12).
#'
#' @return An object of class \code{bssd_post}; see \code{\link{bdribs_post}}.
#'   The element \code{model} is \code{"LL_BERN"}, \code{rr} summarizes
#'   \eqn{RR(t^*)}{RR(t_star)}, \code{estimates} contains the posterior means
#'   of \eqn{\mu_p}{mu_p}, \eqn{\sigma_p}{sigma_p}, \eqn{\mu_d}{mu_d} and
#'   \eqn{\sigma_d}{sigma_d}, and \code{integration} contains the nodes used
#'   and the error estimate.
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{ll_cdf}}, \code{\link{ll_mu}},
#'   \code{\link{expbern_post}}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' t <- runif(n, 4, 104)
#' y <- rbinom(n, 1, ll_cdf(t, mu = 6.5, sigma = 0.8))
#' fit <- llbern_post(y, t, k = 1, mu_hp = 6.5, sigma_hp = 0.8, pt_hp = 200,
#'                    t_star = 104, control = list(check = FALSE))
#' fit
#'
#' @export
llbern_post <- function(y, t, k = 1, mu_hp, sigma_hp, pt_hp, df_p = 1,
                        mu_d_prior = c(mu_hp, 1),
                        sigma_d_prior = c(sigma_hp, 1), t_star, rr0 = 1,
                        cutoff = 0.95, level = 0.95, details = TRUE,
                        control = list()) {
  dat <- check_subject_data(y, t)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(mu_hp, "mu_hp", lower = 0, lower_open = TRUE)
  check_num(sigma_hp, "sigma_hp", lower = 0, lower_open = TRUE)
  check_num(pt_hp, "pt_hp", lower = 0, lower_open = TRUE)
  check_num(df_p, "df_p", lower = 0, lower_open = TRUE)
  check_num(mu_d_prior, "mu_d_prior", len = 2L, lower = 0, lower_open = TRUE)
  check_num(sigma_d_prior, "sigma_d_prior", len = 2L, lower = 0,
            lower_open = TRUE)
  if (missing(t_star)) {
    stop("'t_star' must be specified.", call. = FALSE)
  }
  check_num(t_star, "t_star", lower = 0, lower_open = TRUE)
  check_num(rr0, "rr0", lower = 0, lower_open = TRUE)
  check_num(cutoff, "cutoff", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  check_num(level, "level", lower = 0, upper = 1,
            lower_open = TRUE, upper_open = TRUE)
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("'details' must be TRUE or FALSE.", call. = FALSE)
  }
  ctrl <- merge_control(control, list(n_coarse = 12L, n_outer = 16L,
                                      n_inner = 121L, check = TRUE,
                                      eps = 1e-8, max_iter = 12L))

  wp <- 1 / (k + 1)
  shape <- c(pt_hp * mu_hp / df_p, pt_hp * sigma_hp / df_p,
             mu_d_prior[1L], sigma_d_prior[1L])
  rate <- c(pt_hp / df_p, pt_hp / df_p, mu_d_prior[2L], sigma_d_prior[2L])
  lt_star <- log(t_star)

  # log prior density of exp(u) times the Jacobian exp(u)
  lw_log <- function(u, j) {
    shape[j] * log(rate[j]) - lgamma(shape[j]) + shape[j] * u - rate[j] * exp(u)
  }
  # sigma_d from its prior distribution function
  sd_of <- function(v) {
    pmax(qgamma(v, shape = shape[4L], rate = rate[4L]), 1e-300)
  }
  logpost <- function(nodes) {
    n <- lengths(nodes)
    pp <- expand.grid(mp = exp(nodes[[1L]]), sp = exp(nodes[[2L]]))
    dd <- expand.grid(md = exp(nodes[[3L]]), sd = sd_of(nodes[[4L]]))
    ll <- llbern_loglik_cpp(dat$t, dat$y, wp, pp$mp, pp$sp, dd$md, dd$sd)
    lp <- outer(outer(lw_log(nodes[[1L]], 1L), lw_log(nodes[[2L]], 2L), "+"),
                lw_log(nodes[[3L]], 3L), "+")
    lp <- outer(lp, numeric(n[4L]), "+")
    array(ll, n) + lp
  }

  q_log <- function(p, j) {
    log(max(qgamma(p, shape[j], rate[j]), .Machine$double.xmin))
  }
  lo <- c(q_log(1e-10, 1L), q_log(1e-10, 2L), q_log(1e-10, 3L), 0)
  hi <- c(q_log(1 - 1e-10, 1L), q_log(1 - 1e-10, 2L), q_log(1 - 1e-10, 3L), 1)
  box <- refine_box(logpost, lo, hi, ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
                    lim_lo = c(-Inf, -Inf, -Inf, 0),
                    lim_hi = c(Inf, Inf, Inf, 1))
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  fine <- function(n_outer, n_inner, quantiles) {
    g <- gl_nodes(n_outer)
    map_gl <- function(j) {
      list(x = (box$hi[j] - box$lo[j]) / 2 * g$x + (box$hi[j] + box$lo[j]) / 2,
           w = g$w * (box$hi[j] - box$lo[j]) / 2)
    }
    d1 <- map_gl(1L)
    d2 <- map_gl(2L)
    d4 <- map_gl(4L)
    half <- (box$hi[3L] - box$lo[3L]) / 2
    mid <- (box$hi[3L] + box$lo[3L]) / 2
    u3 <- half * cheb_nodes(n_inner) + mid
    L <- logpost(list(d1$x, d2$x, u3, d4$x))
    lgw <- outer(outer(outer(log(d1$w), log(d2$w), "+"), numeric(n_inner),
                       "+"), log(d4$w), "+")
    L <- L + lgw
    L[is.na(L)] <- -Inf
    W <- exp(L - max(L))
    # rows: mu_d nodes; columns: (mu_p, sigma_p, sigma_d) with mu_p fastest
    Wm <- matrix(aperm(W, c(3L, 1L, 2L, 4L)), nrow = n_inner)
    ccw <- cc_weights(n_inner) * half
    norm <- sum(ccw * Wm)
    anti <- cheb_antideriv(cheb_coef_matrix(n_inner) %*% Wm)
    tr <- expand.grid(mp = exp(d1$x), sp = exp(d2$x), sd = sd_of(d4$x))
    fp_star <- plogis((lt_star - tr$mp) / tr$sp)

    # P(RR(t_star) > q | data): RR > q if and only if mu_d < mu_star
    p_exceed <- function(q) {
      arg <- q * fp_star
      mstar <- lt_star - tr$sd * qlogis(pmin(arg, 1))
      valid <- arg < 1 & mstar > 0
      xs <- rep(-1, length(arg))
      xs[valid] <- (log(mstar[valid]) - mid) / half
      ci <- cheb_cumint(anti, xs) * half
      ci[!valid] <- 0
      min(max(sum(ci) / norm, 0), 1)
    }
    prob <- p_exceed(rr0)

    n_tr <- nrow(tr)
    wq <- Wm * ccw / norm
    md <- exp(u3)
    fd <- plogis((lt_star - matrix(md, n_inner, n_tr)) /
                   matrix(tr$sd, n_inner, n_tr, byrow = TRUE))
    fp <- matrix(fp_star, n_inner, n_tr, byrow = TRUE)
    rr_g <- fd / pmax(fp, 1e-300)
    est <- c(mu_p = sum(wq * matrix(tr$mp, n_inner, n_tr, byrow = TRUE)),
             sigma_p = sum(wq * matrix(tr$sp, n_inner, n_tr, byrow = TRUE)),
             mu_d = sum(wq * md),
             sigma_d = sum(wq * matrix(tr$sd, n_inner, n_tr, byrow = TRUE)))
    inc <- c(placebo = sum(wq * fp), drug = sum(wq * fd),
             difference = sum(wq * (fd - fp)))
    rr_mean <- sum(wq * rr_g)
    rq <- rep(NA_real_, 3L)
    if (quantiles) {
      # starting brackets from the weighted grid distribution of RR(t_star)
      ord <- order(rr_g)
      cw <- cumsum(wq[ord])
      q_rr <- function(level_q) {
        start <- rr_g[ord][min(which(cw >= level_q), length(ord))]
        rng <- log(max(start, 1e-300)) + c(-0.5, 0.5)
        tryCatch(exp(uniroot(function(lq) p_exceed(exp(lq)) - (1 - level_q),
                             rng, extendInt = "downX", tol = 1e-8)$root),
                 error = function(e) NA_real_)
      }
      rq <- vapply(probs, q_rr, numeric(1))
    }
    list(prob = prob, est = est, inc = inc, rr_mean = rr_mean, rq = rq)
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_outer, ctrl$n_inner, quantiles = FALSE)
    n_o <- as.integer(round(1.25 * ctrl$n_outer))
    n_i <- as.integer(round(1.5 * (ctrl$n_inner - 1L))) + 1L
    res <- fine(n_o, n_i, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_o <- ctrl$n_outer
    n_i <- ctrl$n_inner
    res <- fine(n_o, n_i, quantiles = details)
    err <- NA_real_
  }

  out <- list(
    model = "LL_BERN",
    prob = res$prob,
    signal = res$prob > cutoff,
    rr = c(mean = res$rr_mean, median = res$rq[2L], lower = res$rq[1L],
           upper = res$rq[3L]),
    estimates = res$est,
    incidence = res$inc,
    rr_label = "RR(t_star)",
    rr0 = rr0,
    cutoff = cutoff,
    level = level,
    t_star = t_star,
    data = list(n = length(dat$y), n_event = sum(dat$y),
                exposure = sum(dat$t), k = k),
    prior = list(mu_hp = mu_hp, sigma_hp = sigma_hp, pt_hp = pt_hp,
                 df_p = df_p, shape = shape, rate = rate),
    integration = list(method = "Gauss-Legendre and Chebyshev quadrature",
                       nodes = c(outer = n_o, inner = n_i),
                       box = rbind(lower = box$lo, upper = box$hi),
                       error = err)
  )
  class(out) <- "bssd_post"
  out
}
