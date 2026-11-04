#' Posterior probability of a safety signal under the EXP_BERN model
#'
#' Computes the posterior distribution of the rate ratio under the
#' exponential Bernoulli (EXP_BERN) mixture model of Duncan (2026), which
#' models whether each subject had the event during the individual follow-up
#' time and does not assume a Poisson distribution for the number of subjects
#' with the event.
#'
#' @details
#' For subject \eqn{i}{i} with follow-up time \eqn{t_i}{t_i}, the model is
#' \deqn{P(Y_i = 1) = w_p [1 - \exp(-\lambda_p t_i)] + w_d [1 - \exp(-r \lambda_p t_i)],}{P(Y_i = 1) = w_p * (1 - exp(-lambda_p * t_i)) + w_d * (1 - exp(-r * lambda_p * t_i)),}
#' with \eqn{w_p = 1 / (k + 1)}{w_p = 1 / (k + 1)} and
#' \eqn{w_d = k / (k + 1)}{w_d = k / (k + 1)}, and
#' \deqn{\lambda_p \sim \mathrm{Gamma}(\lambda_{HP} PT_{HP} / DF, PT_{HP} / DF), \quad r = p / [k (1 - p)], \quad p \sim \mathrm{Beta}(a, b).}{lambda_p ~ Gamma(lambda_hp * pt_hp / df, pt_hp / df), r = p / (k * (1 - p)), p ~ Beta(a, b).}
#' The Gamma distribution is parameterized by shape and rate.
#'
#' The two-dimensional posterior is integrated numerically on the
#' coordinates \eqn{\log \lambda_p}{log(lambda_p)} and the prior
#' distribution function of \eqn{p}{p}. The latter keeps the tails of the
#' posterior of \eqn{p}{p} under control: the likelihood does not vanish as
#' \eqn{r \to 0}{r -> 0} or \eqn{r \to \infty}{r -> Inf}, so these tails
#' follow the prior. A coarse grid first locates the region that holds all but
#' \code{control$eps} of the posterior mass in each marginal. On that region,
#' Gauss-Legendre quadrature is used for \eqn{\log \lambda_p}{log(lambda_p)},
#' and the marginal density of the second coordinate is interpolated at
#' Chebyshev points, which gives its distribution function at any point (used
#' for \code{prob} and the quantiles of \eqn{r}{r}). When
#' \code{control$check} is \code{TRUE}, the computation is repeated with about
#' 1.5 times as many nodes, the finer result is returned, and the absolute
#' difference in \code{prob} is reported as the integration error estimate.
#'
#' Because the likelihood stays positive as \eqn{r \to \infty}{r -> Inf},
#' the posterior means of \eqn{r}{r} and \eqn{r \lambda_p}{r * lambda_p} are
#' infinite when the second shape parameter of the Beta prior is 1 or less,
#' as for the default \code{c(0.5, 0.5)}. They are then returned as
#' \code{Inf}, and the posterior median of \eqn{r}{r} is the summary to use.
#'
#' A signal is declared when \eqn{P(r > rr_0 \mid \mathrm{data}) > }{P(r > rr0 | data) > }
#' \code{cutoff}.
#'
#' @param y Vector of event indicators (1 if the subject had at least one
#'   occurrence of the event during follow-up, 0 otherwise).
#' @param t Vector of follow-up times (positive), in the same time unit as
#'   \code{lambda_hp}, \code{pt_hp} and \code{t_star}.
#' @inheritParams bdribs_post
#' @param t_star Optional time at which the cumulative incidences are
#'   summarized.
#' @param details Logical; if \code{FALSE}, the quantiles of \eqn{r}{r} are
#'   not computed (faster; used in \code{\link{sim_oc}}).
#' @param control List of settings for the numerical integration:
#'   \code{n_coarse} (cells per coordinate of the coarse grid, default 41),
#'   \code{n_lambda} (Gauss-Legendre nodes, default 48), \code{n_p}
#'   (Chebyshev points, default 49), \code{check} (default \code{TRUE}),
#'   \code{eps} (default \code{1e-8}) and \code{max_iter} (default 12).
#'
#' @return An object of class \code{bssd_post}; see \code{\link{bdribs_post}}.
#'   The element \code{model} is \code{"EXP_BERN"}, \code{estimates} contains
#'   the posterior means of \eqn{p}{p}, \eqn{r}{r}, \eqn{\lambda_p}{lambda_p}
#'   and \eqn{r \lambda_p}{r * lambda_p} (\code{lambda_d}), and
#'   \code{integration} contains the nodes used and the error estimate.
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{llbern_post}}
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' t <- runif(n, 10, 52)
#' y <- rbinom(n, 1, 1 - exp(-0.004 * t))
#' fit <- expbern_post(y, t, k = 1, lambda_hp = 0.003, pt_hp = 2000, df = 2,
#'                     t_star = 52)
#' fit
#'
#' @export
expbern_post <- function(y, t, k = 1, lambda_hp, pt_hp, df = 1,
                         p_prior = c(0.5, 0.5), rr0 = 1, cutoff = 0.95,
                         t_star = NULL, level = 0.95, details = TRUE,
                         control = list()) {
  dat <- check_subject_data(y, t)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(lambda_hp, "lambda_hp", lower = 0, lower_open = TRUE)
  check_num(pt_hp, "pt_hp", lower = 0, lower_open = TRUE)
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
  ctrl <- merge_control(control, list(n_coarse = 41L, n_lambda = 48L,
                                      n_p = 49L, check = TRUE, eps = 1e-8,
                                      max_iter = 12L))

  alpha <- lambda_hp * pt_hp / df
  beta <- pt_hp / df
  a <- p_prior[1L]
  b <- p_prior[2L]
  wp <- 1 / (k + 1)

  # log prior density of exp(u) times the Jacobian exp(u)
  lw_lambda <- function(u) {
    alpha * log(beta) - lgamma(alpha) + alpha * u - beta * exp(u)
  }
  # rate ratio from the prior distribution function v of p; 1 - p is
  # computed from the Beta(b, a) distribution for accuracy near p = 1
  r_of <- function(v) {
    qbeta(v, a, b) / (k * pmax(qbeta(1 - v, b, a), 1e-300))
  }
  logpost <- function(nodes) {
    ll <- expbern_loglik_cpp(dat$t, dat$y, wp, exp(nodes[[1L]]),
                             r_of(nodes[[2L]]))
    ll + lw_lambda(nodes[[1L]])
  }

  lo <- c(log(max(qgamma(1e-10, alpha, beta), .Machine$double.xmin)), 0)
  hi <- c(log(qgamma(1 - 1e-10, alpha, beta)), 1)
  box <- refine_box(logpost, lo, hi, ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
                    lim_lo = c(-Inf, 0), lim_hi = c(Inf, 1))

  v0 <- pbeta(k * rr0 / (k * rr0 + 1), a, b)
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  fine <- function(n_lambda, n_p, quantiles) {
    g <- gl_nodes(n_lambda)
    u <- (box$hi[1L] - box$lo[1L]) / 2 * g$x + (box$hi[1L] + box$lo[1L]) / 2
    gw <- g$w * (box$hi[1L] - box$lo[1L]) / 2
    half <- (box$hi[2L] - box$lo[2L]) / 2
    mid <- (box$hi[2L] + box$lo[2L]) / 2
    v <- half * cheb_nodes(n_p) + mid
    L <- logpost(list(u, v)) + log(gw)
    L[is.na(L)] <- -Inf
    W <- exp(L - max(L))
    mv <- colSums(W)
    ccw <- cc_weights(n_p) * half
    norm <- sum(ccw * mv)
    anti <- cheb_antideriv(cheb_coef_matrix(n_p) %*% mv)
    cdf_v <- function(vq) {
      cheb_cumint(anti, (vq - mid) / half) * half / norm
    }
    prob <- min(max(1 - cdf_v(v0), 0), 1)

    wq <- W * rep(ccw, each = n_lambda) / norm
    lam <- exp(u)
    r <- r_of(v)
    lam_d <- outer(lam, r)
    est <- c(p = sum(wq * rep(qbeta(v, a, b), each = n_lambda)),
             r = Inf,
             lambda_p = sum(wq * lam),
             lambda_d = Inf)
    if (b > 1) {
      # r is unbounded as v -> 1 (integrable when b > 1); integrate the
      # Chebyshev series of the marginal densities with stats::integrate,
      # which avoids evaluating r at the end point
      coef_m <- cheb_coef_matrix(n_p) %*% cbind(mv, colSums(W * lam))
      series <- function(vv, j) {
        x <- pmin(pmax((vv - mid) / half, -1), 1)
        drop(cos(outer(acos(x), 0:(n_p - 1L))) %*% coef_m[, j])
      }
      mean_r <- function(j) {
        integrate(function(vv) r_of(vv) * series(vv, j), box$lo[2L],
                  box$hi[2L], rel.tol = 1e-8, subdivisions = 1000L)$value / norm
      }
      est["r"] <- mean_r(1L)
      est["lambda_d"] <- mean_r(2L)
    }
    inc <- NULL
    if (!is.null(t_star)) {
      fp <- -expm1(-lam * t_star)
      fd <- -expm1(-lam_d * t_star)
      inc <- c(placebo = sum(wq * fp), drug = sum(wq * fd),
               difference = sum(wq * (fd - fp)))
    }
    rq <- rep(NA_real_, 3L)
    if (quantiles) {
      q_v <- function(q) {
        tryCatch(uniroot(function(vv) cdf_v(vv) - q, c(box$lo[2L], box$hi[2L]),
                         tol = 1e-12)$root,
                 error = function(e) NA_real_)
      }
      rq <- r_of(vapply(probs, q_v, numeric(1)))
    }
    list(prob = prob, est = est, inc = inc, rq = rq)
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_lambda, ctrl$n_p, quantiles = FALSE)
    n_l <- as.integer(round(1.5 * ctrl$n_lambda))
    n_p <- as.integer(round(1.5 * (ctrl$n_p - 1L))) + 1L
    res <- fine(n_l, n_p, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_l <- ctrl$n_lambda
    n_p <- ctrl$n_p
    res <- fine(n_l, n_p, quantiles = details)
    err <- NA_real_
  }

  out <- list(
    model = "EXP_BERN",
    prob = res$prob,
    signal = res$prob > cutoff,
    rr = c(mean = unname(res$est["r"]), median = res$rq[2L],
           lower = res$rq[1L], upper = res$rq[3L]),
    estimates = res$est,
    incidence = res$inc,
    rr_label = "r",
    rr0 = rr0,
    cutoff = cutoff,
    level = level,
    t_star = t_star,
    data = list(n = length(dat$y), n_event = sum(dat$y),
                exposure = sum(dat$t), k = k),
    prior = list(lambda_hp = lambda_hp, pt_hp = pt_hp, df = df,
                 lambda_shape = alpha, lambda_rate = beta, p_prior = p_prior),
    integration = list(method = "Gauss-Legendre and Chebyshev quadrature",
                       nodes = c(lambda = n_l, p = n_p),
                       box = rbind(lower = box$lo, upper = box$hi),
                       error = err)
  )
  class(out) <- "bssd_post"
  out
}
