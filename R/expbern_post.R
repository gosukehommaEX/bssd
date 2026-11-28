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
#' The two-dimensional posterior is integrated numerically. The first
#' coordinate is \eqn{\log \lambda_p}{log(lambda_p)}, integrated by
#' Gauss-Legendre quadrature on the region that a coarse grid locates. The
#' second coordinate describes \eqn{p}{p} through its prior distribution
#' function \eqn{F}{F}: either \eqn{v = F(p)}{v = F(p)}, in which the
#' posterior density is the likelihood times the prior density of
#' \eqn{\lambda_p}{lambda_p}, or
#' \eqn{z = \Phi^{-1}(F(p))}{z = qnorm(F(p))}, computed on the logarithmic
#' scale, which also represents the far tails of an informative prior of
#' \eqn{p}{p}. The likelihood does not vanish as \eqn{r \to 0}{r -> 0} or
#' \eqn{r \to \infty}{r -> Inf}, so the tails of the posterior of
#' \eqn{p}{p} follow the prior. For each Gauss-Legendre node, the interval of
#' \eqn{z}{z} that holds all but \code{control$eps} of the conditional
#' posterior mass on each side is located separately, so the intervals follow
#' the posterior when it lies on a ridge of the blinded likelihood or in a
#' tail of the prior of \eqn{p}{p} (for example, when the data conflict with
#' an informative prior of \eqn{p}{p}). Each row is interpolated at Chebyshev
#' points in \eqn{v}{v} or in \eqn{z}{z}, whichever gives the faster decay of
#' the Chebyshev coefficients of its density, which gives the distribution
#' function of \eqn{p}{p}, and hence of \eqn{r}{r}, at any point (used for
#' \code{prob} and the quantiles of \eqn{r}{r}). When \code{control$check}
#' is \code{TRUE}, the computation is repeated with about 1.5 times as many
#' nodes, the finer result is returned, and the absolute difference in
#' \code{prob} is reported as the integration error estimate.
#'
#' Because the likelihood stays positive as \eqn{r \to \infty}{r -> Inf},
#' the posterior means of \eqn{r}{r} and \eqn{r \lambda_p}{r * lambda_p} are
#' infinite when the second shape parameter of the Beta prior is 1 or less,
#' as for the default \code{c(0.5, 0.5)}. They are then returned as
#' \code{Inf}, and the posterior median of \eqn{r}{r} is the summary to use.
#' Otherwise the means are integrated over the intervals of the rows, which
#' leave out a posterior mass of at most \code{control$eps} on each side;
#' when the second shape parameter is close to 1, the part of the means that
#' comes from the far right tail of \eqn{r}{r} is not negligible, and the
#' posterior median is again the summary to use.
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
#'   \code{n_coarse} (cells per coordinate of the coarse grid and per row,
#'   default 41), \code{n_lambda} (Gauss-Legendre nodes, default 48),
#'   \code{n_p} (Chebyshev points per row, default 49), \code{check}
#'   (default \code{TRUE}), \code{eps} (default \code{1e-8}) and
#'   \code{max_iter} (default 12).
#'
#' @return An object of class \code{bssd_post}; see \code{\link{bdribs_post}}.
#'   The element \code{model} is \code{"EXP_BERN"}, \code{estimates} contains
#'   the posterior means of \eqn{p}{p}, \eqn{r}{r}, \eqn{\lambda_p}{lambda_p}
#'   and \eqn{r \lambda_p}{r * lambda_p} (\code{lambda_d}), and
#'   \code{integration} contains the nodes used, the region of
#'   \eqn{(\log \lambda_p, z)}{(log(lambda_p), z)} located by the coarse
#'   grid and the error estimate.
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

  # subjects with the same (y, t) enter the likelihood once, with weights
  pat <- group_patterns(dat$y, dat$t)
  alpha <- lambda_hp * pt_hp / df
  beta <- pt_hp / df
  a <- p_prior[[1L]]
  b <- p_prior[[2L]]
  wp <- 1 / (k + 1)

  # log prior density of exp(u) times the Jacobian exp(u)
  lw_lambda <- function(u) {
    alpha * log(beta) - lgamma(alpha) + alpha * u - beta * exp(u)
  }
  # p and 1 - p at z = qnorm(F(p)), from the logarithm of the smaller tail
  # probability, so that both are accurate in the far tails of the prior
  pq_of_z <- function(z) {
    p <- z
    q <- z
    lt <- pnorm(-abs(z), log.p = TRUE)
    neg <- !is.na(z) & z <= 0
    pos <- !is.na(z) & z > 0
    p[neg] <- qbeta(lt[neg], a, b, log.p = TRUE)
    q[neg] <- 1 - p[neg]
    q[pos] <- qbeta(lt[pos], b, a, log.p = TRUE)
    p[pos] <- 1 - q[pos]
    list(p = p, q = q)
  }
  # p and 1 - p at v = F(p); 1 - p is computed from the Beta(b, a)
  # distribution for accuracy near p = 1
  pq_of_v <- function(v) {
    p <- v
    q <- v
    low <- v <= 0.5
    p[low] <- qbeta(v[low], a, b)
    q[low] <- 1 - p[low]
    q[!low] <- qbeta(1 - v[!low], b, a)
    p[!low] <- 1 - q[!low]
    list(p = p, q = q)
  }
  # log-likelihood of each row: placebo rate lam[i] with the rate ratios
  # P[i, ] / (k * Q[i, ])
  loglik <- function(lam, P, Q) {
    R <- P / (k * pmax(Q, 1e-300))
    out <- matrix(0, nrow(R), ncol(R))
    for (i in seq_len(nrow(R))) {
      out[i, ] <- expbern_loglik_cpp(pat$t, pat$y, pat$w, wp, lam[i], R[i, ])
    }
    out
  }
  logpost_z <- function(u, Z) {
    pq <- pq_of_z(Z)
    loglik(exp(u), pq$p, pq$q) + lw_lambda(u) + dnorm(Z, log = TRUE)
  }

  lo_u <- log(max(qgamma(1e-10, alpha, beta), .Machine$double.xmin))
  hi_u <- log(qgamma(1 - 1e-10, alpha, beta))
  box <- refine_box(function(nodes) {
    n1 <- length(nodes[[1L]])
    n2 <- length(nodes[[2L]])
    # logpost_z() of the matrix of z; z is the same in every row, so p and
    # 1 - p are computed once per column
    pq <- pq_of_z(nodes[[2L]])
    loglik(exp(nodes[[1L]]), matrix(pq$p, n1, n2, byrow = TRUE),
           matrix(pq$q, n1, n2, byrow = TRUE)) + lw_lambda(nodes[[1L]]) +
      dnorm(matrix(nodes[[2L]], n1, n2, byrow = TRUE), log = TRUE)
  }, c(lo_u, -8), c(hi_u, 8), ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
  lim_lo = c(-Inf, -Inf), lim_hi = c(Inf, Inf))

  # threshold p0 = k * rr0 / (k * rr0 + 1) of p on both scales, with
  # 1 - p0 = 1 / (k * rr0 + 1)
  lf0 <- pbeta(k * rr0 / (k * rr0 + 1), a, b, log.p = TRUE)
  ls0 <- pbeta(1 / (k * rr0 + 1), b, a, log.p = TRUE)
  v0 <- exp(lf0)
  z0 <- if (lf0 <= log(0.5)) {
    qnorm(lf0, log.p = TRUE)
  } else {
    -qnorm(ls0, log.p = TRUE)
  }
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  # grid with an interval of z located for each Gauss-Legendre node
  build <- function(n_l, n_p) {
    gq <- gl_panels(box$lo[1L], box$hi[1L], n_l)
    u <- gq$x
    rb <- row_refine(logpost_z, u, rep(box$lo[2L], n_l),
                     rep(box$hi[2L], n_l), ctrl$n_coarse, ctrl$eps,
                     ctrl$max_iter, rep(-Inf, n_l), rep(Inf, n_l))
    zl <- rb$lo
    zh <- rb$hi
    vl <- pnorm(zl)
    vh <- pnorm(zh)
    lam <- exp(u)
    # densities of each row at Chebyshev points in z and in v = pnorm(z)
    Zn <- row_nodes(zl, zh, n_p)
    # (rounding could move the end points out of [0, 1])
    Vn <- pmin(pmax(row_nodes(vl, vh, n_p), 0), 1)
    PZ <- pq_of_z(Zn)
    PV <- pq_of_v(Vn)
    Lz <- loglik(lam, PZ$p, PZ$q) + dnorm(Zn, log = TRUE)
    Lv <- loglik(lam, PV$p, PV$q)
    # v is used only where it is not 0 or 1 in double precision
    vrow <- zh <= 8 & zl >= -37 & row_tail_share(Lv) < row_tail_share(Lz)
    L <- Lz
    L[vrow, ] <- Lv[vrow, ]
    P <- PZ$p
    P[vrow, ] <- PV$p[vrow, ]
    Q <- PZ$q
    Q[vrow, ] <- PV$q[vrow, ]
    g <- post2d(L + lw_lambda(u), u, gq$w, ifelse(vrow, vl, zl),
                ifelse(vrow, vh, zh))
    list(g = g, lam = lam, P = P, Q = Q, vrow = vrow, zl = zl, zh = zh)
  }

  fine <- function(n_l, n_p, quantiles) {
    bg <- build(n_l, n_p)
    g <- bg$g
    prob <- 1 - g$cdf_rows(ifelse(bg$vrow, v0, z0))
    wq <- g$wq
    lam <- bg$lam
    R <- bg$P / (k * pmax(bg$Q, 1e-300))
    lam_d <- R * lam
    est <- c(p = sum(wq * bg$P), r = Inf, lambda_p = sum(wq * lam),
             lambda_d = Inf)
    if (b > 1) {
      est["r"] <- sum(wq * R)
      est["lambda_d"] <- sum(wq * lam_d)
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
      # P(z <= x) with the coordinate of each row
      cdf_z <- function(x) g$cdf_rows(ifelse(bg$vrow, pnorm(x), x))
      pq <- pq_of_z(inv_cdf(cdf_z, probs, min(bg$zl), max(bg$zh)))
      rq <- pq$p / (k * pmax(pq$q, 1e-300))
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
