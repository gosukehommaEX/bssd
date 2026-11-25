#' Meta-analytic predictive prior of an event probability
#'
#' Derives the meta-analytic predictive (MAP) prior of the event probability
#' per unit time in a new control arm from historical control arms with the
#' hierarchical Beta model of Brock et al. (2023), and returns it as a
#' mixture of Beta distributions for \code{\link{mapbern_post}}.
#'
#' @details
#' Historical arm \eqn{j}{j} of \eqn{J}{J} has \eqn{N_j}{N_j} subjects, and
#' \eqn{N = \sum_j N_j}{N = sum(N_j)}. With the event probability
#' \eqn{\pi_j}{pi_j} per unit time and a constant hazard, a subject with time
#' at risk \eqn{t}{t} has the event with probability
#' \eqn{1 - (1 - \pi_j)^t}{1 - (1 - pi_j)^t}. The model is
#' \deqn{\pi_j \mid a, b \sim \mathrm{Beta}(a b N, (1 - a) b N), \quad a \sim \mathrm{Beta}(a_1, a_2), \quad b \sim \mathrm{Beta}(b_1, b_2),}{pi_j | a, b ~ Beta(a * b * N, (1 - a) * b * N), a ~ Beta(a1, a2), b ~ Beta(b1, b2),}
#' so \eqn{a}{a} is the mean event probability and \eqn{b}{b} discounts the
#' historical information. The MAP prior of the probability
#' \eqn{\pi^*}{pi*} in a new arm is
#' \deqn{p(\pi^*) = \int\int \mathrm{Beta}(\pi^*; a b N, (1 - a) b N) \, p(a, b \mid \mathrm{data}) \, da \, db.}{p(pi*) = integral of Beta(pi*; a * b * N, (1 - a) * b * N) * p(a, b | data) da db.}
#' The posterior of \eqn{(a, b)}{(a, b)} is integrated by Gauss-Legendre
#' quadrature on \eqn{(\mathrm{logit}(a), \mathrm{logit}(b))}{(logit(a), logit(b))}
#' after a coarse grid has located the region of the posterior mass; each
#' quadrature node gives one Beta component of the mixture, with the
#' quadrature weight times the posterior density as its weight. The
#' components with the smallest weights, which together have a weight of at
#' most \code{tol}, are dropped. Brock et al. (2023) sample the same model by
#' MCMC and approximate the predictive distribution by a mixture fitted to the
#' draws.
#'
#' The likelihood of an arm is the Beta-binomial probability when only the
#' numbers of events (\code{events}, \code{n}) are given, which assumes one
#' unit of time at risk for every subject. With subject-level data
#' (\code{y}, \code{t}, \code{study}), \eqn{\pi_j}{pi_j} is integrated out by
#' Gauss-Legendre quadrature on \eqn{\mathrm{logit}(\pi_j)}{logit(pi_j)}.
#'
#' @param events Numbers of subjects with the event in the historical arms
#'   (summary data with one unit of time at risk per subject).
#' @param n Numbers of subjects in the historical arms (summary data).
#' @param y,t,study Subject-level data instead of \code{events} and
#'   \code{n}: event indicators, times at risk, and arm identifiers.
#' @param a_prior Shape parameters of the Beta prior of \eqn{a}{a}; the
#'   default is that of Brock et al. (2023).
#' @param b_prior Shape parameters of the Beta prior of \eqn{b}{b}; the
#'   default is that of Brock et al. (2023).
#' @param tol Total weight of the components that may be dropped.
#' @param control List of settings of the numerical integration:
#'   \code{n_coarse} (default 41), \code{n_ab} (Gauss-Legendre nodes per
#'   coordinate, default 32), \code{n_pi} (nodes per panel for
#'   subject-level data, default 16), \code{eps} (default 1e-8),
#'   \code{max_iter} (default 12) and \code{check} (default \code{TRUE}).
#'
#' @return A data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}, one row per component of the mixture, ordered by
#'   decreasing weight. Its attributes are \code{n_total} (\eqn{N}{N}),
#'   \code{hyper} (posterior means of \eqn{a}{a} and \eqn{b}{b}),
#'   \code{moments} (mean and standard deviation of the MAP prior) and
#'   \code{error} (with \code{control$check = TRUE}, the largest absolute
#'   difference of the MAP prior distribution function at its 2.5, 50 and
#'   97.5 percent points between the default and a finer grid).
#'
#' @references
#' Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
#' Penard, N. (2023). A Bayesian method for safety signal detection in
#' ongoing blinded randomised controlled trials. \emph{Pharmaceutical
#' Statistics}, 22(2), 378--395. \doi{10.1002/pst.2278}
#'
#' @seealso \code{\link{mapbern_post}}
#'
#' @examples
#' prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
#' attr(prior, "moments")
#' nrow(prior)
#'
#' @export
map_beta <- function(events = NULL, n = NULL, y = NULL, t = NULL,
                     study = NULL, a_prior = c(1 / 3, 1 / 3),
                     b_prior = c(5, 5), tol = 1e-10, control = list()) {
  ctrl <- merge_control(control, list(n_coarse = 41L, n_ab = 32L,
                                      n_pi = 16L, eps = 1e-8,
                                      max_iter = 12L, check = TRUE))
  check_num(a_prior, "a_prior", len = 2L, lower = 0, lower_open = TRUE)
  check_num(b_prior, "b_prior", len = 2L, lower = 0, lower_open = TRUE)
  check_num(tol, "tol", lower = 0, upper = 0.01)
  summary_data <- !is.null(events) || !is.null(n)
  if (summary_data) {
    if (!is.null(y) || !is.null(t) || !is.null(study)) {
      stop("Give either 'events' and 'n' or 'y', 't' and 'study'.",
           call. = FALSE)
    }
    check_num(n, "n", len = NULL, lower = 1, integer = TRUE)
    check_num(events, "events", len = length(n), lower = 0, integer = TRUE)
    if (any(events > n)) {
      stop("'events' must not exceed 'n'.", call. = FALSE)
    }
    n_j <- round(n)
    arm_ll <- NULL
  } else {
    dat <- check_subject_data(y, t)
    if (is.null(study) || length(study) != length(dat$y)) {
      stop("'study' must have the length of 'y'.", call. = FALSE)
    }
    arms <- split(seq_along(dat$y), as.character(study))
    n_j <- vapply(arms, length, numeric(1))
    arm_ll <- lapply(arms, function(idx) {
      map_arm_grid(dat$y[idx], dat$t[idx], ctrl$n_pi)
    })
  }
  n_total <- sum(n_j)

  # log marginal likelihood of all arms for vectors of shape parameters
  log_marg <- function(s1, s2) {
    out <- 0
    if (summary_data) {
      for (j in seq_along(n_j)) {
        out <- out + lbeta(s1 + events[j], s2 + n_j[j] - events[j]) -
          lbeta(s1, s2)
      }
    } else {
      for (g in arm_ll) {
        L <- outer(s1, g$l1) + outer(s2, g$l0) +
          rep(g$lw, each = length(s1))
        m <- apply(L, 1L, max)
        out <- out + m + log(rowSums(exp(L - m))) - lbeta(s1, s2)
      }
    }
    out
  }
  lp_logit <- function(x, p) {
    p[1L] * plogis(x, log.p = TRUE) + p[2L] * plogis(-x, log.p = TRUE)
  }
  logpost <- function(nodes) {
    a <- plogis(nodes[[1L]])
    b <- plogis(nodes[[2L]])
    na <- length(a)
    s1 <- as.vector(outer(a, b)) * n_total
    s2 <- as.vector(outer(1 - a, b)) * n_total
    L <- matrix(log_marg(s1, s2), na, length(b))
    L + outer(lp_logit(nodes[[1L]], a_prior), lp_logit(nodes[[2L]], b_prior),
              "+")
  }
  box <- refine_box(logpost, c(-15, -15), c(15, 15), ctrl$n_coarse,
                    ctrl$eps, ctrl$max_iter, lim_lo = c(-40, -40),
                    lim_hi = c(40, 40))

  build <- function(n_ab) {
    g <- gl_nodes(n_ab)
    ua <- (box$hi[1L] - box$lo[1L]) / 2 * g$x + (box$hi[1L] + box$lo[1L]) / 2
    ub <- (box$hi[2L] - box$lo[2L]) / 2 * g$x + (box$hi[2L] + box$lo[2L]) / 2
    wa <- g$w * (box$hi[1L] - box$lo[1L]) / 2
    wb <- g$w * (box$hi[2L] - box$lo[2L]) / 2
    L <- logpost(list(ua, ub)) + outer(log(wa), log(wb), "+")
    L[is.na(L)] <- -Inf
    w <- exp(L - max(L))
    w <- as.vector(w / sum(w))
    a <- rep(plogis(ua), times = n_ab)
    b <- rep(plogis(ub), each = n_ab)
    o <- order(w, decreasing = TRUE)
    keep <- o[seq_len(which(cumsum(w[o]) >= 1 - tol)[1L])]
    mix <- data.frame(w = w[keep] / sum(w[keep]),
                      shape1 = a[keep] * b[keep] * n_total,
                      shape2 = (1 - a[keep]) * b[keep] * n_total)
    attr(mix, "hyper") <- c(a = sum(w * a), b = sum(w * b))
    mix
  }
  mix <- build(ctrl$n_ab)
  err <- NA_real_
  if (isTRUE(ctrl$check)) {
    coarse <- mix
    mix <- build(as.integer(round(1.5 * ctrl$n_ab)))
    pts <- qmixbeta(c(0.025, 0.5, 0.975), mix)
    err <- max(abs(pmixbeta(pts, coarse) - pmixbeta(pts, mix)))
  }
  m1 <- sum(mix$w * mix$shape1 / (mix$shape1 + mix$shape2))
  m2 <- sum(mix$w * mix$shape1 * (mix$shape1 + 1) /
              ((mix$shape1 + mix$shape2) * (mix$shape1 + mix$shape2 + 1)))
  rownames(mix) <- NULL
  attr(mix, "n_total") <- n_total
  attr(mix, "moments") <- c(mean = m1, sd = sqrt(max(m2 - m1^2, 0)))
  attr(mix, "error") <- err
  mix
}
