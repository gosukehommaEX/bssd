#' Prior parameters from two quantiles
#'
#' Finds the parameters of a Gamma, Beta or lognormal distribution with two
#' given quantiles, or of a lognormal distribution with a given mean and a
#' given probability content of an interval.
#'
#' @details
#' With \code{q = c(q1, q2)} and \code{p = c(p1, p2)}, the parameters solve
#' \eqn{F(q_1) = p_1}{F(q1) = p1} and \eqn{F(q_2) = p_2}{F(q2) = p2}. For the
#' Gamma distribution, the rate is \eqn{G^{-1}(p_1; s) / q_1}{qgamma(p1, s) / q1}
#' for shape \eqn{s}{s}, and the shape solves the second equation; for the
#' Beta distribution, the second shape parameter solves the first equation
#' for each value of the first one, which then solves the second equation;
#' for the lognormal distribution the solution is explicit. The equations
#' are solved by \code{\link[stats]{uniroot}}, so the quantiles are matched
#' to about 10 significant digits. Schnell and Ball (2016) and Gould and
#' Wang (2017, Table I) elicit their priors in this way; their own software
#' minimizes the squared differences of the probabilities instead.
#'
#' For \code{family = "lnorm"} with \code{mean} given, \code{q} is an
#' interval \code{c(lower, upper)} and \code{p} is its probability content:
#' the lognormal distribution has mean \code{mean} and probability \code{p}
#' between \code{lower} and \code{upper}. This is the prior of the risk ratio
#' of Gould and Wang (2017).
#'
#' @param q Two increasing quantiles, or the interval when \code{mean} is
#'   given.
#' @param p Two increasing probabilities, or the probability content of the
#'   interval when \code{mean} is given.
#' @param family One of \code{"gamma"}, \code{"beta"} and \code{"lnorm"}.
#' @param mean Optional mean of the lognormal distribution.
#'
#' @return A named numeric vector: \code{c(shape, rate)} for the Gamma
#'   distribution, \code{c(shape1, shape2)} for the Beta distribution and
#'   \code{c(meanlog, sdlog)} for the lognormal distribution.
#'
#' @references
#' Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
#' rate differences using data from blinded trials: the canary in the coal
#' mine. \emph{Statistics in Medicine}, 36(1), 92--104.
#' \doi{10.1002/sim.7129}
#'
#' Schnell, P. M. and Ball, G. (2016). A Bayesian exposure-time method for
#' clinical trial safety monitoring with blinded data. \emph{Therapeutic
#' Innovation & Regulatory Science}, 50(6), 833--838.
#' \doi{10.1177/2168479016656702}
#'
#' @examples
#' # 50 percent sure that the rate is below 0.5 and 75 percent sure that it
#' # is below 2 (Schnell and Ball, 2016)
#' elicit_prior(c(0.5, 2), c(0.5, 0.75), "gamma")
#' # median 0.02 and 95th percentile 0.04 of a proportion
#' elicit_prior(c(0.02, 0.04), c(0.5, 0.95), "beta")
#' # lognormal prior of the risk ratio of Gould and Wang (2017)
#' elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)
#'
#' @export
elicit_prior <- function(q, p, family = c("gamma", "beta", "lnorm"),
                         mean = NULL) {
  family <- match.arg(family)
  if (!is.null(mean)) {
    if (family != "lnorm") {
      stop("'mean' can only be given for family = \"lnorm\".", call. = FALSE)
    }
    check_num(mean, "mean", lower = 0, lower_open = TRUE)
    check_num(q, "q", len = 2L, lower = 0, lower_open = TRUE)
    check_num(p, "p", lower = 0, upper = 1, lower_open = TRUE,
              upper_open = TRUE)
    if (q[1L] >= q[2L] || mean <= q[1L] || mean >= q[2L]) {
      stop("'q' must be an increasing interval that contains 'mean'.",
           call. = FALSE)
    }
    content <- function(s) {
      mu <- log(mean) - s^2 / 2
      plnorm(q[2L], mu, s) - plnorm(q[1L], mu, s)
    }
    # the content tends to 1 as s -> 0 and to 0 as s -> Inf
    s <- uniroot(function(ls) content(exp(ls)) - p, c(-10, 5),
                 tol = 1e-12, extendInt = "downX")$root
    s <- exp(s)
    return(c(meanlog = log(mean) - s^2 / 2, sdlog = s))
  }
  check_num(q, "q", len = 2L)
  check_num(p, "p", len = 2L, lower = 0, upper = 1, lower_open = TRUE,
            upper_open = TRUE)
  if (q[1L] >= q[2L] || p[1L] >= p[2L]) {
    stop("'q' and 'p' must be increasing.", call. = FALSE)
  }
  if (family == "lnorm") {
    check_num(q, "q", len = 2L, lower = 0, lower_open = TRUE)
    z <- qnorm(p)
    s <- (log(q[2L]) - log(q[1L])) / (z[2L] - z[1L])
    return(c(meanlog = log(q[1L]) - s * z[1L], sdlog = s))
  }
  if (family == "gamma") {
    check_num(q, "q", len = 2L, lower = 0, lower_open = TRUE)
    rate_of <- function(s) qgamma(p[1L], s) / q[1L]
    f_gamma <- function(ls) {
      s <- exp(ls)
      pgamma(q[2L] * rate_of(s), s) - p[2L]
    }
    ls <- uniroot(f_gamma, c(-5, 5), tol = 1e-12, extendInt = "upX")$root
    return(c(shape = exp(ls), rate = rate_of(exp(ls))))
  }
  check_num(q, "q", len = 2L, lower = 0, upper = 1, lower_open = TRUE,
            upper_open = TRUE)
  # shape2 that gives F(q1) = p1 for shape1 = s; F(q1) increases in shape2
  b_of <- function(s) {
    exp(uniroot(function(lb) pbeta(q[1L], s, exp(lb)) - p[1L], c(-5, 5),
                tol = 1e-12, extendInt = "upX")$root)
  }
  f_beta <- function(la) {
    s <- exp(la)
    pbeta(q[2L], s, b_of(s)) - p[2L]
  }
  grid <- seq(-6, 8, by = 0.5)
  vals <- vapply(grid, function(la) {
    tryCatch(f_beta(la), error = function(e) NA_real_)
  }, numeric(1))
  i <- which(!is.na(vals[-1L]) & !is.na(vals[-length(vals)]) &
               sign(vals[-1L]) != sign(vals[-length(vals)]))
  if (length(i) == 0L) {
    stop("No Beta distribution has these quantiles.", call. = FALSE)
  }
  la <- uniroot(f_beta, grid[c(i[1L], i[1L] + 1L)], tol = 1e-12)$root
  c(shape1 = exp(la), shape2 = b_of(exp(la)))
}
