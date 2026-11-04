#' Operating characteristics of the signal detection models by simulation
#'
#' Simulates blinded trials with \code{\link{sim_blinded}}, applies the
#' selected models at each analysis, and summarizes the probability of
#' declaring a signal (power, or the false positive rate when the study drug
#' and placebo have the same distribution) and the posterior mean of the
#' difference in cumulative incidence.
#'
#' @details
#' The models are \code{"BDRIBS1"} (\code{\link{bdribs_post}} with the total
#' exposure time), \code{"BDRIBS2"} (\code{\link{bdribs_post}} with the total
#' time at risk), \code{"EXP_BERN"} (\code{\link{expbern_post}}) and
#' \code{"LL_BERN"} (\code{\link{llbern_post}}). The posterior functions are
#' called with \code{details = FALSE}. A fit that fails is recorded as
#' \code{NA} and counted in \code{n_failed}, and a warning with the first
#' error message is given.
#'
#' @param n_sim Number of simulated trials.
#' @param design List of arguments of \code{\link{sim_blinded}}.
#' @param models Character vector of models to apply.
#' @param priors List of prior settings: \code{bdribs} (arguments
#'   \code{lambda_hp}, \code{pt_hp} and optionally \code{df} and
#'   \code{p_prior} of \code{\link{bdribs_post}}, used by \code{"BDRIBS1"} and
#'   \code{"BDRIBS2"}), \code{expbern} (the same arguments for
#'   \code{\link{expbern_post}}; \code{bdribs} is used when it is not given)
#'   and \code{llbern} (arguments \code{mu_hp}, \code{sigma_hp}, \code{pt_hp}
#'   and optionally \code{df_p}, \code{mu_d_prior} and \code{sigma_d_prior} of
#'   \code{\link{llbern_post}}).
#' @param rr0 Critical value of the relative risk.
#' @param cutoff Posterior probability threshold for a signal.
#' @param t_star Time at which the relative risk of \code{"LL_BERN"} and the
#'   cumulative incidences are evaluated. Required for \code{"LL_BERN"}.
#' @param control List with elements \code{EXP_BERN} and \code{LL_BERN}, the
#'   \code{control} arguments of \code{\link{expbern_post}} and
#'   \code{\link{llbern_post}}. By default the convergence check is turned
#'   off.
#' @param seed Optional seed passed to \code{\link{set.seed}}.
#' @param keep Logical; if \code{TRUE}, the results of each simulated trial
#'   are returned.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{summary}}{Data frame with one row per analysis and model and
#'     the columns \code{look}, \code{analysis_time}, \code{model},
#'     \code{n_ok} (successful fits), \code{n_failed}, \code{prob_signal}
#'     (proportion of trials with a signal), \code{se} (its Monte Carlo
#'     standard error), \code{mean_prob} (mean posterior probability),
#'     \code{mean_ip_difference} (mean of the posterior means of the
#'     difference in cumulative incidence at \code{t_star}), \code{mean_n}
#'     and \code{mean_events}.}
#'   \item{\code{truth}}{True cumulative incidences at \code{t_star} for
#'     placebo and study drug and their difference.}
#'   \item{\code{results}}{Data frame of the results of each trial (only when
#'     \code{keep = TRUE}).}
#' }
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{sim_blinded}}
#'
#' @examples
#' design <- list(n = 100, k = 1, accrual = 0, follow_up = 52,
#'                dist = "exponential", par_p = -log(0.9) / 52,
#'                par_d = -log(0.75) / 52)
#' oc <- sim_oc(n_sim = 20, design = design, models = c("BDRIBS1", "BDRIBS2"),
#'              priors = list(bdribs = list(lambda_hp = -log(0.9) / 52,
#'                                          pt_hp = 2000)),
#'              t_star = 52, seed = 1)
#' oc$summary
#'
#' @export
sim_oc <- function(n_sim, design,
                   models = c("BDRIBS1", "BDRIBS2", "EXP_BERN", "LL_BERN"),
                   priors, rr0 = 1, cutoff = 0.95, t_star = NULL,
                   control = list(), seed = NULL, keep = FALSE) {
  check_num(n_sim, "n_sim", lower = 1, integer = TRUE)
  if (!is.list(design)) {
    stop("'design' must be a list of arguments of sim_blinded().",
         call. = FALSE)
  }
  all_models <- c("BDRIBS1", "BDRIBS2", "EXP_BERN", "LL_BERN")
  if (!is.character(models) || length(models) == 0L ||
      any(!models %in% all_models)) {
    stop("'models' must be a subset of ",
         paste(all_models, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(priors)) {
    stop("'priors' must be a list.", call. = FALSE)
  }
  if (is.null(priors$expbern)) {
    priors$expbern <- priors$bdribs
  }
  need <- c(BDRIBS1 = "bdribs", BDRIBS2 = "bdribs", EXP_BERN = "expbern",
            LL_BERN = "llbern")[models]
  missing_prior <- unique(need[vapply(need, function(p) is.null(priors[[p]]),
                                      logical(1))])
  if (length(missing_prior) > 0L) {
    stop("'priors' must contain: ", paste(missing_prior, collapse = ", "),
         ".", call. = FALSE)
  }
  if ("LL_BERN" %in% models && is.null(t_star)) {
    stop("'t_star' must be given for LL_BERN.", call. = FALSE)
  }
  ctrl <- list(EXP_BERN = list(check = FALSE), LL_BERN = list(check = FALSE))
  ctrl[names(control)] <- control
  k <- if (is.null(design$k)) 1 else design$k
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_look <- length(if (is.null(design$analysis_time)) 1 else
    design$analysis_time)
  n_row <- n_sim * n_look * length(models)
  res <- data.frame(sim = integer(n_row), look = integer(n_row),
                    analysis_time = numeric(n_row),
                    model = character(n_row), n = integer(n_row),
                    n_event = integer(n_row), prob = numeric(n_row),
                    ip_difference = numeric(n_row),
                    stringsAsFactors = FALSE)
  i <- 0L
  first_error <- NULL
  for (s in seq_len(n_sim)) {
    dat <- do.call(sim_blinded, design)
    for (j in seq_len(n_look)) {
      d <- dat[dat$look == j, , drop = FALSE]
      for (m in models) {
        i <- i + 1L
        fit <- tryCatch(fit_oc_model(m, d, k, priors, rr0, cutoff, t_star,
                                     ctrl),
                        error = function(e) {
                          if (is.null(first_error)) {
                            first_error <<- conditionMessage(e)
                          }
                          NULL
                        })
        res$sim[i] <- s
        res$look[i] <- j
        res$analysis_time[i] <- if (nrow(d) > 0L) d$analysis_time[1L] else
          NA_real_
        res$model[i] <- m
        res$n[i] <- nrow(d)
        res$n_event[i] <- sum(d$event)
        res$prob[i] <- if (is.null(fit)) NA_real_ else fit$prob
        res$ip_difference[i] <- if (is.null(fit) || is.null(fit$incidence)) {
          NA_real_
        } else {
          unname(fit$incidence["difference"])
        }
      }
    }
  }
  res$signal <- res$prob > cutoff
  if (!is.null(first_error)) {
    warning(sum(is.na(res$prob)), " of ", n_row, " fits failed; the first ",
            "error was: ", first_error, call. = FALSE)
  }

  keys <- unique(res[, c("look", "model")])
  keys <- keys[order(keys$look, match(keys$model, all_models)), ]
  summ <- lapply(seq_len(nrow(keys)), function(r) {
    x <- res[res$look == keys$look[r] & res$model == keys$model[r], ]
    ok <- !is.na(x$prob)
    p <- mean(x$signal[ok])
    data.frame(look = keys$look[r],
               analysis_time = mean(x$analysis_time, na.rm = TRUE),
               model = keys$model[r], n_ok = sum(ok), n_failed = sum(!ok),
               prob_signal = p, se = sqrt(p * (1 - p) / sum(ok)),
               mean_prob = mean(x$prob[ok]),
               mean_ip_difference = mean(x$ip_difference[ok]),
               mean_n = mean(x$n), mean_events = mean(x$n_event),
               stringsAsFactors = FALSE)
  })
  summ <- do.call(rbind, summ)
  rownames(summ) <- NULL

  truth <- NULL
  if (!is.null(t_star)) {
    dist <- if (is.null(design$dist)) "exponential" else
      match.arg(design$dist, c("exponential", "loglogistic"))
    if (dist == "exponential") {
      fp <- -expm1(-design$par_p * t_star)
      fd <- -expm1(-design$par_d * t_star)
    } else {
      fp <- ll_cdf(t_star, design$par_p[1L], design$par_p[2L])
      fd <- ll_cdf(t_star, design$par_d[1L], design$par_d[2L])
    }
    truth <- c(placebo = fp, drug = fd, difference = fd - fp)
  }

  out <- list(summary = summ, truth = truth)
  if (keep) {
    out$results <- res
  }
  out
}
