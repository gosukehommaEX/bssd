#' Fit one model to one simulated blinded data set
#'
#' @param model One of \code{"BDRIBS1"}, \code{"BDRIBS2"}, \code{"EXP_BERN"}
#'   and \code{"LL_BERN"}.
#' @param d Data frame of one analysis from \code{sim_blinded}.
#' @param k Randomization ratio.
#' @param priors List of prior settings (see \code{sim_oc}).
#' @param rr0,cutoff,t_star Settings passed to the posterior function.
#' @param control List of \code{control} arguments for \code{"EXP_BERN"} and
#'   \code{"LL_BERN"}.
#'
#' @return An object of class \code{bssd_post}.
#' @keywords internal
#' @noRd
fit_oc_model <- function(model, d, k, priors, rr0, cutoff, t_star, control) {
  common <- list(k = k, rr0 = rr0, cutoff = cutoff, t_star = t_star,
                 details = FALSE)
  ctrl_of <- function(m) if (is.null(control[[m]])) list() else control[[m]]
  switch(model,
         BDRIBS1 = do.call(bdribs_post,
                           c(list(y = sum(d$event),
                                  exposure = sum(d$follow_up)),
                             common, priors$bdribs)),
         BDRIBS2 = do.call(bdribs_post,
                           c(list(y = sum(d$event),
                                  exposure = sum(d$time_at_risk)),
                             common, priors$bdribs)),
         EXP_BERN = do.call(expbern_post,
                            c(list(y = d$event, t = d$follow_up), common,
                              list(control = ctrl_of("EXP_BERN")),
                              priors$expbern)),
         LL_BERN = do.call(llbern_post,
                           c(list(y = d$event, t = d$follow_up), common,
                             list(control = ctrl_of("LL_BERN")),
                             priors$llbern)))
}
