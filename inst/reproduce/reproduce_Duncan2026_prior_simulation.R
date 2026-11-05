# Simulations for the sensitivity of the simulation results of Duncan (2026)
# to the prior of the models with a constant hazard rate (BDRIBS and
# EXP_BERN):
#   Duncan B (2026). Robust safety monitoring and signal detection using
#   alternatives to the standard Poisson distribution. Journal of
#   Biopharmaceutical Statistics 36(2):227-244.
#   doi:10.1080/10543406.2024.2395532
#
# This script only runs the simulations and saves the results; the
# judgements are made by reproduce_Duncan2026_prior.R.
#
# The article does not report the historical placebo data of the priors.
# reproduce_Duncan2026_simulation.R uses a historical placebo group of the
# same size as the current placebo group with the planned follow-up,
# discounted by 2 (df = 2), with the placebo incidence at the end of
# follow-up. With these priors, the statements D_Q_F5a and D_Q_F6a of
# reproduce_Duncan2026_figures.R fail. This script repeats the simulations
# of Figures 4 to 6 with other priors. All other settings are those of
# reproduce_Duncan2026_simulation.R.
# - Figure 4 (scenarios 3.1 to 3.4, final analysis at week 78): df = 1, that
#   is, a prior twice as strong. The final analysis needs no assumption on
#   the times of the interim analyses.
# - Figures 5 and 6 (scenarios 4.1 to 4.3, interim 3 at week 52): lambda_hp
#   is the constant hazard rate that gives the placebo incidence of 4% at
#   week 104 stated in Section 4.2.1, -log(0.96) / 104, instead of the
#   week 52 incidence of the log-logistic distribution (5.0%, which Table 1
#   gives as the week 104 incidence; see item D_T1_04P of
#   reproduce_Duncan2026.R). pt_hp = n / 2 * 52 and df = 2 as in the main
#   simulation.
# These priors were chosen after the main simulation by comparison with the
# article. Agreement with the article therefore explains the differences of
# the main simulation but is not an independent reproduction.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Duncan2026_prior_simulation.R")
# Output: inst/reproduce/results/Duncan2026_prior.rds
# Run time with n_sim_prior = 2000: about 30 minutes on one core (EXP_BERN
# takes most of it). n_sim_prior can be set before sourcing this script.

if (!exists("n_sim_prior")) n_sim_prior <- 2000L

accrual <- 26
rate_of <- function(ip, t) -log(1 - ip) / t
cutoff <- 0.95

tasks <- list()
for (i in 1:4) {
  ip_d <- c(0.30, 0.40, 0.45, 0.50)[i]
  tasks[[i]] <- list(
    figure = "4", point = i, main_id = sprintf("F4_S3_P%d", i),
    prior_label = "df = 1",
    design = list(n = 200, k = 1, accrual = accrual, follow_up = 52,
                  analysis_time = accrual + 52, dist = "exponential",
                  par_p = rate_of(0.30, 52), par_d = rate_of(ip_d, 52)),
    models = c("BDRIBS1", "BDRIBS2", "EXP_BERN"),
    priors = list(bdribs = list(lambda_hp = rate_of(0.30, 52),
                                pt_hp = 200 / 2 * 52, df = 1)))
}
for (i in 1:3) {
  mu_d <- c(9.85, 8.35, 7.4)[i]
  tasks[[4L + i]] <- list(
    figure = "5", point = i, main_id = sprintf("F5_S4_P%d", i),
    prior_label = "4% at week 104",
    design = list(n = 375, k = 1, accrual = accrual, follow_up = 52,
                  analysis_time = 3 * (52 / 3), dist = "loglogistic",
                  par_p = c(9.85, 2.0), par_d = c(mu_d, 2.0)),
    models = c("BDRIBS2", "EXP_BERN"),
    priors = list(bdribs = list(lambda_hp = rate_of(0.04, 104),
                                pt_hp = 375 / 2 * 52, df = 2)))
}

out <- vector("list", length(tasks))
for (j in seq_along(tasks)) {
  tk <- tasks[[j]]
  t0 <- proc.time()[["elapsed"]]
  out[[j]] <- sim_oc(n_sim_prior, tk$design, models = tk$models,
                     priors = tk$priors, rr0 = 1, cutoff = cutoff,
                     t_star = 52, seed = 20271105L + j, keep = TRUE)$results
  message(sprintf("[%d/%d] Figure %s, point %d: %.0f seconds", j,
                  length(tasks), tk$figure, tk$point,
                  proc.time()[["elapsed"]] - t0))
}

dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(n_sim = n_sim_prior, tasks = tasks, out = out),
        "inst/reproduce/results/Duncan2026_prior.rds")
