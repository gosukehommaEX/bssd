# Simulations for the approximate reproduction of Figures 2 to 10 of
#   Duncan B (2026). Robust safety monitoring and signal detection using
#   alternatives to the standard Poisson distribution. Journal of
#   Biopharmaceutical Statistics 36(2):227-244.
#   doi:10.1080/10543406.2024.2395532
#
# This script only runs the simulations and saves the results; the figures
# and the judgements are made by reproduce_Duncan2026_figures.R.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Duncan2026_simulation.R")
# With n_cores > 1 the scenarios run in parallel, which needs the installed
# package (run devtools::install() first).
#
# Output (inst/reproduce/results/):
#   Duncan2026_simulation.rds         results of every simulated trial
#   Duncan2026_simulation_design.csv  settings of every scenario
# Run time with n_sim = 1000: about 2 to 3 hours on one core.
#
# The article does not report the accrual, the times of the interim analyses,
# the historical patient time, the discount factors or the number of
# simulated trials. The following assumptions are used (time in weeks):
# - k = 1 and uniform accrual over 26 weeks.
# - Analyses are equally spaced in calendar time: every 52 / 3 weeks (4
#   months) for 52-week follow-up and every 39 weeks (9 months) for 104-week
#   follow-up. With these spacings the risk differences in the simulated data
#   at interim 4 of Figure 3, interim 3 of Figure 5 and interim 3 of Figure 9
#   agree with the x-axes of these figures. The final analysis is at the end
#   of follow-up of the last subject (week 78 or 130). Interim 4 of scenario
#   8.1 would fall after the final analysis and is replaced by the final
#   analysis.
# - Scenario 4 uses 52-week follow-up, because the incidences in Table 1 are
#   the week 52 values of the stated log-logistic parameters (see
#   reproduce_Duncan2026.R).
# - Priors: a historical placebo group of the same size as the current
#   placebo group with the planned follow-up, discounted by 2 (df = 2). For
#   BDRIBS and EXP_BERN, lambda_hp is the constant hazard rate that gives the
#   assumed placebo incidence at the end of follow-up and pt_hp is in
#   patient-weeks; for LL_BERN, mu_hp and sigma_hp are the assumed placebo
#   parameters and pt_hp is in patient-years.
# - Scenario 5.2 (constant hazard) uses the LL_BERN prior of scenario 5.1,
#   that is, LL_BERN is misspecified.
# - The relative risk of LL_BERN and the incidences are evaluated at the end
#   of follow-up (t_star = 52 or 104).
# - The posterior probabilities of every trial are kept, so that the signal
#   rates can be computed for any threshold (0.925 and 0.95 are used).

# n_sim and n_cores can be set before sourcing this script, for example
#   n_sim <- 20; n_cores <- 1   (a quick run of a few minutes to check)
#   n_sim <- 1000; n_cores <- 8 (needs devtools::install() first)
if (!exists("n_sim")) n_sim <- 1000L
if (!exists("n_cores")) n_cores <- 1L

accrual <- 26
step_52 <- 52 / 3
step_104 <- 39
rate_of <- function(ip, t) -log(1 - ip) / t

prior_const <- function(n, follow_up, ip_p) {
  list(lambda_hp = rate_of(ip_p, follow_up), pt_hp = n / 2 * follow_up,
       df = 2)
}
prior_ll <- function(n, follow_up, par) {
  list(mu_hp = par[1], sigma_hp = par[2], pt_hp = n / 2 * follow_up / 52,
       df_p = 2)
}

tasks <- list()
add_task <- function(figure, scenario, point, label, n, follow_up, look,
                     analysis_time, dist, par_p, par_d, models, priors) {
  tasks[[length(tasks) + 1L]] <<- list(
    id = sprintf("F%s_S%s_P%d", figure, scenario, point), figure = figure,
    scenario = scenario, point = point, label = label, look = look,
    design = list(n = n, k = 1, accrual = accrual, follow_up = follow_up,
                  analysis_time = analysis_time, dist = dist, par_p = par_p,
                  par_d = par_d),
    models = models, priors = priors, t_star = follow_up,
    seed = 20261104L + length(tasks))
}
m_const <- c("BDRIBS1", "BDRIBS2", "EXP_BERN")
m_ll <- c("BDRIBS2", "EXP_BERN", "LL_BERN")

# Figure 2 (scenarios 1.1 to 1.4): 0.5% and 4.0% at week 52, interim 2
for (i in 1:4) {
  n <- c(150, 225, 300, 375)[i]
  add_task("2", "1", i, sprintf("N = %d", n), n, 52, 2, 2 * step_52,
           "exponential", rate_of(0.005, 52), rate_of(0.040, 52), m_const,
           list(bdribs = prior_const(n, 52, 0.005)))
}
# Figure 3 (scenarios 2.1 to 2.4): 10% and 10% to 25% at week 52, interim 4
for (i in 1:4) {
  ip_d <- c(0.10, 0.15, 0.20, 0.25)[i]
  add_task("3", "2", i, sprintf("drug %.0f%%", 100 * ip_d), 350, 52, 4,
           4 * step_52, "exponential", rate_of(0.10, 52), rate_of(ip_d, 52),
           m_const, list(bdribs = prior_const(350, 52, 0.10)))
}
# Figure 4 (scenarios 3.1 to 3.4): 30% and 30% to 50% at week 52, final
for (i in 1:4) {
  ip_d <- c(0.30, 0.40, 0.45, 0.50)[i]
  add_task("4", "3", i, sprintf("drug %.0f%%", 100 * ip_d), 200, 52, 0,
           accrual + 52, "exponential", rate_of(0.30, 52), rate_of(ip_d, 52),
           m_const, list(bdribs = prior_const(200, 52, 0.30)))
}
# Figures 5 and 6 (scenarios 4.1 to 4.3): decreasing hazard, interim 3
for (i in 1:3) {
  mu_d <- c(9.85, 8.35, 7.4)[i]
  add_task("5", "4", i, sprintf("mu_d = %g", mu_d), 375, 52, 3, 3 * step_52,
           "loglogistic", c(9.85, 2.0), c(mu_d, 2.0), m_ll,
           list(bdribs = prior_const(375, 52, ll_cdf(52, 9.85, 2.0)),
                llbern = prior_ll(375, 52, c(9.85, 2.0))))
}
# Figure 7 (scenario 5.1): decreasing hazard, equal groups, interims 1 to 3
add_task("7", "5.1", 1, "null", 400, 52, 1:3, (1:3) * step_52,
         "loglogistic", c(6.17, 1.6), c(6.17, 1.6), m_ll,
         list(bdribs = prior_const(400, 52, ll_cdf(52, 6.17, 1.6)),
              llbern = prior_ll(400, 52, c(6.17, 1.6))))
# Figure 8 (scenario 5.2): constant hazard with the same week 52 incidence
add_task("8", "5.2", 1, "null", 400, 52, 1:3, (1:3) * step_52,
         "exponential", rate_of(ll_cdf(52, 6.17, 1.6), 52),
         rate_of(ll_cdf(52, 6.17, 1.6), 52), m_ll,
         list(bdribs = prior_const(400, 52, ll_cdf(52, 6.17, 1.6)),
              llbern = prior_ll(400, 52, c(6.17, 1.6))))
# Figure 9 (scenarios 6.1 to 6.4): increasing hazard, interim 3
for (i in 1:4) {
  mu_d <- c(7.3, 6.5, 6.1, 5.85)[i]
  add_task("9", "6", i, sprintf("mu_d = %g", mu_d), 150, 104, 3,
           3 * step_104, "loglogistic", c(7.3, 0.8), c(mu_d, 0.8), m_ll,
           list(bdribs = prior_const(150, 104, ll_cdf(104, 7.3, 0.8)),
                llbern = prior_ll(150, 104, c(7.3, 0.8))))
}
# Figure 10 (scenarios 7.1 to 7.4): increasing hazard, final analysis
par_d7 <- list(c(5.85, 0.7), c(5.5, 0.7), c(5.15, 0.6), c(5.07, 0.7))
for (i in 1:4) {
  add_task("10", "7", i, sprintf("(%g, %g)", par_d7[[i]][1], par_d7[[i]][2]),
           150, 104, 0, accrual + 104, "loglogistic", c(5.85, 0.7),
           par_d7[[i]], m_ll,
           list(bdribs = prior_const(150, 104, ll_cdf(104, 5.85, 0.7)),
                llbern = prior_ll(150, 104, c(5.85, 0.7))))
}
# Scenario 8.1: placebo rate underestimated by the prior, equal groups
add_task("text_8.1", "8.1", 1, "null", 150, 104, 0, accrual + 104,
         "loglogistic", c(6.5, 0.8), c(6.5, 0.8), m_ll,
         list(bdribs = prior_const(150, 104, ll_cdf(104, 6.85, 0.8)),
              llbern = prior_ll(150, 104, c(6.85, 0.8))))
# Scenario 8.2: placebo rate overestimated by the prior, interim 3
add_task("text_8.2", "8.2", 1, "drug 18.1%", 200, 104, 3, 3 * step_104,
         "loglogistic", c(6.845, 0.8), c(5.85, 0.8), m_ll,
         list(bdribs = prior_const(200, 104, ll_cdf(104, 6.6, 0.8)),
              llbern = prior_ll(200, 104, c(6.6, 0.8))))

run_task <- function(task, n_sim) {
  t0 <- proc.time()[["elapsed"]]
  oc <- sim_oc(n_sim, task$design, models = task$models,
               priors = task$priors, rr0 = 1, cutoff = 0.95,
               t_star = task$t_star, seed = task$seed, keep = TRUE)
  list(id = task$id, summary = oc$summary, truth = oc$truth,
       results = oc$results, seconds = proc.time()[["elapsed"]] - t0)
}

if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, library(bssd)))
  out <- parallel::parLapplyLB(cl, tasks, run_task, n_sim = n_sim)
  parallel::stopCluster(cl)
} else {
  out <- vector("list", length(tasks))
  for (i in seq_along(tasks)) {
    out[[i]] <- run_task(tasks[[i]], n_sim)
    message(sprintf("[%d/%d] %s: %.0f seconds", i, length(tasks),
                    tasks[[i]]$id, out[[i]]$seconds))
  }
}

dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(n_sim = n_sim, tasks = tasks, out = out),
        "inst/reproduce/results/Duncan2026_simulation.rds")
design <- do.call(rbind, lapply(tasks, function(tk) {
  data.frame(id = tk$id, figure = tk$figure, scenario = tk$scenario,
             point = tk$point, label = tk$label, n = tk$design$n,
             follow_up = tk$design$follow_up, accrual = tk$design$accrual,
             analysis_time = paste(round(tk$design$analysis_time, 2),
                                   collapse = " "),
             dist = tk$design$dist,
             par_p = paste(signif(tk$design$par_p, 6), collapse = " "),
             par_d = paste(signif(tk$design$par_d, 6), collapse = " "),
             models = paste(tk$models, collapse = " "),
             lambda_hp = signif(tk$priors$bdribs$lambda_hp, 6),
             pt_hp_weeks = tk$priors$bdribs$pt_hp,
             df = tk$priors$bdribs$df,
             mu_hp = if (is.null(tk$priors$llbern)) NA else
               tk$priors$llbern$mu_hp,
             sigma_hp = if (is.null(tk$priors$llbern)) NA else
               tk$priors$llbern$sigma_hp,
             pt_hp_years = if (is.null(tk$priors$llbern)) NA else
               tk$priors$llbern$pt_hp,
             t_star = tk$t_star, seed = tk$seed, n_sim = n_sim,
             stringsAsFactors = FALSE)
}))
utils::write.csv(design,
                 "inst/reproduce/results/Duncan2026_simulation_design.csv",
                 row.names = FALSE)
message(sprintf("Total: %.0f minutes",
                sum(vapply(out, function(o) o$seconds, numeric(1))) / 60))
