# Sensitivity of the simulation results of Duncan (2026) to the prior of the
# models with a constant hazard rate (BDRIBS and EXP_BERN):
#   Duncan B (2026). Robust safety monitoring and signal detection using
#   alternatives to the standard Poisson distribution. Journal of
#   Biopharmaceutical Statistics 36(2):227-244.
#   doi:10.1080/10543406.2024.2395532
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
# Published values: Duncan2026_digitized.csv. Tolerance: four Monte Carlo
# standard errors plus the reading error of the figures, taken as a fifth of
# the grid interval (0.02 for the power axes with grid lines every 10% and
# 0.004 for the risk difference axis of Figure 6 with grid lines every 2%).
# The note of each item gives the value of the main simulation.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Duncan2026_prior.R")
# Output: inst/reproduce/results/reproduce_Duncan2026_prior.csv
# Run time with n_sim_prior = 2000: about 30 minutes on one core (EXP_BERN
# takes most of it). n_sim_prior can be set before sourcing this script.

if (!exists("n_sim_prior")) n_sim_prior <- 2000L

source("inst/reproduce/judge.R")
explained <- read_explained()
dig <- utils::read.csv("inst/reproduce/Duncan2026_digitized.csv",
                       stringsAsFactors = FALSE)
main <- readRDS("inst/reproduce/results/Duncan2026_simulation.rds")
main_by_id <- stats::setNames(main$out,
                              vapply(main$out, function(o) o$id, ""))

accrual <- 26
rate_of <- function(ip, t) -log(1 - ip) / t
cutoff <- 0.95
read_power <- 0.02
read_rd <- 0.004

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

power_of <- function(r, model) {
  x <- r$prob[r$model == model]
  x <- x[!is.na(x)]
  p <- mean(x > cutoff)
  c(value = p, se = sqrt(p * (1 - p) / length(x)))
}
mean_of <- function(x) {
  x <- x[!is.na(x)]
  c(value = mean(x), se = stats::sd(x) / sqrt(length(x)))
}
rd_of <- function(r, model) mean_of(r$ip_difference[r$model == model])
obs_of <- function(r) mean_of(r$observed_difference[r$model == r$model[1]])
main_res <- function(id) {
  r <- main_by_id[[id]]$results
  r[r$look == 1, ]
}
published <- function(figure, point, model, quantity) {
  v <- dig$value[dig$figure == figure & dig$point == point &
                   dig$look == 1 & dig$model == model &
                   dig$quantity == quantity]
  stopifnot(length(v) == 1L)
  v
}

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}
note_of <- function(v, v0, digits) {
  sprintf(paste0("Monte Carlo SE %.4f; main simulation %.", digits, "f"),
          v[["se"]], v0[["value"]])
}

for (j in seq_along(tasks)) {
  tk <- tasks[[j]]
  r <- out[[j]]
  m0 <- main_res(tk$main_id)
  src <- sprintf("Figure %s (prior: %s)", tk$figure, tk$prior_label)
  for (m in tk$models) {
    v <- power_of(r, m)
    add(sprintf("D_PRIOR_%s_P%d_%s_power", tk$figure, tk$point, m), src,
        sprintf("power of %s (point %d, threshold %g)", m, tk$point, cutoff),
        published(tk$figure, tk$point, m, "power"), v[["value"]],
        4 * v[["se"]] + read_power, note_of(v, power_of(m0, m), 3))
  }
  if (tk$figure == "5") {
    src6 <- sprintf("Figure 6 (prior: %s)", tk$prior_label)
    for (m in tk$models) {
      v <- rd_of(r, m)
      add(sprintf("D_PRIOR_6_P%d_%s_rd", tk$point, m), src6,
          sprintf("posterior mean of the risk difference of %s (point %d)",
                  m, tk$point),
          published("6", tk$point, m, "rd"), v[["value"]],
          4 * v[["se"]] + read_rd, note_of(v, rd_of(m0, m), 4))
    }
    v <- obs_of(r)
    add(sprintf("D_PRIOR_6_P%d_SIM_expected_rd", tk$point), src6,
        sprintf("risk difference in the simulated data (point %d)",
                tk$point),
        published("6", tk$point, "SIM", "expected_rd"), v[["value"]],
        4 * v[["se"]] + read_rd, note_of(v, obs_of(m0), 4))
  }
}

# statements D_Q_F5a and D_Q_F6a of reproduce_Duncan2026_figures.R with the
# prior of this script
statement <- function(id, source, text, ok) {
  add(id, source, text, 1, as.numeric(isTRUE(ok)), 0.5)
}
res5 <- function(point) out[[4L + point]]
pw5 <- function(point, m) power_of(res5(point), m)[["value"]]
rd5 <- function(point, m) rd_of(res5(point), m)[["value"]]
ob5 <- function(point) obs_of(res5(point))[["value"]]
statement("D_PRIOR_Q_F5a", "Section 4.2.1, Figure 5 (prior: 4% at week 104)",
          paste("decreasing hazard: the Type I errors of EXP_BERN and",
                "BDRIBS(2) are inflated (> 0.10), criterion of D_Q_F5a"),
          pw5(1, "EXP_BERN") > 0.10 && pw5(1, "BDRIBS2") > 0.10)
statement("D_PRIOR_Q_F6a", "Section 4.2.1, Figure 6 (prior: 4% at week 104)",
          paste("decreasing hazard: EXP_BERN and BDRIBS(2) overestimate the",
                "risk difference (by >= 0.01), criterion of D_Q_F6a"),
          all(vapply(1:3, function(p) {
            min(rd5(p, "EXP_BERN"), rd5(p, "BDRIBS2")) - ob5(p) >= 0.01
          }, logical(1))))

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Duncan2026_prior")
self_check(tab)
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Duncan2026_prior.csv",
                 row.names = FALSE)
