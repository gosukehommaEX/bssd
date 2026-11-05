# Judgements for the sensitivity of the simulation results of Duncan (2026)
# to the prior of the models with a constant hazard rate (BDRIBS and
# EXP_BERN), from the results of reproduce_Duncan2026_prior_simulation.R.
# The priors and the reason for this script are described in the header of
# that script.
#
# Published values: Duncan2026_digitized.csv. Tolerance: four Monte Carlo
# standard errors plus the reading error of the figures, taken as a fifth of
# the grid interval (0.02 for the power axes with grid lines every 10% and
# 0.004 for the risk difference axis of Figure 6 with grid lines every 2%).
# The note of each item gives the value of the main simulation
# (reproduce_Duncan2026_simulation.R).
#
# Run from the package root after the simulation script:
#   source("inst/reproduce/reproduce_Duncan2026_prior.R")
# Output: inst/reproduce/results/reproduce_Duncan2026_prior.csv

source("inst/reproduce/judge.R")
explained <- read_explained()
dig <- utils::read.csv("inst/reproduce/Duncan2026_digitized.csv",
                       stringsAsFactors = FALSE)
main <- readRDS("inst/reproduce/results/Duncan2026_simulation.rds")
main_by_id <- stats::setNames(main$out,
                              vapply(main$out, function(o) o$id, ""))
prior <- readRDS("inst/reproduce/results/Duncan2026_prior.rds")
tasks <- prior$tasks
out <- prior$out

cutoff <- 0.95
read_power <- 0.02
read_rd <- 0.004

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
# priors of reproduce_Duncan2026_prior_simulation.R
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
