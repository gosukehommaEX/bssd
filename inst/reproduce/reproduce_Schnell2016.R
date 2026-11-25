# Reproduction of the numerical results of
#   Schnell PM, Ball G (2016). A Bayesian exposure-time method for clinical
#   trial safety monitoring with blinded data. Therapeutic Innovation &
#   Regulatory Science 50(6):833-838. doi:10.1177/2168479016656702
# with elicit_prior() and exptime_post() of the bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Schnell2016.R")
# Table 1 (simulation) is judged by reproduce_Schnell2016_table1.R.
# Output: inst/reproduce/results/reproduce_Schnell2016.csv
# Run time: about one minute.
#
# The simulated case study (Table 2) is generated with the data generation
# of the authors' code (schnell2016_data.R) and seed 1. The paper gives the
# true rates 0.011 (control) and 0.016 (treatment) per person-month; the
# authors' script for the example (example/3-data-analysis.R) uses 0.02 for
# the treatment. Both are generated; the counts of the published table
# decide which one is compared (the other is INFO).
#
# The posterior summaries of the paper are Monte Carlo estimates from 1000
# draws of the authors' Gibbs sampler, whose Metropolis step for the
# treatment assignments uses Poisson probabilities of the arm totals
# instead of the subject-level likelihood of the model (see
# Rpackage/tools/external/schnell2016_author_code.R, outside the package).
# bssd computes the posterior of the subject-level model exactly. The
# tolerances are 0.001 for the posterior means (four decimals published) and
# 0.03 for the probabilities (about three Monte Carlo standard errors of a
# proportion from 1000 dependent draws). Counts are integers and are compared
# with the tolerance 0.5 (equal counts), so that judge.R's self-check (a shift
# of five tolerances) detects a change.

source("inst/reproduce/judge.R")
source("inst/reproduce/schnell2016_data.R")
explained <- read_explained()

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

# Prior elicitation example and the priors of the case study ---------------
pe <- elicit_prior(c(0.5, 2), c(0.5, 0.75), "gamma")
add("S_PE_1", "Prior Elicitation", "shape of the elicited Gamma prior",
    0.346, pe[["shape"]], 0.0005)
add("S_PE_2", "Prior Elicitation", "rate of the elicited Gamma prior",
    0.209, pe[["rate"]], 0.0005)
prior_c <- c(339, 26064)
prior_d <- c(339 / 100 + 4, 26064 / 100 + 167)
critical <- 1.2 * 0.013
add("S_T2_P0_MC", "Table 2", "prior mean of the control rate", 0.0130,
    prior_c[1] / prior_c[2], 0.00005)
add("S_T2_P0_MD", "Table 2", "prior mean of the treatment rate", 0.0173,
    prior_d[1] / prior_d[2], 0.00005)
add("S_T2_P0_PR", "Table 2", "prior P(treatment rate > 0.0156)", 0.551,
    pgamma(critical, prior_d[1], prior_d[2], lower.tail = FALSE), 0.0005)

# Table 2: interim analyses of the simulated case study --------------------
interim <- c(4, 8, 12, 16)
pub <- list(enrolled = c(998, 1500, 1500, 1500),
            follow = c(2.1, 7.6, 13.6, 19.6),
            events = c(34, 114, 222, 336),
            mean_c = c(0.0131, 0.0131, 0.0133, 0.0136),
            mean_d = c(0.0175, 0.0158, 0.0172, 0.0179),
            prob = c(0.710, 0.518, 0.872, 0.968))
for (rate_d in c(0.016, 0.02)) {
  dat <- schnell_generate(1, 1500, 2 / 3, 6, 18, 0.011, rate_d)
  snap <- lapply(interim, function(tm) schnell_data_at(dat, tm))
  counts_ok <- all(vapply(snap, function(d) sum(d$y), numeric(1)) ==
                     pub$events)
  tol <- function(x) if (counts_ok) x else NA
  tag <- if (rate_d == 0.016) "A" else "B"
  lab <- sprintf("data generated with treatment rate %g", rate_d)
  note <- if (counts_ok) "" else
    "the event counts of these data differ from Table 2 (INFO)"
  for (j in seq_along(interim)) {
    d <- snap[[j]]
    add(sprintf("S_T2_%s_N%d", tag, j), "Table 2",
        sprintf("patients enrolled at month %d, %s", interim[j], lab),
        pub$enrolled[j], nrow(d), tol(0.5), note)
    add(sprintf("S_T2_%s_F%d", tag, j), "Table 2",
        sprintf("total follow-up (1000 person-months), month %d, %s",
                interim[j], lab), pub$follow[j], sum(d$t) / 1000, tol(0.05),
        note)
    add(sprintf("S_T2_%s_E%d", tag, j), "Table 2",
        sprintf("total SAE at month %d, %s", interim[j], lab),
        pub$events[j], sum(d$y), tol(0.5), note)
    if (counts_ok) {
      fit <- exptime_post(d$y, d$t, k = 2, prior_c = prior_c,
                          prior_d = prior_d, rate0 = critical,
                          cutoff = 0.62)
      add(sprintf("S_T2_%s_MC%d", tag, j), "Table 2",
          sprintf("posterior mean of the control rate, month %d",
                  interim[j]), pub$mean_c[j],
          fit$estimates[["lambda_c"]], 0.001)
      add(sprintf("S_T2_%s_MD%d", tag, j), "Table 2",
          sprintf("posterior mean of the treatment rate, month %d",
                  interim[j]), pub$mean_d[j],
          fit$estimates[["lambda_d"]], 0.001)
      add(sprintf("S_T2_%s_PR%d", tag, j), "Table 2",
          sprintf("P(treatment rate > 0.0156), month %d", interim[j]),
          pub$prob[j], fit$prob, 0.03)
      add(sprintf("S_T2_%s_ER%d", tag, j), "Table 2",
          sprintf("integration error estimate of bssd, month %d",
                  interim[j]), NA, fit$integration$error, NA)
    }
  }
  if (counts_ok) {
    # end of the study: every subject followed for 18 months
    add(sprintf("S_T2_%s_END_F", tag), "Simulated Case Study",
        "total follow-up at the end of the study (1000 person-months)", 27,
        1500 * 18 / 1000, 0.5)
    add(sprintf("S_T2_%s_END_E", tag), "Simulated Case Study",
        "total SAE at the end of the study", 490, length(dat$ev_id), 0.5)
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Schnell2016")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Schnell2016.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
