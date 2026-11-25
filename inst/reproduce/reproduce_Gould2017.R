# Reproduction of the numerical results of
#   Gould AL, Wang WB (2017). Monitoring potential adverse event rate
#   differences using data from blinded trials: the canary in the coal mine.
#   Statistics in Medicine 36(1):92-104. doi:10.1002/sim.7129
# with elicit_prior(), signal_count(), canary_post() and count_oc() of the
# bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Gould2017.R")
# Output: inst/reproduce/results/reproduce_Gould2017.csv
# Run time: about ten minutes (the critical counts for N = 1000).
#
# Table I: parameters of the Gamma prior from the median and an upper
# quantile (two decimals, tolerance 0.005, or 0.0005 for three decimals).
# Table III: critical counts X_crit = min(X: P(R > R_crit | X) > gamma) for
# the binomial model with the ratio metric. The lognormal prior of R has mean
# 1 and probability 0.95 between 0.1 and 8, the default of the authors' code
# (get.Xcrit.fn with ML = 0.1 and MU = 8; the text says 10). The counts are
# integers, so the tolerance is 0.5 (equal counts).
# Table V: critical counts of the Poisson model with the Gamma(75.71, 0.85)
# prior of the expected number of control events (10,000 patients per arm,
# so 20,000 patient-years in bssd and a Gamma(75.71, 8500) prior of the rate)
# and the default lognormal prior (0.1, 10) of the authors' Poisson code; and
# the probabilities that X exceeds the published X_crit under the true
# control counts 70, 90, 110 and rate ratios 1, 2, 4. The table does not say
# whether "exceeding" means X >= X_crit or X > X_crit; the exact Poisson
# probability of the two that is closer to the published value is compared
# with tolerance 0.005 (two decimals), and the note says which one.
#
# The published critical counts were computed by the authors' code on fixed
# grids. A transcription of the code (Rpackage/tools/external/
# gould2017_author_grid.py, outside the package because the code has no
# stated license; its counts are in inst/validation/python/
# ref_gould2017_author_grid.csv) reproduces all of them; differences from
# the exact counts of bssd are registered in explained.csv.

source("inst/reproduce/judge.R")
explained <- read_explained()

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

# Table I ------------------------------------------------------------------
t1 <- data.frame(
  ratio = rep(c(2, 4), each = 3),
  zeta = rep(c(0.8, 0.9, 0.95), 2),
  a = c(1.42, 2.96, 4.69, 0.45, 0.77, 1.11),
  b1 = c(1.11, 2.63, 4.36, 0.182, 0.47, 0.8),
  b2 = c(0.55, 1.32, 2.18, 0.091, 0.24, 0.4),
  b5 = c(0.22, 0.53, 0.87, 0.036, 0.094, 0.16))
dec_tol <- function(x) {
  # half a unit of the last published decimal
  d <- nchar(sub("^[^.]*\\.?", "", format(x, scientific = FALSE)))
  0.5 * 10^(-d)
}
for (i in seq_len(nrow(t1))) {
  for (med in c(1, 2, 5)) {
    p <- elicit_prior(c(med, t1$ratio[i] * med), c(0.5, t1$zeta[i]),
                      "gamma")
    lab <- sprintf("median %g, ratio %g, zeta %g", med, t1$ratio[i],
                   t1$zeta[i])
    if (med == 1) {
      add(sprintf("G_T1_a_%d", i), "Table I", paste("a,", lab), t1$a[i],
          p[["shape"]], dec_tol(t1$a[i]))
    }
    b_pub <- t1[[paste0("b", med)]][i]
    add(sprintf("G_T1_b%d_%d", med, i), "Table I", paste("b,", lab), b_pub,
        p[["rate"]], dec_tol(b_pub))
  }
}

# Table III -----------------------------------------------------------------
mp8 <- elicit_prior(c(0.1, 8), 0.95, "lnorm", mean = 1)
t3 <- list(
  list(ab = c(1, 49), x = c(12, 11, 18, 90, 90, 144)),
  list(ab = c(20, 980), x = c(6, 6, 8, 31, 36, 44)),
  list(ab = c(5, 45), x = c(21, 23, 31, 182, 206, 275)),
  list(ab = c(100, 900), x = c(16, 18, 22, 118, 157, 173)))
rules <- list(c(1, 0.8), c(2, 0.5), c(2, 0.8))
for (s in seq_along(t3)) {
  j <- 0L
  for (nn in c(100, 1000)) {
    for (rl in rules) {
      j <- j + 1L
      sc <- signal_count(n = nn, method = "canary", k = 1,
                         prior_c = t3[[s]]$ab, metric_prior = mp8,
                         crit = rl[1], cutoff = rl[2])
      add(sprintf("G_T3_S%d_%d", s, j), "Table III",
          sprintf("X_crit, prior (%g, %g), N = %d, R_crit = %g, gamma = %g",
                  t3[[s]]$ab[1], t3[[s]]$ab[2], nn, rl[1], rl[2]),
          t3[[s]]$x[j], sc$y_crit, 0.5)
    }
  }
}

# Table V -------------------------------------------------------------------
mp10 <- elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)
t5_rules <- list(c(1, 0.8), c(2, 0.5), c(2, 0.8), c(4, 0.5))
t5_x <- c(215, 297, 316, 498)
for (j in seq_along(t5_rules)) {
  sc <- signal_count(n = 20000, method = "canary", k = 1,
                     prior_c = c(75.71, 8500), likelihood = "poisson",
                     metric_prior = mp10, crit = t5_rules[[j]][1],
                     cutoff = t5_rules[[j]][2], y_max = 2000)
  add(sprintf("G_T5_X%d", j), "Table V",
      sprintf("X_crit, Poisson, R_crit = %g, gamma = %g", t5_rules[[j]][1],
              t5_rules[[j]][2]), t5_x[j], sc$y_crit, 0.5)
}
t5_p <- list(
  "70" = list("1" = c(0, 0, 0, 0), "2" = c(0.35, 0, 0, 0),
              "4" = c(1, 1, 0.97, 0)),
  "90" = list("1" = c(0.01, 0, 0, 0), "2" = c(1, 0.05, 0, 0),
              "4" = c(1, 1, 1, 0.01)),
  "110" = list("1" = c(0.62, 0, 0, 0), "2" = c(1, 0.97, 0.77, 0),
               "4" = c(1, 1, 1, 0.99)))
for (lc in names(t5_p)) {
  for (rr in names(t5_p[[lc]])) {
    # X ~ Poisson((1 + R) lambda_C) with lambda_C the expected control count:
    # exposure 20000 and control rate lambda_C / 10000
    p_ge <- count_oc(t5_x, 20000, rate_c = as.numeric(lc) / 10000,
                     rate_d = as.numeric(rr) * as.numeric(lc) / 10000,
                     likelihood = "poisson")
    p_gt <- count_oc(t5_x + 1, 20000, rate_c = as.numeric(lc) / 10000,
                     rate_d = as.numeric(rr) * as.numeric(lc) / 10000,
                     likelihood = "poisson")
    for (j in seq_along(t5_x)) {
      pub <- t5_p[[lc]][[rr]][j]
      use_ge <- abs(p_ge[j] - pub) <= abs(p_gt[j] - pub)
      add(sprintf("G_T5_P%s_%s_%d", lc, rr, j), "Table V",
          sprintf("P(X exceeds %d), true lambda_C = %s, R = %s", t5_x[j], lc,
                  rr), pub, if (use_ge) p_ge[j] else p_gt[j], 0.005,
          sprintf("P(X >= X_crit) = %.4f, P(X > X_crit) = %.4f; compared: %s",
                  p_ge[j], p_gt[j], if (use_ge) ">=" else ">"))
    }
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Gould2017")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Gould2017.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
