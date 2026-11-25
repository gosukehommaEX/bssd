# Reproduction of the numerical results of
#   Waterhouse B, Hartford A, Mukhopadhyay S, Ferguson R, Hendrickson BA
#   (2022). Using the Bayesian detection of potential risk using inference
#   on blinded safety data (BDRIBS) method to support the decision to refer
#   an event for unblinded evaluation. Pharmaceutical Statistics
#   21(2):372-385. doi:10.1002/pst.2175
# with bdribs_post(), bdribs_grid() and pooled_ratio() of the bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Waterhouse2022.R")
# Output: inst/reproduce/results/reproduce_Waterhouse2022.csv
# Run time: seconds.
#
# The examples use a fixed background rate (pt_hp = Inf) and "a
# non-informative prior on p centered on 0.5"; the Jeffreys prior Beta(0.5,
# 0.5), the default named in the documentation of the bdribs package, is
# used (Beta(1, 1) gives probabilities that differ by up to 0.02). The
# published values were read from an R Shiny application that samples by
# MCMC, so the tolerance is 0.01 for values given in percent and 0.02 for
# values given as "about" or "approximately". Statements of the form "greater
# than" are checked as indicators (published 1, tolerance 0.5).

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
prob <- function(y, exposure, k, rate) {
  bdribs_post(y = y, exposure = exposure, k = k, lambda_hp = rate,
              pt_hp = Inf, p_prior = c(0.5, 0.5), details = FALSE)$prob
}

# Example 1: MACE, 3.1 events per 100 patient-years, randomization 1:1
add("W_EX1_1", "Example 1", "P(RR > 1), 6 events, 80 patient-years", 0.90,
    prob(6, 80, 1, 0.031), 0.02, "published as about 90%")
cross6 <- uniroot(function(e) prob(6, e, 1, 0.031) - 0.80, c(80, 200),
                  tol = 1e-8)$root
add("W_EX1_2", "Example 1",
    "patient-years at which P(RR > 1) falls to 0.80 with 6 events", 100,
    cross6, 5,
    "published as 'if the exposure had been more than 100 PYs, this probability would be less than 80%'")
cross7 <- uniroot(function(e) prob(7, e, 1, 0.031) - 0.80, c(80, 200),
                  tol = 1e-8)$root
add("W_EX1_3", "Example 1",
    "patient-years at which P(RR > 1) falls to 0.80 with 7 events", 125,
    cross7, 5, "published as 'until after 125 patient-years'")
add("W_EX1_4", "Example 1", "P(RR > 1), 6 events, 125 patient-years", 0.60,
    prob(6, 125, 1, 0.031), NA, "published as 'just above 60%' (INFO)")
add("W_EX1_5", "Example 1, Figure 4",
    "P(RR > 1), 7 events, 140 patient-years", 0.74, prob(7, 140, 1, 0.031),
    0.01)
add("W_EX1_6", "Example 1",
    "indicator P(RR > 1) > 0.80, 10 events, 200 patient-years", 1,
    as.numeric(prob(10, 200, 1, 0.031) > 0.80), 0.5)
# Figure 6: the 90% credible interval of the background rate from 496
# events in 15,730 patient-years (Figure 5b), and the probability at its
# limits
bg <- qgamma(c(0.05, 0.95), 496, 15730)
add("W_EX1_7", "Example 1, Figure 6",
    "P(RR > 1), 10 events, 200 PY, lower limit of the background rate",
    0.86, prob(10, 200, 1, bg[1]), 0.02, "published as about 86%")
add("W_EX1_8", "Example 1, Figure 6",
    "P(RR > 1), 10 events, 200 PY, upper limit of the background rate",
    0.76, prob(10, 200, 1, bg[2]), 0.02, "published as about 76%")

# Example 2: urinary tract infection, 3.7 per 100 PY, randomization 2:1
add("W_EX2_1", "Example 2, Figure 7",
    "P(RR > 1), 22 events, 504 patient-years, k = 2", 0.70,
    prob(22, 504, 2, 0.037), 0.02, "published as approximately 0.70")

# Example 3: pooled trials A (250 PY, 1:1) and B (80 PY, 2:1)
kp <- pooled_ratio(c(250, 80), c(1, 2))
add("W_EX3_1", "Example 3", "pooled randomization ratio", 1.18, kp, 0.005)
add("W_EX3_2", "Example 3", "P(RR > 1), 16 events, 330 PY, pooled ratio",
    0.90, prob(16, 330, kp, 0.031), 0.01)
add("W_EX3_3", "Example 3, Figure 9",
    "P(RR > 1), 16 events, 330 PY, lower limit of the background rate",
    0.85, prob(16, 330, kp, bg[2]), 0.02,
    paste("published as 'about 0.85' for the lower bound of the",
          "sensitivity analysis, which is the upper limit of the",
          "background rate"))
add("W_EX3_4", "Example 3", "P(RR > 1) before pooling (trial A alone)",
    0.80, prob(12, 250, 1, 0.031), NA,
    paste("the probability 'had risen from 0.80'; the paper does not say",
          "which data gave 0.80 (INFO)"))

# Figure 3: the probabilities on a grid agree with single calls
g <- bdribs_grid(y = c(6, 7, 8), exposure = c(80, 100, 125, 140),
                 lambda_hp = 0.031, p_prior = c(0.5, 0.5))
add("W_F3_1", "Figure 3",
    "max |bdribs_grid - bdribs_post| on the grid of Figure 3", 0,
    max(abs(g$prob - mapply(prob, g$y, g$exposure, 1, 0.031))), 1e-12)

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Waterhouse2022")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Waterhouse2022.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
