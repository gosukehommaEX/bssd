# Reproduction of the numerical results of
#   Goren E, Lin L-A, Blaustein RO, Ball G (2020). Bayesian meta-analysis of
#   safety outcomes using blinded clinical trial data. Therapeutic
#   Innovation & Regulatory Science 54:1557-1565.
#   doi:10.1007/s43441-020-00183-7
# with multibin_post() of the bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Goren2020.R")
# Table 2 (simulation) is judged by reproduce_Goren2020_table2.R.
# Output: inst/reproduce/results/reproduce_Goren2020.csv
# Run time: seconds.
#
# Table 4 (real data of eight studies, Table 3): the posterior probabilities
# that the rate ratio exceeds 1, 2 and 3 with homogeneous control priors
# (20 events among 1687 control patients: Beta(20, 1667)) and heterogeneous
# priors (studies in Japan: mean 0.009 and prior sample size 219, that is
# Beta(2, 217); others: mean 0.012 and 1468, that is Beta(18, 1450)), with all
# studies and without study 5. The published values are MCMC estimates with
# standard errors; the tolerance is three standard errors (at least 0.0005)
# plus 0.0005 for rounding.

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

# Table 4 ---------------------------------------------------------------
y <- c(8, 4, 5, 6, 20, 2, 7, 39)
n <- c(354, 305, 307, 582, 244, 68, 307, 1612)
k <- c(2, 2, 1, 1, 2, 1, 2, 1)
japan <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)
hom <- matrix(c(20, 1667), 8, 2, byrow = TRUE)
het <- t(sapply(japan, function(j) if (j) c(2, 217) else c(18, 1450)))
pub <- list(
  full_hom = c(1.000, 0.984, 0.462), full_hom_se = c(0, 0.001, 0.004),
  sens_hom = c(0.998, 0.561, 0.034), sens_hom_se = c(0, 0.004, 0.001),
  full_het = c(1.000, 0.836, 0.158), full_het_se = c(0, 0.003, 0.003),
  sens_het = c(0.994, 0.486, 0.031), sens_het_se = c(0.001, 0.004, 0.001))
for (an in c("full", "sens")) {
  use <- if (an == "full") seq_len(8) else setdiff(seq_len(8), 5)
  for (pr in c("hom", "het")) {
    pc <- if (pr == "hom") hom[use, , drop = FALSE] else het[use, , drop = FALSE]
    key <- paste(an, pr, sep = "_")
    for (j in 1:3) {
      fit <- multibin_post(y[use], n[use], k = k[use], prior_c = pc,
                           rr0 = j, details = FALSE)
      se <- max(pub[[paste0(key, "_se")]][j], 0.0005)
      add(sprintf("GO_T4_%s_%d", key, j), "Table 4",
          sprintf("P(delta > %d), %s analysis, %s priors", j,
                  if (an == "full") "full" else "sensitivity",
                  if (pr == "hom") "homogeneous" else "heterogeneous"),
          pub[[key]][j], fit$prob, 3 * se + 0.0005)
    }
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Goren2020")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Goren2020.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
