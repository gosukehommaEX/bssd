# Reproduction of the numerical results of
#   Mukhopadhyay S, Waterhouse B, Hartford A (2018). Bayesian detection of
#   potential risk using inference on blinded safety data. Pharmaceutical
#   Statistics 17(6):823-834. doi:10.1002/pst.1898
# with bdribs_post() of the bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Mukhopadhyay2018.R")
# Output: inst/reproduce/results/reproduce_Mukhopadhyay2018.csv
# Run time: about one minute.
#
# The published values were computed by MCMC with the bdribs package; the
# tolerance for posterior quantities is 0.01 (0.02 for the posterior mean),
# which allows for Monte Carlo error. All examples use the screening analysis
# with a known background rate (pt_hp = Inf) and the uniform prior for p
# (p_prior = c(1, 1)), which is the default of the bdribs package.

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
post_prob <- function(y, exposure, k, rate, c, level = 0.95,
                      details = FALSE) {
  bdribs_post(y = y, exposure = exposure, k = k, lambda_hp = rate,
              pt_hp = Inf, p_prior = c(1, 1), rr0 = c, level = level,
              details = details)
}

# Table 1: limits on the posterior probability scale that correspond to
# Bayes factors 1, 3, 10, 30 and 100. The prior probability Q of r > c is
# obtained from bdribs_post() with no events and negligible exposure (1e-30),
# where the posterior equals the prior; the limit for Bayes factor B is
# B * O / (1 + B * O) with O = Q / (1 - Q).
table1 <- list(
  list(prior = c(1, 1), k = 1, c = 1.0, v = c(0.50, 0.75, 0.91, 0.97, 0.99)),
  list(prior = c(1, 1), k = 1, c = 1.2, v = c(0.45, 0.71, 0.89, 0.96, 0.99)),
  list(prior = c(1, 1), k = 1, c = 1.5, v = c(0.40, 0.67, 0.87, 0.95, 0.99)),
  list(prior = c(1, 1), k = 2, c = 1.0, v = c(0.33, 0.60, 0.83, 0.94, 0.98)),
  list(prior = c(1, 1), k = 2, c = 1.2, v = c(0.29, 0.56, 0.81, 0.93, 0.98)),
  list(prior = c(1, 1), k = 2, c = 1.5, v = c(0.25, 0.50, 0.77, 0.91, 0.97)),
  list(prior = c(0.5, 0.5), k = 1, c = 1.0,
       v = c(0.50, 0.75, 0.91, 0.97, 0.99)),
  list(prior = c(0.5, 0.5), k = 1, c = 1.2,
       v = c(0.47, 0.73, 0.90, 0.96, 0.99)),
  list(prior = c(0.5, 0.5), k = 1, c = 1.5,
       v = c(0.44, 0.70, 0.89, 0.96, 0.99)),
  list(prior = c(0.5, 0.5), k = 2, c = 1.0,
       v = c(0.39, 0.66, 0.87, 0.95, 0.98)),
  list(prior = c(0.5, 0.5), k = 2, c = 1.2,
       v = c(0.36, 0.63, 0.85, 0.95, 0.98)),
  list(prior = c(0.5, 0.5), k = 2, c = 1.5,
       v = c(0.33, 0.60, 0.83, 0.94, 0.98))
)
bf <- c(1, 3, 10, 30, 100)
for (i in seq_along(table1)) {
  s <- table1[[i]]
  q <- bdribs_post(y = 0, exposure = 1e-30, k = s$k, lambda_hp = 1,
                   pt_hp = Inf, p_prior = s$prior, rr0 = s$c,
                   details = FALSE)$prob
  o <- q / (1 - q)
  lim <- bf * o / (1 + bf * o)
  for (j in seq_along(bf)) {
    add(sprintf("M_T1_%02d_%d", i, j), "Table 1",
        sprintf("limit for BF = %g, p ~ Beta(%g, %g), k = %g, c = %g",
                bf[j], s$prior[1], s$prior[2], s$k, s$c),
        s$v[j], lim[j], 0.005 + 1e-9, "published values rounded to 0.01")
  }
}

# Section 3.4 and 4.4: the boundary 3 / (3 + k c) of Equation (1)
add("M_EQ1_1", "Section 4.4", "3 / (3 + k c), k = 2, c = 1", 0.60,
    {
      q <- bdribs_post(y = 0, exposure = 1e-30, k = 2, lambda_hp = 1,
                       pt_hp = Inf, p_prior = c(1, 1), rr0 = 1,
                       details = FALSE)$prob
      3 * (q / (1 - q)) / (1 + 3 * q / (1 - q))
    }, 0.005)
add("M_EQ1_2", "Section 4.4", "3 / (3 + k c), k = 2, c = 1.3", 0.536,
    {
      q <- bdribs_post(y = 0, exposure = 1e-30, k = 2, lambda_hp = 1,
                       pt_hp = Inf, p_prior = c(1, 1), rr0 = 1.3,
                       details = FALSE)$prob
      3 * (q / (1 - q)) / (1 + 3 * q / (1 - q))
    }, 0.0005)

# Section 4.3, hypothetical example 1: k = 1, background rate 0.0045, c = 1.2
add("M_EX1_1", "Section 4.3",
    "P(r > 1.2 | y = 4, E = 570), read from the contour plot", 0.55,
    post_prob(4, 570, 1, 0.0045, 1.2)$prob, 0.05,
    "published value is read from Figure 2 (about 0.55)")
add("M_EX1_2", "Section 4.3", "P(r > 1.2 | y = 7, E = 800)", 0.80,
    post_prob(7, 800, 1, 0.0045, 1.2)$prob, 0.01)
ci1 <- qgamma(c(0.05, 0.95), 18.5, 4000.5)
add("M_EX1_3", "Section 4.3", "90% interval of Gamma(18.5, 4000.5), lower",
    0.0030, ci1[1], 0.00005)
add("M_EX1_4", "Section 4.3", "90% interval of Gamma(18.5, 4000.5), upper",
    0.0065, ci1[2], 0.00005)
add("M_EX1_5", "Section 4.3",
    "P(r > 1.2 | y = 7, E = 800), background rate 0.0065", 0.52,
    post_prob(7, 800, 1, 0.0065, 1.2)$prob, 0.01)
add("M_EX1_6", "Section 4.3",
    "P(r > 1.2 | y = 7, E = 800), background rate 0.0030", 0.95,
    post_prob(7, 800, 1, 0.0030, 1.2)$prob, 0.01)
add("M_EX1_7", "Section 4.3",
    "upper end of the range quoted later in the same paragraph", 0.92,
    post_prob(7, 800, 1, 0.0030, 1.2)$prob, 0.01)

# Section 4.4, hypothetical example 2: k = 2, background rate 0.011, E = 540
add("M_EX2_1", "Section 4.4", "P(r > 1 | y = 8)", 0.622,
    post_prob(8, 540, 2, 0.011, 1)$prob, 0.01)
add("M_EX2_2", "Section 4.4", "P(r > 1.3 | y = 8)", 0.434,
    post_prob(8, 540, 2, 0.011, 1.3)$prob, 0.01)
add("M_EX2_3", "Section 4.4", "P(r > 1 | y = 9)", 0.756,
    post_prob(9, 540, 2, 0.011, 1)$prob, 0.01)
add("M_EX2_4", "Section 4.4", "P(r > 1.3 | y = 9)", 0.584,
    post_prob(9, 540, 2, 0.011, 1.3)$prob, 0.01)
f95 <- post_prob(9, 540, 2, 0.011, 1, level = 0.95, details = TRUE)
f90 <- post_prob(9, 540, 2, 0.011, 1, level = 0.90, details = TRUE)
add("M_EX2_5", "Section 4.4", "posterior mean of r | y = 9", 1.53,
    unname(f95$rr["mean"]), 0.02)
add("M_EX2_6", "Section 4.4", "95% interval of r | y = 9, lower", 0.52,
    unname(f95$rr["lower"]), 0.02)
add("M_EX2_7", "Section 4.4", "95% interval of r | y = 9, upper", 2.83,
    unname(f95$rr["upper"]), 0.02)
add("M_EX2_8", "Section 4.4",
    "published '95%' interval compared with the 90% interval, lower", 0.52,
    unname(f90$rr["lower"]), 0.02)
add("M_EX2_9", "Section 4.4",
    "published '95%' interval compared with the 90% interval, upper", 2.83,
    unname(f90$rr["upper"]), 0.02)
ci2 <- qgamma(c(0.05, 0.95), 33.5, 3000.5)
add("M_EX2_10", "Section 4.4", "90% interval of Gamma(33.5, 3000.5), lower",
    0.0082, ci2[1], 0.00005)
add("M_EX2_11", "Section 4.4", "90% interval of Gamma(33.5, 3000.5), upper",
    0.0145, ci2[2], 0.00005)
add("M_EX2_12", "Section 4.4", "P(r > 1 | y = 9), background rate 0.0145",
    0.48, post_prob(9, 540, 2, 0.0145, 1)$prob, 0.01)
add("M_EX2_13", "Section 4.4", "P(r > 1 | y = 9), background rate 0.0082",
    0.92, post_prob(9, 540, 2, 0.0082, 1)$prob, 0.01)
add("M_EX2_14", "Section 4.4", "P(r > 1.3 | y = 9), background rate 0.0145",
    0.28, post_prob(9, 540, 2, 0.0145, 1.3)$prob, 0.01)
add("M_EX2_15", "Section 4.4", "P(r > 1.3 | y = 9), background rate 0.0082",
    0.83, post_prob(9, 540, 2, 0.0082, 1.3)$prob, 0.01)

# Section 5, operating characteristics: true rate 0.005, true r = 2 (null
# r = 1), c = 1.2, E = 1200, k = 2, background rate 0.0045, threshold 0.60.
# Y is Poisson with mean E * rate * (r k + 1) / (k + 1), and P(r > c | y)
# increases with y, so the rates of the decision rule are Poisson tails.
pr <- vapply(0:40, function(y) post_prob(y, 1200, 2, 0.0045, 1.2)$prob,
             numeric(1))
y_star <- min(which(pr > 0.60)) - 1L
mean_null <- 1200 * 0.005 * (1 * 2 + 1) / 3
mean_alt <- 1200 * 0.005 * (2 * 2 + 1) / 3
add("M_S5_1", "Section 5", "true positive rate for threshold 0.60", 0.67,
    ppois(y_star - 1, mean_alt, lower.tail = FALSE), 0.01)
add("M_S5_2", "Section 5",
    "false positive rate for threshold 0.60 ('within 15%')", 0.15,
    ppois(y_star - 1, mean_null, lower.tail = FALSE), 0.0025)
add("M_S5_3", "Section 5", "smallest number of events with a signal", NA,
    y_star, NA, "events needed for P(r > 1.2 | y) > 0.60")

tab <- judge(do.call(rbind, rows), explained)
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Mukhopadhyay2018.csv",
                 row.names = FALSE)
print_counts(tab, "reproduce_Mukhopadhyay2018")
self_check(tab)
