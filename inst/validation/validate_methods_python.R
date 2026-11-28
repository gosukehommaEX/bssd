# Validation of the functions added in bssd 0.2.0 (exptime_post(),
# binbeta_post(), canary_post(), mapbern_post(), map_beta() and
# multibin_post()), and of expbern_post() for data that conflict with an
# informative prior of p, against independent Python computations
# (inst/validation/python/ref_*.py): composite Gauss-Legendre rules split
# exactly at the thresholds, an exact finite mixture for binbeta_post(), and
# nested adaptive quadrature for multibin_post(). The cases include data
# that conflict with informative priors (E3, B3 to B5, M3, G3, G4, X1 to X3),
# all subjects with the event (C5), MAP priors from historical arms with few
# events (MS1 to MS7) and control priors with a shape parameter below 1
# (BS1, ES1). The reference values are collected in python/ref_methods.csv
# and the data in python/data_methods.csv.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_methods_python.R")
# Output: inst/validation/results/validate_methods_python.csv
# Run time: about two minutes.
#
# Tolerances: 2e-6 for posterior probabilities and the MAP prior moments,
# 1e-5 relative (at least 1e-7) for posterior means, quantiles and the MAP
# prior distribution function.

source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")
ref <- utils::read.csv("inst/validation/python/ref_methods.csv",
                       stringsAsFactors = FALSE)
dat <- utils::read.csv("inst/validation/python/data_methods.csv",
                       stringsAsFactors = FALSE)
get_data <- function(name) dat[dat$dataset == name, ]

rows <- list()
add <- function(id, value) {
  r <- ref[ref$id == id, ]
  if (nrow(r) != 1L) {
    stop("Reference value not found: ", id)
  }
  tol <- if (grepl("_prob|_mean$|_sd$", id)) 2e-6 else
    max(1e-5 * abs(r$reference), 1e-7)
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = "Python", quantity = r$quantity,
    published = r$reference, reproduced = value, tolerance = tol,
    note = "", stringsAsFactors = FALSE)
}
add_q <- function(prefix, v) {
  add(paste0(prefix, "_lower"), v[["lower"]])
  add(paste0(prefix, "_median"), v[["median"]])
  add(paste0(prefix, "_upper"), v[["upper"]])
}

# EXPTIME ---------------------------------------------------------------
d <- get_data("E1")
f <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500), prior_d = c(2, 70),
                  rate0 = 0.03, rr0 = 1.2)
add("PY_EXPTIME_E1_prob_rate", f$probs[["drug"]])
add("PY_EXPTIME_E1_prob_ratio", f$probs[["ratio"]])
add("PY_EXPTIME_E1_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_E1_mean_lambda_d", f$estimates[["lambda_d"]])
add("PY_EXPTIME_E1_mean_r", f$estimates[["r"]])
add_q("PY_EXPTIME_E1_q_lambda_d", f$rate)
add_q("PY_EXPTIME_E1_q_r", f$rr)
d <- get_data("E2")
f <- exptime_post(d$y, d$t, k = 1, prior_c = c(50, 300), prior_d = c(1, 5))
add("PY_EXPTIME_E2_prob_ratio", f$prob)
add("PY_EXPTIME_E2_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_E2_mean_lambda_d", f$estimates[["lambda_d"]])
add("PY_EXPTIME_E2_mean_r", f$estimates[["r"]])
add_q("PY_EXPTIME_E2_q_lambda_d", f$rate)
add_q("PY_EXPTIME_E2_q_r", f$rr)
# E3: priors of Schnell and Ball (2016), study drug rate far in the upper
# tail of the prior of lambda_d
d <- get_data("E3")
f <- exptime_post(d$y, d$t, k = 2, prior_c = c(340, 26064),
                  prior_d = c(339 / 100 + 4, 26064 / 100 + 167),
                  rate0 = 0.065, rr0 = 5)
add("PY_EXPTIME_E3_prob_rate", f$probs[["drug"]])
add("PY_EXPTIME_E3_prob_ratio", f$probs[["ratio"]])
add("PY_EXPTIME_E3_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_E3_mean_lambda_d", f$estimates[["lambda_d"]])
add("PY_EXPTIME_E3_mean_r", f$estimates[["r"]])
add_q("PY_EXPTIME_E3_q_lambda_d", f$rate)
add_q("PY_EXPTIME_E3_q_r", f$rr)

# BINBETA -----------------------------------------------------------------
f <- binbeta_post(40, 300, k = 1, prior_c = c(31, 271), prior_d = c(1, 1),
                  p0 = 0.15)
add("PY_BINBETA_B1_prob_rate", f$probs[["drug"]])
add("PY_BINBETA_B1_prob_ratio", f$probs[["ratio"]])
add("PY_BINBETA_B1_mean_theta_c", f$estimates[["theta_c"]])
add("PY_BINBETA_B1_mean_theta_d", f$estimates[["theta_d"]])
add("PY_BINBETA_B1_mean_rr", f$estimates[["rr"]])
add_q("PY_BINBETA_B1_q_theta_d", f$rate)
add_q("PY_BINBETA_B1_q_rr", f$rr)
f <- binbeta_post(12, 80, k = 2, prior_c = c(5, 45), prior_d = c(0.5, 0.5),
                  rr0 = 1.5)
add("PY_BINBETA_B2_prob_ratio", f$prob)
add("PY_BINBETA_B2_mean_theta_c", f$estimates[["theta_c"]])
add("PY_BINBETA_B2_mean_theta_d", f$estimates[["theta_d"]])
add("PY_BINBETA_B2_mean_rr", f$estimates[["rr"]])
add_q("PY_BINBETA_B2_q_theta_d", f$rate)
add_q("PY_BINBETA_B2_q_rr", f$rr)
# B3 to B5: informative priors that conflict with the data (the posterior
# mean of the ratio is infinite in B3 and B4, whose control prior has
# shape1 = 1)
cases <- list(
  B3 = list(y = 12, n = 100, k = 1, prior_c = c(1, 49), prior_d = c(1, 49),
            p0 = 0.05, rr0 = 1),
  B4 = list(y = 100, n = 100, k = 1, prior_c = c(1, 99), prior_d = c(1, 99),
            p0 = 0.3, rr0 = 2),
  B5 = list(y = 5, n = 1000, k = 1, prior_c = c(5, 495),
            prior_d = c(30, 270), p0 = 0.04, rr0 = 8))
for (id in names(cases)) {
  f <- do.call(binbeta_post, cases[[id]])
  add(sprintf("PY_BINBETA_%s_prob_rate", id), f$probs[["drug"]])
  add(sprintf("PY_BINBETA_%s_prob_ratio", id), f$probs[["ratio"]])
  add(sprintf("PY_BINBETA_%s_mean_theta_c", id), f$estimates[["theta_c"]])
  add(sprintf("PY_BINBETA_%s_mean_theta_d", id), f$estimates[["theta_d"]])
  if (id == "B5") {
    add("PY_BINBETA_B5_mean_rr", f$estimates[["rr"]])
  }
  add_q(sprintf("PY_BINBETA_%s_q_theta_d", id), f$rate)
  add_q(sprintf("PY_BINBETA_%s_q_rr", id), f$rr)
}

# CANARY ------------------------------------------------------------------
cases <- list(
  C1 = list(y = 12, n = 100, k = 1, prior_c = c(1, 49), metric = "ratio",
            crit = 1),
  C2 = list(y = 12, n = 100, k = 1, prior_c = c(1, 49),
            metric = "difference", crit = 0.02),
  C3 = list(y = 20, n = 120, k = 2, prior_c = c(5, 45), metric = "odds",
            crit = 1.2),
  C4 = list(y = 250, n = 20000, k = 1, prior_c = c(75.71, 8500),
            metric = "ratio", crit = 2, likelihood = "poisson"),
  C5 = list(y = 100, n = 100, k = 1, prior_c = c(1, 49), metric = "ratio",
            crit = 2))
for (id in names(cases)) {
  f <- do.call(canary_post, cases[[id]])
  add(sprintf("PY_CANARY_%s_prob", id), f$prob)
  add(sprintf("PY_CANARY_%s_mean_metric", id), f$rr[["mean"]])
  if (id == "C4") {
    add("PY_CANARY_C4_mean_lambda_c", f$estimates[["lambda_c"]])
  } else {
    add(sprintf("PY_CANARY_%s_mean_p_c", id), f$estimates[["p_c"]])
  }
  # the reference has no quantiles for C5 (ref_canary.py)
  if (id != "C5") {
    add_q(sprintf("PY_CANARY_%s_q", id), f$rr)
  }
}

# MAPBERN and MAP prior ---------------------------------------------------
d <- get_data("M1")
mix <- data.frame(w = c(0.7, 0.3), shape1 = c(20, 2), shape2 = c(180, 18))
for (cs in list(list(id = "M1", rr0 = 1, prior_d = c(1 / 3, 1 / 3)),
                list(id = "M2", rr0 = 1.5, prior_d = c(1 / 3, 1 / 3)),
                list(id = "M3", rr0 = 0.25, prior_d = c(2, 98)))) {
  f <- mapbern_post(d$y, d$t, k = 3, prior_c = mix, prior_d = cs$prior_d,
                    rr0 = cs$rr0)
  add(sprintf("PY_MAPBERN_%s_prob_ratio", cs$id), f$prob)
  add(sprintf("PY_MAPBERN_%s_mean_pi_c", cs$id), f$estimates[["pi_c"]])
  add(sprintf("PY_MAPBERN_%s_mean_pi_d", cs$id), f$estimates[["pi_d"]])
  add_q(sprintf("PY_MAPBERN_%s_q_pi_d", cs$id), f$rate)
  add_q(sprintf("PY_MAPBERN_%s_q_rr", cs$id), f$rr)
}
for (cs in list(list(id = "P1", e = c(9, 11, 10), n = c(100, 100, 100)),
                list(id = "P2", e = c(3, 9, 1, 6), n = c(40, 60, 30, 80)))) {
  p <- map_beta(events = cs$e, n = cs$n)
  add(sprintf("PY_MAP_%s_mean", cs$id), attr(p, "moments")[["mean"]])
  add(sprintf("PY_MAP_%s_sd", cs$id), attr(p, "moments")[["sd"]])
  for (x in c(0.05, 0.08, 0.1, 0.12, 0.15)) {
    add(sprintf("PY_MAP_%s_cdf_%g", cs$id, x), pmixbeta(x, p))
  }
}
# MS1 to MS7: MAP priors from map_beta() for historical control arms with
# few events (python/ref_mapbern_sparse.py, which integrates the MAP model
# over (a, b, pi_c, pi_d) without the mixture of map_beta()); complete
# follow-up, k = 3 and prior_d = Beta(1/3, 1/3). MS2 and MS4 (no
# historical events) use n_u = 96: with the default 48 nodes the errors are
# about 8e-6 for the probabilities of MS2 and 2e-7 for the posterior mean
# of pi_c of MS4 (relative error 6e-5), below the error estimates of about
# 1e-4 and 9e-6 that mapbern_post() reports.
ms <- list(MS1 = list(e = c(1, 1, 0), n = c(40, 40, 40), y = 1, m = 50),
           MS2 = list(e = c(0, 0, 0), n = c(40, 40, 40), y = 1, m = 50),
           MS3 = list(e = c(2, 1, 1), n = c(40, 40, 40), y = 2, m = 100),
           MS4 = list(e = c(0, 0, 0), n = c(40, 40, 40), y = 10, m = 50),
           MS5 = list(e = c(2, 5, 7, 5, 10, 14, 10, 19, 28),
                      n = rep(c(50, 100, 200), each = 3), y = 95, m = 1000),
           MS6 = list(e = c(1, 2, 3), n = c(100, 100, 100), y = 3, m = 200),
           MS7 = list(e = c(9, 11, 10), n = c(100, 100, 100), y = 30,
                      m = 200))
for (id in names(ms)) {
  cs <- ms[[id]]
  p <- map_beta(events = cs$e, n = cs$n)
  yy <- rep(c(1, 0), c(cs$y, cs$m - cs$y))
  ctrl <- if (id %in% c("MS2", "MS4")) list(n_u = 96L) else list()
  f1 <- mapbern_post(yy, rep(1, cs$m), k = 3, prior_c = p, rate0 = 0.05,
                     details = FALSE, control = ctrl)
  f2 <- mapbern_post(yy, rep(1, cs$m), k = 3, prior_c = p, rr0 = 2,
                     details = FALSE, control = ctrl)
  pre <- sprintf("PY_MAPBERN_%s_", id)
  add(paste0(pre, "prob_ratio_1"), f1$probs[["ratio"]])
  add(paste0(pre, "prob_ratio_2"), f2$prob)
  add(paste0(pre, "prob_rate"), f1$probs[["drug"]])
  add(paste0(pre, "mean_pc"), f1$estimates[["pi_c"]])
  add(paste0(pre, "mean_pd"), f1$estimates[["pi_d"]])
}
d <- get_data("MAP_SUBJ")
p <- map_beta(y = d$y, t = d$t, study = d$study)
# the reference for subject-level data uses a coarser grid (relative error
# of about 1e-6), so a relative tolerance of 1e-5 is used for both moments
add("PY_MAP_SUBJ_mean", attr(p, "moments")[["mean"]])
add("PY_MAP_SUBJ_sd", attr(p, "moments")[["sd"]])
rows[[length(rows) - 1L]]$tolerance <- 1e-5 * rows[[length(rows) - 1L]]$published
rows[[length(rows)]]$tolerance <- 1e-5 * rows[[length(rows)]]$published

# BS1 and ES1: control priors with a shape parameter below 1, integrated
# over the normal quantile of the prior (python/ref_shape_below1.py: the
# exact finite mixture of Beta products for BS1 and a dense rule on the
# normal quantiles of both priors for ES1)
f <- binbeta_post(12, 100, k = 1, prior_c = c(0.5, 5), prior_d = c(1, 1),
                  p0 = 0.15, details = FALSE)
add("PY_BINBETA_BS1_prob_rate", f$probs[["drug"]])
add("PY_BINBETA_BS1_prob_ratio", f$probs[["ratio"]])
add("PY_BINBETA_BS1_mean_theta_c", f$estimates[["theta_c"]])
add("PY_BINBETA_BS1_mean_theta_d", f$estimates[["theta_d"]])
d <- get_data("E1")
f <- exptime_post(d$y, d$t, k = 2, prior_c = c(0.7, 30), prior_d = c(2, 70),
                  rate0 = 0.04, details = FALSE)
add("PY_EXPTIME_ES1_prob_rate", f$probs[["drug"]])
add("PY_EXPTIME_ES1_prob_ratio", f$probs[["ratio"]])
add("PY_EXPTIME_ES1_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_ES1_mean_lambda_d", f$estimates[["lambda_d"]])

# MULTIBIN ----------------------------------------------------------------
f <- multibin_post(c(6, 10, 25, 47), c(200, 400, 800, 1600), k = 1,
                   prior_c = c(2, 98))
add("PY_MULTIBIN_G1_prob", f$prob)
add("PY_MULTIBIN_G1_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G1_q", f$rr)
f <- multibin_post(c(3, 9), c(150, 300), k = c(1, 2),
                   prior_c = rbind(c(2, 98), c(4, 196)), lrr_prior = c(0, 2),
                   rr0 = 1.5)
add("PY_MULTIBIN_G2_prob", f$prob)
add("PY_MULTIBIN_G2_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G2_q", f$rr)
# G3 and G4: an informative prior of log(delta) that conflicts with the data
f <- multibin_post(c(60, 120), c(1000, 2000), k = 1,
                   prior_c = rbind(c(20, 1980), c(40, 3960)),
                   lrr_prior = c(0, 0.1), rr0 = 1.5)
add("PY_MULTIBIN_G3_prob", f$prob)
add("PY_MULTIBIN_G3_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G3_q", f$rr)
f <- multibin_post(c(2, 3), c(1000, 2000), k = 1,
                   prior_c = rbind(c(200, 19800), c(400, 39600)),
                   lrr_prior = c(0, 0.1), rr0 = 0.8)
add("PY_MULTIBIN_G4_prob", f$prob)
add("PY_MULTIBIN_G4_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G4_q", f$rr)

# EXP_BERN ----------------------------------------------------------------
# X1 to X3: informative priors of p that conflict with the data (the event
# rate of data set X1 is about 2.7 times the historical placebo rate, that
# of data set X3 about half of it)
for (cs in list(list(id = "X1", data = "X1", k = 1, p_prior = c(10, 10),
                     rr0 = 1),
                list(id = "X2", data = "X1", k = 1, p_prior = c(30, 30),
                     rr0 = 1.5),
                list(id = "X3", data = "X3", k = 2, p_prior = c(10, 10),
                     rr0 = 1))) {
  d <- get_data(cs$data)
  f <- expbern_post(d$y, d$t, k = cs$k, lambda_hp = -log(0.9) / 52,
                    pt_hp = nrow(d) / 2 * 52, df = 2, p_prior = cs$p_prior,
                    rr0 = cs$rr0, t_star = 52)
  pre <- sprintf("PY_EXPBERN_%s_", cs$id)
  add(paste0(pre, "prob"), f$prob)
  add(paste0(pre, "mean_p"), f$estimates[["p"]])
  add(paste0(pre, "mean_r"), f$estimates[["r"]])
  add(paste0(pre, "mean_lambda_p"), f$estimates[["lambda_p"]])
  add(paste0(pre, "mean_lambda_d"), f$estimates[["lambda_d"]])
  add(paste0(pre, "ip_placebo"), f$incidence[["placebo"]])
  add(paste0(pre, "ip_drug"), f$incidence[["drug"]])
  add_q(paste0(pre, "q"), f$rr)
}

missing <- setdiff(ref$id, vapply(rows, function(r) r$id, character(1)))
if (length(missing) > 0L) {
  stop("Reference values not compared: ", paste(missing, collapse = ", "))
}
tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_methods_python")
out <- tab
names(out)[names(out) == "published"] <- "reference"
names(out)[names(out) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "inst/validation/results/validate_methods_python.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
