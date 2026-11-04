# Validation of expbern_post() and llbern_post() against MCMC with JAGS.
# The JAGS models below are written directly from the model specifications
# of Duncan (2026, Sections 3.3 and 3.4) and share no code with bssd.
#
# Requirements: JAGS 4.x (https://mcmc-jags.sourceforge.io/) and the R
# packages rjags and coda. The rjags 4-x package on CRAN does not work with
# JAGS 5.
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_jags.R")
# Output: inst/validation/results/validate_jags.csv
# Run time: about 10 to 20 minutes.
#
# The data sets are simulated with sim_blinded() with fixed seeds. For each
# data set, the posterior probability of a signal, the MCMC proportion of
# draws below the bssd posterior median of the relative risk, and the
# posterior means of the cumulative incidences at t_star are compared. An
# item is PASS when the bssd value is within four Monte Carlo standard errors
# of the MCMC value (standard errors from coda::effectiveSize).

for (pkg in c("rjags", "coda")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed (rjags also needs JAGS).")
  }
}
source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")

n_burn <- 5000
n_iter <- 50000
n_chain <- 4

model_exp <- "
model {
  for (i in 1:N) {
    pp[i] <- 1 - exp(-lambda * t[i])
    pd[i] <- 1 - exp(-r * lambda * t[i])
    y[i] ~ dbern(wp * pp[i] + (1 - wp) * pd[i])
  }
  lambda ~ dgamma(alpha, beta)
  p ~ dbeta(a, b)
  r <- p / (k * (1 - p))
}"
model_ll <- "
model {
  for (i in 1:N) {
    pp[i] <- ilogit((log(t[i]) - mu_p) / sigma_p)
    pd[i] <- ilogit((log(t[i]) - mu_d) / sigma_d)
    y[i] ~ dbern(wp * pp[i] + (1 - wp) * pd[i])
  }
  mu_p ~ dgamma(s1, r1)
  sigma_p ~ dgamma(s2, r2)
  mu_d ~ dgamma(s3, r3)
  sigma_d ~ dgamma(s4, r4)
}"

run_jags <- function(model, data, inits, vars) {
  jm <- rjags::jags.model(textConnection(model), data = data, inits = inits,
                          n.chains = n_chain, quiet = TRUE)
  stats::update(jm, n_burn, progress.bar = "none")
  smp <- rjags::coda.samples(jm, vars, n.iter = n_iter,
                             progress.bar = "none")
  as.data.frame(do.call(rbind, lapply(smp, as.matrix)))
}
mc_prop <- function(ind) {
  ess <- coda::effectiveSize(ind)
  p <- mean(ind)
  c(value = p, se = sqrt(p * (1 - p) / ess), ess = unname(ess))
}
mc_mean <- function(x) {
  ess <- coda::effectiveSize(x)
  c(value = mean(x), se = stats::sd(x) / sqrt(ess), ess = unname(ess))
}

rows <- list()
add <- function(id, source, quantity, reference, bssd, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = reference,
    reproduced = bssd, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}
add_mc <- function(id, model, quantity, mc, bssd) {
  add(id, model, quantity, unname(mc["value"]), bssd, 4 * unname(mc["se"]),
      sprintf("MCMC SE %.4f, ESS %.0f", mc["se"], mc["ess"]))
}

# EXP_BERN data sets: k, design and prior settings
exp_cases <- list(
  list(id = "E1", seed = 101, k = 1, n = 150, ip_p = 0.10, ip_d = 0.20,
       tau = 52, ip_prior = 0.10, a = 0.5, b = 0.5, rr0 = 1),
  list(id = "E2", seed = 102, k = 1, n = 300, ip_p = 0.30, ip_d = 0.30,
       tau = 78, ip_prior = 0.30, a = 0.5, b = 0.5, rr0 = 1),
  list(id = "E3", seed = 103, k = 2, n = 200, ip_p = 0.05, ip_d = 0.10,
       tau = 40, ip_prior = 0.05, a = 1, b = 1, rr0 = 1.5)
)
for (cs in exp_cases) {
  set.seed(cs$seed)
  d <- sim_blinded(n = cs$n, k = cs$k, accrual = 26, follow_up = 52,
                   analysis_time = cs$tau, par_p = -log(1 - cs$ip_p) / 52,
                   par_d = -log(1 - cs$ip_d) / 52)
  lam_hp <- -log(1 - cs$ip_prior) / 52
  pt_hp <- cs$n / 2 * 52
  fit <- expbern_post(d$event, d$follow_up, k = cs$k, lambda_hp = lam_hp,
                      pt_hp = pt_hp, df = 2, p_prior = c(cs$a, cs$b),
                      rr0 = cs$rr0, t_star = 52)
  dr <- run_jags(model_exp,
                 list(N = nrow(d), t = d$follow_up, y = d$event,
                      wp = 1 / (cs$k + 1), k = cs$k,
                      alpha = lam_hp * pt_hp / 2, beta = pt_hp / 2,
                      a = cs$a, b = cs$b),
                 lapply(seq_len(n_chain), function(j) {
                   list(lambda = lam_hp * c(0.8, 1, 1.2, 1.1)[j],
                        p = c(0.3, 0.5, 0.6, 0.4)[j],
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = cs$seed * 10 + j)
                 }),
                 c("lambda", "r"))
  lab <- sprintf("%s: n = %d, events = %d", cs$id, nrow(d), sum(d$event))
  add_mc(sprintf("J_%s_PROB", cs$id), "EXP_BERN",
         sprintf("P(r > %g | data), %s", cs$rr0, lab),
         mc_prop(as.numeric(dr$r > cs$rr0)), fit$prob)
  add_mc(sprintf("J_%s_MED", cs$id), "EXP_BERN",
         sprintf("MCMC proportion of r below the bssd median, %s", lab),
         c(value = 0.5, mc_prop(as.numeric(dr$r <= fit$rr[["median"]]))[-1]),
         mean(dr$r <= fit$rr[["median"]]))
  add_mc(sprintf("J_%s_IPP", cs$id), "EXP_BERN",
         sprintf("posterior mean of the placebo incidence at week 52, %s",
                 lab), mc_mean(1 - exp(-dr$lambda * 52)),
         fit$incidence[["placebo"]])
  add_mc(sprintf("J_%s_IPD", cs$id), "EXP_BERN",
         sprintf("posterior mean of the study drug incidence at week 52, %s",
                 lab), mc_mean(1 - exp(-dr$r * dr$lambda * 52)),
         fit$incidence[["drug"]])
}

# LL_BERN data sets
ll_cases <- list(
  list(id = "L1", seed = 201, k = 1, n = 150, par_p = c(7.3, 0.8),
       par_d = c(6.1, 0.8), follow_up = 104, tau = 117, prior = c(7.3, 0.8),
       rr0 = 1),
  list(id = "L2", seed = 202, k = 1, n = 200, par_p = c(6.17, 1.6),
       par_d = c(6.17, 1.6), follow_up = 52, tau = 104 / 3,
       prior = c(6.17, 1.6), rr0 = 1),
  list(id = "L3", seed = 203, k = 1, n = 100, par_p = c(7.3, 0.8),
       par_d = c(7.3, 0.8), follow_up = 104, tau = 78, prior = c(7.3, 0.8),
       rr0 = 1)
)
for (cs in ll_cases) {
  set.seed(cs$seed)
  d <- sim_blinded(n = cs$n, k = cs$k, accrual = 26,
                   follow_up = cs$follow_up, analysis_time = cs$tau,
                   dist = "loglogistic", par_p = cs$par_p, par_d = cs$par_d)
  pt_hp <- cs$n / 2 * cs$follow_up / 52
  fit <- llbern_post(d$event, d$follow_up, k = cs$k, mu_hp = cs$prior[1],
                     sigma_hp = cs$prior[2], pt_hp = pt_hp, df_p = 2,
                     t_star = cs$follow_up, rr0 = cs$rr0)
  shape <- c(pt_hp * cs$prior[1] / 2, pt_hp * cs$prior[2] / 2, cs$prior[1],
             cs$prior[2])
  rate <- c(pt_hp / 2, pt_hp / 2, 1, 1)
  dr <- run_jags(model_ll,
                 list(N = nrow(d), t = d$follow_up, y = d$event,
                      wp = 1 / (cs$k + 1), s1 = shape[1], r1 = rate[1],
                      s2 = shape[2], r2 = rate[2], s3 = shape[3],
                      r3 = rate[3], s4 = shape[4], r4 = rate[4]),
                 lapply(seq_len(n_chain), function(j) {
                   list(mu_p = cs$prior[1], sigma_p = cs$prior[2],
                        mu_d = cs$prior[1] * c(0.8, 0.9, 1, 1.1)[j],
                        sigma_d = cs$prior[2] * c(0.8, 1, 1.2, 0.9)[j],
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = cs$seed * 10 + j)
                 }),
                 c("mu_p", "sigma_p", "mu_d", "sigma_d"))
  ts <- cs$follow_up
  fp <- stats::plogis((log(ts) - dr$mu_p) / dr$sigma_p)
  fd <- stats::plogis((log(ts) - dr$mu_d) / dr$sigma_d)
  rr <- fd / fp
  lab <- sprintf("%s: n = %d, events = %d", cs$id, nrow(d), sum(d$event))
  add_mc(sprintf("J_%s_PROB", cs$id), "LL_BERN",
         sprintf("P(RR(t_star) > %g | data), %s", cs$rr0, lab),
         mc_prop(as.numeric(rr > cs$rr0)), fit$prob)
  add_mc(sprintf("J_%s_MED", cs$id), "LL_BERN",
         sprintf("MCMC proportion of RR below the bssd median, %s", lab),
         c(value = 0.5, mc_prop(as.numeric(rr <= fit$rr[["median"]]))[-1]),
         mean(rr <= fit$rr[["median"]]))
  add_mc(sprintf("J_%s_IPP", cs$id), "LL_BERN",
         sprintf("posterior mean of the placebo incidence at t_star, %s",
                 lab), mc_mean(fp), fit$incidence[["placebo"]])
  add_mc(sprintf("J_%s_IPD", cs$id), "LL_BERN",
         sprintf("posterior mean of the study drug incidence at t_star, %s",
                 lab), mc_mean(fd), fit$incidence[["drug"]])
  add(sprintf("J_%s_ERR", cs$id), "LL_BERN",
      sprintf("integration error estimate of bssd, %s", lab), NA,
      fit$integration$error, NA)
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_jags")
self_check(tab)
names(tab)[names(tab) == "published"] <- "reference"
names(tab)[names(tab) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/validation/results/validate_jags.csv",
                 row.names = FALSE)
