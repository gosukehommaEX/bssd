# Validation of exptime_post(), binbeta_post(), mapbern_post(),
# multibin_post() and canary_post() (Poisson likelihood) against MCMC with
# JAGS. The JAGS models below are written directly from the model
# specifications of Schnell and Ball (2016), Ball and Schnell (2016), Brock
# et al. (2023), Goren et al. (2020) and Gould and Wang (2017), with the
# treatment assignments of the subjects as latent variables where the papers
# define them, and share no code with bssd.
#
# Requirements: JAGS 4.x and the R packages rjags and coda.
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_jags_methods.R")
# Output: inst/validation/results/validate_jags_methods.csv
# Run time: about five to ten minutes.
#
# An item is PASS when the bssd value is within four Monte Carlo standard
# errors of the MCMC value (standard errors from coda::effectiveSize).

for (pkg in c("rjags", "coda")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed (rjags also needs JAGS).")
  }
}
source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")
dat <- utils::read.csv("inst/validation/python/data_methods.csv",
                       stringsAsFactors = FALSE)

n_burn <- 5000
n_iter <- 50000
n_chain <- 4

run_jags <- function(model, data, inits, vars) {
  jm <- rjags::jags.model(textConnection(model), data = data, inits = inits,
                          n.chains = n_chain, quiet = TRUE)
  stats::update(jm, n_burn, progress.bar = "none")
  smp <- rjags::coda.samples(jm, vars, n.iter = n_iter,
                             progress.bar = "none")
  as.data.frame(do.call(rbind, lapply(smp, as.matrix)))
}
mc_prop <- function(ind) {
  ess <- unname(coda::effectiveSize(ind))
  p <- mean(ind)
  c(value = p, se = sqrt(p * (1 - p) / ess), ess = ess)
}
mc_mean <- function(x) {
  ess <- unname(coda::effectiveSize(x))
  c(value = mean(x), se = stats::sd(x) / sqrt(ess), ess = ess)
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
  add(id, model, quantity, mc[["value"]], bssd, 4 * mc[["se"]],
      sprintf("MCMC SE %.4f, ESS %.0f", mc[["se"]], mc[["ess"]]))
}
add_median <- function(id, model, x, med) {
  add_mc(id, model, "MCMC proportion below the bssd posterior median",
         c(value = 0.5, mc_prop(as.numeric(x <= med))[-1]), mean(x <= med))
}
inits_of <- function(seed, ...) {
  lapply(seq_len(n_chain), function(j) {
    c(list(...), list(.RNG.name = "base::Mersenne-Twister",
                      .RNG.seed = seed * 10 + j))
  })
}

# EXPTIME: Schnell and Ball (2016), latent assignments -----------------------
model_exptime <- "
model {
  for (i in 1:N) {
    a[i] ~ dbern(wd)
    y[i] ~ dpois(((1 - a[i]) * lc + a[i] * ld) * t[i])
  }
  lc ~ dgamma(ac, bc)
  ld ~ dgamma(ad, bd)
  r <- ld / lc
}"
d <- dat[dat$dataset == "E1", ]
fit <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                    prior_d = c(2, 70), rate0 = 0.03, rr0 = 1.2)
dr <- run_jags(model_exptime,
               list(N = nrow(d), y = d$y, t = d$t, wd = 2 / 3, ac = 30,
                    bc = 1500, ad = 2, bd = 70),
               inits_of(1, lc = 0.02, ld = 0.03), c("lc", "ld", "r"))
add_mc("J2_EXPTIME_RATE", "EXPTIME", "P(lambda_d > 0.03 | data)",
       mc_prop(as.numeric(dr$ld > 0.03)), fit$probs[["drug"]])
add_mc("J2_EXPTIME_RATIO", "EXPTIME", "P(r > 1.2 | data)",
       mc_prop(as.numeric(dr$r > 1.2)), fit$probs[["ratio"]])
add_mc("J2_EXPTIME_LC", "EXPTIME", "posterior mean of lambda_c",
       mc_mean(dr$lc), fit$estimates[["lambda_c"]])
add_median("J2_EXPTIME_MED", "EXPTIME", dr$r, fit$rr[["median"]])

# BINBETA: Ball and Schnell (2016), latent number on control ---------------
model_binbeta <- "
model {
  m0 ~ dbin(wc, N)
  y0 ~ dbin(tc, m0)
  y1 ~ dbin(td, N - m0)
  z ~ dsum(y0, y1)
  tc ~ dbeta(a, b)
  td ~ dbeta(c, d)
  rr <- td / tc
}"
fit <- binbeta_post(40, 300, k = 1, prior_c = c(31, 271), prior_d = c(1, 1),
                    p0 = 0.15)
# the marginal model Z ~ Binomial(N, w_c theta_c + w_d theta_d) is used if
# JAGS cannot sample the latent counts
model_binbeta_marginal <- "
model {
  z ~ dbin(wc * tc + (1 - wc) * td, N)
  tc ~ dbeta(a, b)
  td ~ dbeta(c, d)
  rr <- td / tc
}"
bb_data <- list(N = 300, z = 40, wc = 0.5, a = 31, b = 271, c = 1, d = 1)
bb_note <- "latent counts"
dr <- tryCatch(
  run_jags(model_binbeta, bb_data,
           inits_of(2, tc = 0.1, td = 0.15, m0 = 150, y0 = 15, y1 = 25),
           c("tc", "td", "rr")),
  error = function(e) {
    bb_note <<- paste("marginal model; the latent model failed:",
                      conditionMessage(e))
    run_jags(model_binbeta_marginal, bb_data,
             inits_of(2, tc = 0.1, td = 0.15), c("tc", "td", "rr"))
  })
add_mc("J2_BINBETA_RATE", "BINBETA", "P(theta_d > 0.15 | data)",
       mc_prop(as.numeric(dr$td > 0.15)), fit$probs[["drug"]])
add_mc("J2_BINBETA_RATIO", "BINBETA", "P(rr > 1 | data)",
       mc_prop(as.numeric(dr$rr > 1)), fit$probs[["ratio"]])
add_median("J2_BINBETA_MED", "BINBETA", dr$rr, fit$rr[["median"]])
add("J2_BINBETA_MODEL", "BINBETA", "JAGS model used", NA, NA, NA, bb_note)

# MAPBERN: Brock et al. (2023), mixture prior and latent assignments -------
model_mapbern <- "
model {
  comp ~ dcat(w[])
  pc ~ dbeta(s1[comp], s2[comp])
  pd ~ dbeta(ad, bd)
  for (i in 1:N) {
    a[i] ~ dbern(wd)
    y[i] ~ dbern((1 - a[i]) * (1 - pow(1 - pc, t[i])) +
                   a[i] * (1 - pow(1 - pd, t[i])))
  }
  rr <- pd / pc
}"
d <- dat[dat$dataset == "M1", ]
mix <- data.frame(w = c(0.7, 0.3), shape1 = c(20, 2), shape2 = c(180, 18))
fit <- mapbern_post(d$y, d$t, k = 3, prior_c = mix)
dr <- run_jags(model_mapbern,
               list(N = nrow(d), y = d$y, t = d$t, wd = 0.75, w = mix$w,
                    s1 = mix$shape1, s2 = mix$shape2, ad = 1 / 3,
                    bd = 1 / 3),
               inits_of(3, pc = 0.1, pd = 0.2), c("pc", "pd", "rr"))
add_mc("J2_MAPBERN_PROB", "MAPBERN", "P(pi_d > pi_c | data)",
       mc_prop(as.numeric(dr$rr > 1)), fit$prob)
add_mc("J2_MAPBERN_PC", "MAPBERN", "posterior mean of pi_c", mc_mean(dr$pc),
       fit$estimates[["pi_c"]])
add_median("J2_MAPBERN_MED", "MAPBERN", dr$rr, fit$rr[["median"]])

# MULTIBIN: Goren et al. (2020), theta_s <= 1 / delta by the ones trick ----
model_multibin <- "
model {
  ld ~ dnorm(m, 1 / (v * v))
  delta <- exp(ld)
  for (s in 1:S) {
    th[s] ~ dbeta(a[s], b[s])
    ones[s] ~ dbern(step(1 - delta * th[s]))
    y[s] ~ dbin(min(th[s] * (1 - r[s] + r[s] * delta), 1), n[s])
  }
}"
y <- c(6, 10, 25, 47)
n <- c(200, 400, 800, 1600)
fit <- multibin_post(y, n, k = 1, prior_c = c(2, 98))
dr <- run_jags(model_multibin,
               list(S = 4, y = y, n = n, r = rep(0.5, 4), a = rep(2, 4),
                    b = rep(98, 4), m = 0, v = 5, ones = rep(1, 4)),
               inits_of(4, ld = 0, th = rep(0.02, 4)), c("delta", "ld"))
add_mc("J2_MULTIBIN_PROB", "MULTIBIN", "P(delta > 1 | data)",
       mc_prop(as.numeric(dr$delta > 1)), fit$prob)
add_mc("J2_MULTIBIN_LD", "MULTIBIN", "posterior mean of log(delta)",
       mc_mean(dr$ld), fit$estimates[["log_delta"]])
add_median("J2_MULTIBIN_MED", "MULTIBIN", dr$delta, fit$rr[["median"]])

# CANARY (Poisson): Gould and Wang (2017) -----------------------------------
model_canary <- "
model {
  y ~ dpois(E * lc * (1 + k * R) / (k + 1))
  lc ~ dgamma(a, b)
  lr ~ dnorm(mu, 1 / (s * s))
  R <- exp(lr)
}"
fit <- canary_post(250, 20000, k = 1, prior_c = c(75.71, 8500),
                   likelihood = "poisson", crit = 2)
mp <- fit$prior$metric_prior
dr <- run_jags(model_canary,
               list(y = 250, E = 20000, k = 1, a = 75.71, b = 8500,
                    mu = mp[["meanlog"]], s = mp[["sdlog"]]),
               inits_of(5, lc = 0.009, lr = 0.5), c("R", "lc"))
add_mc("J2_CANARY_PROB", "CANARY", "P(R > 2 | data), Poisson",
       mc_prop(as.numeric(dr$R > 2)), fit$prob)
add_mc("J2_CANARY_LC", "CANARY", "posterior mean of lambda_c, Poisson",
       mc_mean(dr$lc), fit$estimates[["lambda_c"]])
add_median("J2_CANARY_MED", "CANARY", dr$R, fit$rr[["median"]])

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_jags_methods")
out <- tab
names(out)[names(out) == "published"] <- "reference"
names(out)[names(out) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "inst/validation/results/validate_jags_methods.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
