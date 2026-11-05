# LL_BERN: log-logistic Bernoulli mixture model

``` r

library(bssd)
```

## Model

The log-logistic Bernoulli (LL_BERN) model of Duncan (2026, Section 3.4)
replaces the exponential distribution of EXP_BERN
([`vignette("expbern")`](https://gosukehommaex.github.io/bssd/articles/expbern.md))
by the log-logistic distribution, so that the hazard of the event can
change over time. For subject \\i\\ with the follow-up time \\t_i\\, the
event indicator \\Y_i\\ follows \\P(Y_i = 1) = w_p F(t_i; \mu_p,
\sigma_p) + w_d F(t_i; \mu_d, \sigma_d), \quad F(t; \mu, \sigma) = 1 -
\frac{1}{1 + \exp\\(\log t - \mu) / \sigma\\},\\ where \\w_p = 1 / (k +
1)\\ and \\w_d = k / (k + 1)\\ for the randomization ratio \\k:1\\, and
\\(\mu_p, \sigma_p)\\ and \\(\mu_d, \sigma_d)\\ are the location and
scale parameters on placebo and on the study drug. The hazard is
decreasing when \\\sigma \ge 1\\ and first increases and then decreases
when \\\sigma \< 1\\. With the values \\\mu\\ and \\\sigma\\ estimated
from historical placebo data, the priors are \\\mu_p \sim
\mathrm{Gamma}(PT\_{HP} \mu / DF_p, PT\_{HP} / DF_p), \quad \sigma_p
\sim \mathrm{Gamma}(PT\_{HP} \sigma / DF_p, PT\_{HP} / DF_p),\\ \\\mu_d
\sim \mathrm{Gamma}(\mu, 1), \quad \sigma_d \sim \mathrm{Gamma}(\sigma,
1),\\ with the shape and rate parameterization. The priors of the
placebo parameters have the means \\\mu\\ and \\\sigma\\, and \\PT\_{HP}
/ DF_p\\ sets their precision; the priors of the study drug parameters
are weakly informative.

Duncan (2026) does not state which relative risk is used for the signal
criterion of this model. bssd uses the ratio of the cumulative
incidences at a time \\t^\*\\, \\RR(t^\*) = F(t^\*; \mu_d, \sigma_d) /
F(t^\*; \mu_p, \sigma_p),\\ and declares a signal when \\P\\RR(t^\*) \>
\theta_c \mid \mathrm{data}\\ \> P\_{\mathrm{cutoff}}\\.

The arguments of
[`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
correspond to the symbols as follows.

| Symbol | Argument | Description |
|----|----|----|
| \\Y_i\\ | `y` | event indicator of each subject (1 if the subject had the event during follow-up) |
| \\t_i\\ | `t` | follow-up time of each subject |
| \\k\\ | `k` | randomization ratio of study drug to placebo |
| \\\mu\\, \\\sigma\\ | `mu_hp`, `sigma_hp` | log-logistic parameters estimated from the historical placebo data |
| \\PT\_{HP}\\ | `pt_hp` | total patient time of the historical placebo data |
| \\DF_p\\ | `df_p` | discount factor |
|  | `mu_d_prior`, `sigma_d_prior` | shape and rate of the priors of \\\mu_d\\ and \\\sigma_d\\ (defaults c(mu_hp, 1) and c(sigma_hp, 1)) |
| \\t^\*\\ | `t_star` | time at which the relative risk and the cumulative incidences are evaluated |
| \\\theta_c\\ | `rr0` | critical value of \\RR(t^\*)\\ (default 1) |
| \\P\_{\mathrm{cutoff}}\\ | `cutoff` | posterior probability threshold (default 0.95) |

The parameters \\\mu\\ and \\\sigma\\ refer to the time unit of `t` and
`t_star`.
[`ll_cdf()`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)
computes \\F(t; \mu, \sigma)\\, and
[`ll_mu()`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md)
gives the location parameter for a cumulative incidence at a given time.

## Example

The following data are simulated with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md):
150 subjects randomized 1:1 with uniform entry over 26 weeks, up to 104
weeks of follow-up, and log-logistic times to the event with \\\sigma =
0.8\\. The cumulative incidence at week 104 is 3.5% on placebo (\\\mu_p
= 7.3\\) and 12% on the study drug, whose location parameter is obtained
with
[`ll_mu()`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md).
The blinded analysis is at week 117.

``` r

set.seed(2026)
mu_d <- ll_mu(0.12, t = 104, sigma = 0.8)
mu_d
#> [1] 6.238335
d <- sim_blinded(n = 150, k = 1, accrual = 26, follow_up = 104,
                 analysis_time = 117, dist = "loglogistic",
                 par_p = c(7.3, 0.8), par_d = c(mu_d, 0.8))
c(subjects = nrow(d), events = sum(d$event))
#> subjects   events 
#>      150       10
```

The priors use the placebo parameters as `mu_hp` and `sigma_hp`, and a
historical placebo group of 75 subjects followed for 104 weeks (150
patient-years) discounted by \\DF_p = 2\\, as in the simulations in the
folder `reproduce` of the installed package. The relative risk is
evaluated at the end of follow-up, \\t^\* = 104\\ weeks.

``` r

fit <- llbern_post(y = d$event, t = d$follow_up, k = 1, mu_hp = 7.3,
                   sigma_hp = 0.8, pt_hp = 75 * 104 / 52, df_p = 2,
                   t_star = 104)
fit
#> Model: LL_BERN 
#> Relative risk: ratio of cumulative incidences at t_star = 104 
#> P(RR(t_star) > 1 | data) = 0.6967
#> Signal (posterior probability > 0.95): no
#> Posterior of RR(t_star): mean 3.205, median 2.097, 95% interval (5.798e-28, 13.52)
#> Cumulative incidence at t_star = 104: placebo 0.04739, study drug 0.08356, difference 0.03616
summary(fit)
#>     model      prob signal  rr_mean rr_median     rr_lower rr_upper ip_placebo
#> 1 LL_BERN 0.6966964  FALSE 3.204699  2.096947 5.798076e-28 13.52066  0.0473919
#>      ip_drug ip_difference       error
#> 1 0.08355556    0.03616366 0.001012859
```

The posterior probabilities for other critical values of \\RR(104)\\ are
obtained with `rr0`; `details = FALSE` skips the quantiles and
`control = list(check = FALSE)` skips the repeated computation for the
error estimate:

``` r

vapply(c(1, 1.5, 2), function(q) {
  llbern_post(y = d$event, t = d$follow_up, k = 1, mu_hp = 7.3,
              sigma_hp = 0.8, pt_hp = 75 * 104 / 52, df_p = 2,
              t_star = 104, rr0 = q, details = FALSE,
              control = list(check = FALSE))$prob
}, numeric(1))
#> [1] 0.6972028 0.6041112 0.5161246
```

The posterior distribution of \\RR(t^\*)\\ can have a heavy right tail,
so the posterior median is a more stable summary than the posterior
mean.

## Numerical integration

[`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
integrates the four-dimensional posterior numerically on the coordinates
\\\log \mu_p\\, \\\log \sigma_p\\, \\\log \mu_d\\ and the prior
distribution function of \\\sigma_d\\. A coarse grid locates the region
that holds almost all posterior mass. Gauss-Legendre quadrature is used
for three coordinates and Chebyshev interpolation for \\\log \mu_d\\;
because \\RR(t^\*) \> \theta_c\\ holds exactly when \\\mu_d\\ is below a
value determined by the other three parameters, the posterior
probability is obtained from the cumulative integral of the interpolant
(see
[`?llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)).
With `control$check = TRUE`, the default, the computation is repeated on
a finer grid, the finer result is returned, and twice the absolute
difference in `prob` is reported as the error estimate
(`fit$integration$error`, here 1.0e-03). The factor two is explained in
[`?llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md):
when there are few events, the error decreases slowly as the grid is
refined. The error estimate can be reduced by increasing
`control$n_inner`.

The script `inst/validation/validate_convergence.R` compares the default
grid with a much finer grid for the following data sets. For L1, whose
posterior of \\\sigma_d\\ has much mass near zero, both are also
compared with an independent reference computed in Python
(`inst/validation/python/ref_llbern_breakpoints.py`). The column Error
estimate is the error estimate reported with the default grid.

| Data | Comparison | Quantity | Reference | bssd | Difference | Error estimate |
|:---|:---|:---|---:|---:|---:|---:|
| L1: n = 150, events = 9 | finer grid vs default | P(RR(t_star) \> 1) | 0.610905 | 0.610078 | 8.3e-04 | 1.1e-03 |
| L1: n = 150, events = 9 | finer grid vs default | study drug incidence | 0.068460 | 0.068525 | 6.6e-05 |  |
| L1: n = 150, events = 9 | Python reference vs default | P(RR(t_star) \> 1) | 0.610957 | 0.610078 | 8.8e-04 | 1.1e-03 |
| L1: n = 150, events = 9 | Python reference vs default | study drug incidence | 0.068470 | 0.068525 | 5.6e-05 |  |
| L1: n = 150, events = 9 | Python reference vs finer grid | P(RR(t_star) \> 1) | 0.610957 | 0.610905 | 5.1e-05 |  |
| L2: n = 200, events = 30 | finer grid vs default | P(RR(t_star) \> 1) | 0.711812 | 0.711820 | 8.1e-06 | 2.5e-04 |
| L2: n = 200, events = 30 | finer grid vs default | study drug incidence | 0.269389 | 0.269390 | 3.2e-07 |  |
| L3: n = 100, events = 3 | finer grid vs default | P(RR(t_star) \> 1) | 0.349256 | 0.349250 | 5.7e-06 | 2.7e-04 |
| L3: n = 100, events = 3 | finer grid vs default | study drug incidence | 0.058325 | 0.058272 | 5.3e-05 |  |
| L4: n = 375, events = 20 | finer grid vs default | P(RR(t_star) \> 1) | 0.646184 | 0.646164 | 2.0e-05 | 1.0e-03 |
| L4: n = 375, events = 20 | finer grid vs default | study drug incidence | 0.066850 | 0.066850 | 1.3e-07 |  |

The largest difference is 8.8e-04, and it does not exceed the error
estimate; all 11 items are PASS.

## Comparison with MCMC

The script `inst/validation/validate_jags.R` compares
[`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
with MCMC in JAGS for the simulated data sets below. The JAGS model is
written directly from the model specification of Duncan (2026, Section
3.4) and shares no code with bssd:

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
    } 

The results in this section are read from the output of that script
(`validate_jags.csv` in the folder `validation/results` of the installed
package); JAGS is not run when this article is built.

| Data | N | k | Events | Placebo (mu, sigma) | Study drug (mu, sigma) | Incidences at t\* | t\* (week) | Analysis (week) | Prior (mu_hp, sigma_hp) | rr0 |
|:---|---:|---:|:---|:---|:---|:---|---:|---:|:---|---:|
| L1 | 150 | 1 | 9 | (7.3, 0.8) | (6.1, 0.8) | 3.5% / 13.9% | 104 | 117.0 | (7.3, 0.8) | 1 |
| L2 | 200 | 1 | 30 | (6.17, 1.6) | (6.17, 1.6) | 20.0% / 20.0% | 52 | 34.7 | (6.17, 1.6) | 1 |
| L3 | 100 | 1 | 3 | (7.3, 0.8) | (7.3, 0.8) | 3.5% / 3.5% | 104 | 78.0 | (7.3, 0.8) | 1 |

The data are simulated with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)
with the placebo and study drug parameters, and the relative risk and
the incidences are evaluated at the end of follow-up (\\t^\*\\); the
script gives the other settings. JAGS was run with 4 chains, 5,000
burn-in iterations and 50,000 iterations per chain. The Monte Carlo
standard error (SE) of each MCMC estimate uses the effective sample size
of the draws
([`coda::effectiveSize()`](https://rdrr.io/pkg/coda/man/effectiveSize.html));
for the posterior probabilities, the effective sample sizes are between
38,816 and 75,628. For the median, the table shows the MCMC proportion
of draws below the posterior median of \\RR(t^\*)\\ computed by bssd,
which should be 0.5. The incidences are the posterior means of the
cumulative incidences at \\t^\*\\. The column Ratio is the absolute
difference divided by the SE.

| Data | Quantity              |   bssd |   MCMC | Difference |      SE | Ratio |
|:-----|:----------------------|-------:|-------:|-----------:|--------:|------:|
| L1   | P(RR(t_star) \> 1)    | 0.6101 | 0.6131 |   -3.0e-03 | 2.5e-03 |  1.20 |
| L1   | P(RR \<= bssd median) | 0.5000 | 0.4975 |    2.5e-03 | 2.3e-03 |  1.10 |
| L1   | placebo incidence     | 0.0482 | 0.0481 |    1.3e-04 | 1.2e-04 |  1.04 |
| L1   | study drug incidence  | 0.0685 | 0.0685 |   -1.6e-05 | 2.7e-04 |  0.06 |
| L2   | P(RR(t_star) \> 1)    | 0.7118 | 0.7112 |    6.4e-04 | 1.7e-03 |  0.37 |
| L2   | P(RR \<= bssd median) | 0.5000 | 0.4998 |    2.5e-04 | 2.0e-03 |  0.13 |
| L2   | placebo incidence     | 0.2049 | 0.2048 |    1.2e-04 | 1.6e-04 |  0.71 |
| L2   | study drug incidence  | 0.2694 | 0.2693 |    1.3e-04 | 5.1e-04 |  0.25 |
| L3   | P(RR(t_star) \> 1)    | 0.3492 | 0.3490 |    2.6e-04 | 1.7e-03 |  0.15 |
| L3   | P(RR \<= bssd median) | 0.5000 | 0.4995 |    5.0e-04 | 1.9e-03 |  0.26 |
| L3   | placebo incidence     | 0.0505 | 0.0505 |   -2.5e-05 | 8.6e-05 |  0.29 |
| L3   | study drug incidence  | 0.0583 | 0.0583 |   -1.1e-05 | 5.1e-04 |  0.02 |

In all 12 items, the bssd value is within 1.2 standard errors of the
MCMC value, and 12 items are within two standard errors. The largest
absolute difference in a posterior probability is 0.0030, while the SEs
of the MCMC probabilities are between 0.0017 and 0.0025. The integration
error estimates of bssd for these data sets are between 2.5e-04 and
1.1e-03. Differences of this size are expected from the Monte Carlo
error of MCMC alone, so the difference between the numerical integration
of bssd and MCMC is negligible in practice. The validation script judges
an item PASS when the difference is at most four standard errors; 12 of
the 12 compared items are PASS.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227-244.
<https://doi.org/10.1080/10543406.2024.2395532>
