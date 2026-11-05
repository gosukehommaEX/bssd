# EXP_BERN: exponential Bernoulli mixture model

``` r

library(bssd)
```

## Model

The exponential Bernoulli (EXP_BERN) model of Duncan (2026, Section 3.3)
assumes constant hazard rates, like BDRIBS, but models whether each
subject had the event during the individual follow-up time instead of
the number of subjects with the event. For subject \\i\\ with the
follow-up time \\t_i\\, the event indicator \\Y_i\\ follows the mixture
of two Bernoulli distributions, \\P(Y_i = 1) = w_p \\1 - \exp(-\lambda_p
t_i)\\ + w_d \\1 - \exp(-r \lambda_p t_i)\\,\\ where \\w_p = 1 / (k +
1)\\ and \\w_d = k / (k + 1)\\ are the proportions of subjects on
placebo and on the study drug for the randomization ratio \\k:1\\,
\\\lambda_p\\ is the hazard rate on placebo and \\r\\ is the rate ratio
of the study drug to placebo. The priors are \\\lambda_p \sim
\mathrm{Gamma}(\lambda\_{HP} PT\_{HP} / DF, PT\_{HP} / DF), \quad r =
\frac{p}{k (1 - p)}, \quad p \sim \mathrm{Beta}(a, b),\\ with the shape
and rate parameterization of the Gamma distribution, as in BDRIBS
([`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md)).
Duncan (2026) uses the Beta(0.5, 0.5) prior for \\p\\. A signal is
declared when \\P(r \> \theta_c \mid \mathrm{data}) \>
P\_{\mathrm{cutoff}}\\.

The arguments of
[`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
correspond to the symbols as follows.

| Symbol | Argument | Description |
|----|----|----|
| \\Y_i\\ | `y` | event indicator of each subject (1 if the subject had the event during follow-up) |
| \\t_i\\ | `t` | follow-up time of each subject |
| \\k\\ | `k` | randomization ratio of study drug to placebo |
| \\\lambda\_{HP}\\ | `lambda_hp` | event rate in the historical placebo data |
| \\PT\_{HP}\\ | `pt_hp` | total patient time of the historical placebo data |
| \\DF\\ | `df` | discount factor |
| \\(a, b)\\ | `p_prior` | shape parameters of the Beta prior of \\p\\ (default c(0.5, 0.5)) |
| \\\theta_c\\ | `rr0` | critical value of \\r\\ (default 1) |
| \\P\_{\mathrm{cutoff}}\\ | `cutoff` | posterior probability threshold (default 0.95) |

`t`, `lambda_hp`, `pt_hp` and the optional time `t_star` for the
cumulative incidences must be in the same time unit. Because the
likelihood stays positive as \\r \to \infty\\, the posterior mean of
\\r\\ is infinite when the second shape parameter of the Beta prior is 1
or less, as for the default prior; the posterior median is then the
summary of \\r\\ to use.

## Example

The following data are simulated with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md):
300 subjects randomized 1:1 with uniform entry over 26 weeks, 52 weeks
of planned follow-up, and exponential times to the event with cumulative
incidences of 10% on placebo and 20% on the study drug at week 52. At
the blinded analysis at week 52, the follow-up times are between 26 and
52 weeks.

``` r

set.seed(2026)
rate_of <- function(ip, t) -log(1 - ip) / t
d <- sim_blinded(n = 300, k = 1, accrual = 26, follow_up = 52,
                 analysis_time = 52, dist = "exponential",
                 par_p = rate_of(0.10, 52), par_d = rate_of(0.20, 52))
head(d[, c("id", "follow_up", "event")])
#>   id follow_up event
#> 1  1  44.36873     0
#> 2  2  42.33813     0
#> 3  3  42.56790     0
#> 4  4  42.60844     0
#> 5  5  31.83289     0
#> 6  6  31.65294     0
```

The prior of \\\lambda_p\\ is based on a historical placebo group of 150
subjects followed for 52 weeks with the incidence of 10%, discounted by
\\DF = 2\\.

``` r

fit <- expbern_post(y = d$event, t = d$follow_up, k = 1,
                    lambda_hp = rate_of(0.10, 52), pt_hp = 150 * 52,
                    df = 2, t_star = 52)
fit
#> Model: EXP_BERN 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.751
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean Inf, median 1.626, 95% interval (0.2321, 5.326)
#> Cumulative incidence at t_star = 52: placebo 0.1094, study drug 0.1638, difference 0.05434
summary(fit)
#>      model      prob signal rr_mean rr_median  rr_lower rr_upper ip_placebo
#> 1 EXP_BERN 0.7510482  FALSE     Inf  1.625838 0.2321424 5.326324  0.1094196
#>     ip_drug ip_difference        error
#> 1 0.1637567    0.05433711 5.669898e-11
```

The posterior probabilities for other critical values of \\r\\ are
obtained with `rr0`; `details = FALSE` skips the quantiles of \\r\\:

``` r

vapply(c(1, 1.5, 2), function(q) {
  expbern_post(y = d$event, t = d$follow_up, k = 1,
               lambda_hp = rate_of(0.10, 52), pt_hp = 150 * 52, df = 2,
               rr0 = q, details = FALSE)$prob
}, numeric(1))
#> [1] 0.7510482 0.5476504 0.3748275
```

## Numerical integration

[`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
integrates the two-dimensional posterior numerically on the coordinates
\\\log \lambda_p\\ and the prior distribution function of \\p\\. A
coarse grid locates the region that holds almost all posterior mass;
Gauss-Legendre quadrature is used for \\\log \lambda_p\\ and Chebyshev
interpolation for the second coordinate (see
[`?expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)).
With `control$check = TRUE`, the default, the computation is repeated on
a finer grid, the finer result is returned, and the absolute difference
in `prob` is reported as the error estimate (`fit$integration$error`,
here 5.7e-11).

The script `inst/validation/validate_convergence.R` compares the default
grid with a much finer grid for the following data sets. The column
Error estimate is the error estimate reported with the default grid.

| Data | Comparison | Quantity | Reference | bssd | Difference | Error estimate |
|:---|:---|:---|---:|---:|---:|---:|
| E1: n = 150, events = 15 | finer grid vs default | P(r \> 1) | 0.541398 | 0.541398 | 2.4e-09 | 2.0e-11 |
| E1: n = 150, events = 15 | finer grid vs default | study drug incidence | 0.126419 | 0.126419 | 6.2e-10 |  |
| E2: n = 300, events = 101 | finer grid vs default | P(r \> 1) | 0.714322 | 0.714322 | 1.1e-10 | 5.7e-12 |
| E2: n = 300, events = 101 | finer grid vs default | study drug incidence | 0.365416 | 0.365416 | 3.3e-10 |  |
| E3: n = 200, events = 4 | finer grid vs default | P(r \> 1.5) | 0.121202 | 0.121202 | 5.5e-10 | 1.2e-09 |
| E3: n = 200, events = 4 | finer grid vs default | study drug incidence | 0.026955 | 0.026955 | 3.8e-11 |  |
| E4: n = 60, events = 0 | finer grid vs default | P(r \> 1) | 0.396049 | 0.396050 | 2.1e-07 | 5.7e-06 |
| E4: n = 60, events = 0 | finer grid vs default | study drug incidence | 0.008988 | 0.008988 | 9.9e-09 |  |

The largest difference is 2.1e-07, and all 8 items are PASS.

## Comparison with MCMC

The script `inst/validation/validate_jags.R` compares
[`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
with MCMC in JAGS for the simulated data sets below. The JAGS model is
written directly from the model specification of Duncan (2026, Section
3.3) and shares no code with bssd:

    model {
      for (i in 1:N) {
        pp[i] <- 1 - exp(-lambda * t[i])
        pd[i] <- 1 - exp(-r * lambda * t[i])
        y[i] ~ dbern(wp * pp[i] + (1 - wp) * pd[i])
      }
      lambda ~ dgamma(alpha, beta)
      p ~ dbeta(a, b)
      r <- p / (k * (1 - p))
    } 

The results in this section are read from the output of that script
(`validate_jags.csv` in the folder `validation/results` of the installed
package); JAGS is not run when this article is built.

| Data | N | k | Events | Placebo (week 52) | Study drug (week 52) | Analysis (week) | Prior incidence (week 52) | Prior of p | rr0 |
|:---|---:|---:|:---|:---|:---|---:|:---|:---|---:|
| E1 | 150 | 1 | 15 | 10% | 20% | 52 | 10% | Beta(0.5, 0.5) | 1.0 |
| E2 | 300 | 1 | 101 | 30% | 30% | 78 | 30% | Beta(0.5, 0.5) | 1.0 |
| E3 | 200 | 2 | 4 | 5% | 10% | 40 | 5% | Beta(1, 1) | 1.5 |

Placebo and Study drug are the incidences used to simulate the data with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md),
and Prior incidence is the incidence from which `lambda_hp` is computed;
the script gives the other settings. JAGS was run with 4 chains, 5,000
burn-in iterations and 50,000 iterations per chain. The Monte Carlo
standard error (SE) of each MCMC estimate uses the effective sample size
of the draws
([`coda::effectiveSize()`](https://rdrr.io/pkg/coda/man/effectiveSize.html));
for the posterior probabilities, the effective sample sizes are between
21,753 and 59,895. For the median, the table shows the MCMC proportion
of draws below the posterior median computed by bssd, which should be
0.5. The incidences are the posterior means of the cumulative incidences
at week 52. The column Ratio is the absolute difference divided by the
SE.

| Data | Quantity             |   bssd |   MCMC | Difference |      SE | Ratio |
|:-----|:---------------------|-------:|-------:|-----------:|--------:|------:|
| E1   | P(r \> 1)            | 0.5414 | 0.5363 |    5.1e-03 | 3.4e-03 |  1.50 |
| E1   | P(r \<= bssd median) | 0.5000 | 0.5049 |   -4.9e-03 | 3.4e-03 |  1.45 |
| E1   | placebo incidence    | 0.1174 | 0.1179 |   -4.8e-04 | 3.7e-04 |  1.28 |
| E1   | study drug incidence | 0.1264 | 0.1257 |    6.8e-04 | 5.2e-04 |  1.32 |
| E2   | P(r \> 1)            | 0.7143 | 0.7146 |   -3.2e-04 | 3.0e-03 |  0.11 |
| E2   | P(r \<= bssd median) | 0.5000 | 0.4976 |    2.4e-03 | 3.4e-03 |  0.70 |
| E2   | placebo incidence    | 0.3038 | 0.3037 |    1.4e-04 | 3.9e-04 |  0.36 |
| E2   | study drug incidence | 0.3654 | 0.3656 |   -2.2e-04 | 5.1e-04 |  0.44 |
| E3   | P(r \> 1.5)          | 0.1212 | 0.1205 |    6.7e-04 | 1.3e-03 |  0.51 |
| E3   | P(r \<= bssd median) | 0.5000 | 0.5008 |   -8.2e-04 | 2.0e-03 |  0.41 |
| E3   | placebo incidence    | 0.0533 | 0.0533 |   -5.4e-05 | 1.2e-04 |  0.46 |
| E3   | study drug incidence | 0.0270 | 0.0269 |    2.2e-05 | 7.5e-05 |  0.30 |

In all 12 items, the bssd value is within 1.5 standard errors of the
MCMC value, and 12 items are within two standard errors. The largest
absolute difference in a posterior probability is 0.0051, while the SEs
of the MCMC probabilities are between 0.0013 and 0.0034. Differences of
this size are expected from the Monte Carlo error of MCMC alone, so the
difference between the numerical integration of bssd and MCMC is
negligible in practice. The validation script judges an item PASS when
the difference is at most four standard errors; 12 of the 12 items are
PASS.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227-244.
<https://doi.org/10.1080/10543406.2024.2395532>
