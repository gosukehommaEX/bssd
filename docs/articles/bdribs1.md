# BDRIBS(1): Poisson model with the total exposure time

``` r

library(bssd)
```

## Model

BDRIBS (Bayesian Detection of potential Risk using Inference on Blinded
Safety data) is the model of Mukhopadhyay et al. (2018). Duncan (2026,
Section 3.2.2) calls it BDRIBS(1) when the exposure is the total
exposure time of the blinded population. The number of subjects with the
event in the blinded data, \\Y\\, follows \\Y \sim \mathrm{Poisson}\\E
\delta (r k + 1) / (k + 1)\\,\\ where \\E\\ is the total exposure time,
\\\delta\\ is the event rate on placebo, \\r\\ is the rate ratio of the
study drug to placebo, and the study drug and placebo are randomized in
the ratio \\k:1\\. The priors are \\\delta \sim
\mathrm{Gamma}(\lambda\_{HP} PT\_{HP} / DF, PT\_{HP} / DF), \quad r =
\frac{p}{k (1 - p)}, \quad p \sim \mathrm{Beta}(a, b),\\ where the Gamma
distribution has the shape and rate parameters, \\\lambda\_{HP}\\ and
\\PT\_{HP}\\ are the event rate and the total patient time of historical
placebo data, and \\DF\\ is a discount factor for the historical data.
The parameter \\p\\ is the probability that an event in the blinded data
occurred on the study drug. A signal is declared when \\P(r \> \theta_c
\mid Y) \> P\_{\mathrm{cutoff}}\\.

The arguments of
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
correspond to the symbols as follows.

| Symbol | Argument | Description |
|----|----|----|
| \\Y\\ | `y` | number of subjects with the event |
| \\E\\ | `exposure` | total exposure time, the sum of the follow-up times |
| \\k\\ | `k` | randomization ratio of study drug to placebo |
| \\\lambda\_{HP}\\ | `lambda_hp` | event rate in the historical placebo data |
| \\PT\_{HP}\\ | `pt_hp` | total patient time of the historical placebo data; `Inf` fixes \\\delta\\ at `lambda_hp` |
| \\DF\\ | `df` | discount factor |
| \\(a, b)\\ | `p_prior` | shape parameters of the Beta prior of \\p\\ (default c(0.5, 0.5)) |
| \\\theta_c\\ | `rr0` | critical value of \\r\\ (default 1) |
| \\P\_{\mathrm{cutoff}}\\ | `cutoff` | posterior probability threshold (default 0.95) |

The exposure, `lambda_hp`, `pt_hp` and the optional time `t_star` for
the cumulative incidences must be in the same time unit.

## Example: hypothetical example 2 of Mukhopadhyay et al. (2018)

In Section 4.4 of Mukhopadhyay et al. (2018), the study drug and placebo
are randomized 2:1, the blinded data have 540 patient-years of exposure,
and the background rate on placebo is about 0.011 per patient-year,
estimated from 33 events in about 3000 patient-years of historical data.

### Screening analysis

The screening analysis treats the background rate as known
(`pt_hp = Inf`) and uses the uniform prior for \\p\\. Mukhopadhyay et
al. (2018, Sections 3.4 and 4.4) compare the posterior probability with
the boundary \\3 / (3 + k \theta_c)\\, which is used here as `cutoff`.

``` r

screen <- expand.grid(rr0 = c(1, 1.3), y = c(8, 9))
screen$cutoff <- 3 / (3 + 2 * screen$rr0)
screen$prob <- mapply(function(y, rr0, cutoff) {
  bdribs_post(y = y, exposure = 540, k = 2, lambda_hp = 0.011, pt_hp = Inf,
              p_prior = c(1, 1), rr0 = rr0, cutoff = cutoff,
              details = FALSE)$prob
}, screen$y, screen$rr0, screen$cutoff)
screen$signal <- screen$prob > screen$cutoff
```

| Events | rr0 | Boundary |  bssd | Signal | Mukhopadhyay et al. |
|-------:|----:|---------:|------:|:-------|--------------------:|
|      8 | 1.0 |    0.600 | 0.619 | TRUE   |               0.622 |
|      8 | 1.3 |    0.536 | 0.433 | FALSE  |               0.434 |
|      9 | 1.0 |    0.600 | 0.753 | TRUE   |               0.756 |
|      9 | 1.3 |    0.536 | 0.580 | TRUE   |               0.584 |

The published values were computed by MCMC with the bdribs package, and
they differ from the values of bssd by at most 0.004. With 9 events, the
posterior is

``` r

fit9 <- bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011,
                    pt_hp = Inf, p_prior = c(1, 1))
fit9
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.753
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.522, median 1.438, 95% interval (0.379, 3.143)
```

### Sensitivity analysis

In the second step, the screening analysis is repeated with background
rates at the limits of the 90% credible interval of the background rate.
As in Section 4.3 of the article, the interval is computed from the
Gamma distribution with the number of historical events plus 0.5 as the
shape and the historical patient-years plus 0.5 as the rate.

``` r

rate_ci <- qgamma(c(0.05, 0.95), shape = 33 + 0.5, rate = 3000 + 0.5)
rate_ci
#> [1] 0.008192346 0.014515593
sens <- expand.grid(rr0 = c(1, 1.3), lambda_hp = rate_ci)
sens$prob <- mapply(function(lambda_hp, rr0) {
  bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = lambda_hp,
              pt_hp = Inf, p_prior = c(1, 1), rr0 = rr0,
              details = FALSE)$prob
}, sens$lambda_hp, sens$rr0)
```

| Background rate | rr0 |  bssd | Mukhopadhyay et al. |
|----------------:|----:|------:|--------------------:|
|          0.0082 | 1.0 | 0.920 |                0.92 |
|          0.0082 | 1.3 | 0.833 |                0.83 |
|          0.0145 | 1.0 | 0.478 |                0.48 |
|          0.0145 | 1.3 | 0.280 |                0.28 |

The interval (0.0082, 0.0145) agrees with the interval (0.0082, 0.0145)
of Section 4.4 of the article. The posterior probabilities differ from
the published values, which are given to two decimal places, by at most
0.003.

Instead of fixing the background rate, the historical data can be used
as the Gamma prior of \\\delta\\, which gives the posterior probability
unconditional on the background rate (Mukhopadhyay et al., 2018, Section
4.3):

``` r

fit_gamma <- bdribs_post(y = 9, exposure = 540, k = 2,
                         lambda_hp = (33 + 0.5) / (3000 + 0.5),
                         pt_hp = 3000 + 0.5, df = 1, p_prior = c(1, 1))
fit_gamma
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.7076
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.494, median 1.373, 95% interval (0.3078, 3.371)
```

### Total exposure time from subject-level data

For subject-level data,
[`time_at_risk()`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md)
returns the total exposure time (`exposure`), which BDRIBS(1) uses, and
the total time at risk, which BDRIBS(2) uses
([`vignette("bdribs2")`](https://gosukehommaex.github.io/bssd/articles/bdribs2.md)).

``` r

time_at_risk(follow_up = c(52, 52, 30, 12), event = c(1, 0, 1, 0),
             event_time = c(10, NA, 25, NA))
#>     exposure time_at_risk      n_event 
#>          146           99            2
```

## Numerical integration

Because the Gamma prior of \\\delta\\ is conjugate to the Poisson
likelihood,
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
integrates \\\delta\\ out analytically. The remaining posterior
quantities are one-dimensional integrals over \\\mathrm{logit}(p)\\,
which are computed by adaptive quadrature
([`stats::integrate()`](https://rdrr.io/r/stats/integrate.html)) with a
relative tolerance of \\10^{-10}\\ (see
[`?bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)).
The expected values of the tests in `tests/testthat` were computed
beforehand by independent Python implementations
(`inst/validation/python/ref_bdribs.py` and `ref_bdribs_fixed.py`).

## Comparison with MCMC

The bdribs package of Mukhopadhyay et al. (2018) samples the same model
by MCMC with JAGS. The script
`inst/validation/validate_bdribs_package.R` compares
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with [`bdribs::bdribs()`](https://rdrr.io/pkg/bdribs/man/bdribs.html)
for the settings below. The results in this section are read from the
output of that script (`validate_bdribs_package.csv` in the folder
`validation/results` of the installed package); JAGS is not run when
this article is built.

| Data | y | E | k | Historical events | Historical patient-years | Background rate | Prior of p |
|:---|---:|---:|---:|---:|---:|:---|:---|
| B1 | 5 | 500 | 2 | 18 | 4000 | fixed (inf.type = 1) | Beta(1, 1) |
| B2 | 5 | 500 | 2 | 18 | 4000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B3 | 5 | 500 | 2 | 18 | 4000 | fixed (inf.type = 1) | Beta(0.5, 0.5) |
| B4 | 9 | 540 | 2 | 33 | 3000 | fixed (inf.type = 1) | Beta(1, 1) |
| B5 | 9 | 540 | 2 | 33 | 3000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B6 | 7 | 800 | 1 | 18 | 4000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B7 | 30 | 2000 | 1 | 20 | 2000 | Gamma prior (inf.type = 2) | Beta(0.5, 0.5) |

With `inf.type = 1`,
[`bdribs::bdribs()`](https://rdrr.io/pkg/bdribs/man/bdribs.html) fixes
the background rate at the historical events divided by the historical
patient-years, which is
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with `pt_hp = Inf`. With `inf.type = 2`, it uses the Gamma prior with
the shape (historical events + 0.001) and the rate (historical
patient-years + 0.001), which is
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with these values as `lambda_hp * pt_hp` and `pt_hp` and `df = 1`.
[`bdribs::bdribs()`](https://rdrr.io/pkg/bdribs/man/bdribs.html) was run
with 2 chains, 5,000 burn-in iterations and 100,000 iterations (the
arguments `nc`, `burn` and `iter` of `mc.params`). The Monte Carlo
standard error (SE) of each MCMC estimate uses the effective sample size
of the draws
([`coda::effectiveSize()`](https://rdrr.io/pkg/coda/man/effectiveSize.html));
for the posterior probabilities, the effective sample sizes are between
46,649 and 115,901. For the median, the table shows the MCMC proportion
of draws below the posterior median computed by bssd, which should be
0.5. The column Ratio is the absolute difference divided by the SE.

| Data | Quantity             |   bssd |   MCMC | Difference |      SE | Ratio |
|:-----|:---------------------|-------:|-------:|-----------:|--------:|------:|
| B1   | P(r \> 1)            | 0.8154 | 0.8163 |   -8.9e-04 | 1.6e-03 |  0.55 |
| B1   | P(r \> 1.2)          | 0.7523 | 0.7526 |   -2.8e-04 | 1.7e-03 |  0.16 |
| B1   | P(r \> 1.5)          | 0.6520 | 0.6518 |    1.9e-04 | 1.8e-03 |  0.11 |
| B1   | P(r \<= bssd median) | 0.5000 | 0.5006 |   -5.5e-04 | 1.7e-03 |  0.33 |
| B1   | mean of r            | 2.1834 | 2.1842 |   -7.7e-04 | 4.4e-03 |  0.17 |
| B2   | P(r \> 1)            | 0.7850 | 0.7859 |   -9.3e-04 | 1.7e-03 |  0.55 |
| B2   | P(r \> 1.2)          | 0.7195 | 0.7202 |   -7.4e-04 | 1.8e-03 |  0.41 |
| B2   | P(r \> 1.5)          | 0.6196 | 0.6209 |   -1.3e-03 | 1.8e-03 |  0.69 |
| B2   | P(r \<= bssd median) | 0.5000 | 0.4985 |    1.5e-03 | 1.8e-03 |  0.85 |
| B2   | mean of r            | 2.1907 | 2.1929 |   -2.2e-03 | 5.5e-03 |  0.40 |
| B3   | P(r \> 1)            | 0.8461 | 0.8465 |   -4.6e-04 | 1.6e-03 |  0.28 |
| B3   | P(r \> 1.2)          | 0.7933 | 0.7930 |    3.0e-04 | 1.8e-03 |  0.17 |
| B3   | P(r \> 1.5)          | 0.7056 | 0.7052 |    3.7e-04 | 1.9e-03 |  0.20 |
| B3   | P(r \<= bssd median) | 0.5000 | 0.5005 |   -5.4e-04 | 1.7e-03 |  0.31 |
| B3   | mean of r            | 2.3919 | 2.3897 |    2.2e-03 | 5.0e-03 |  0.44 |
| B4   | P(r \> 1)            | 0.7530 | 0.7550 |   -2.0e-03 | 1.5e-03 |  1.30 |
| B4   | P(r \> 1.2)          | 0.6392 | 0.6411 |   -1.9e-03 | 1.6e-03 |  1.21 |
| B4   | P(r \> 1.5)          | 0.4647 | 0.4670 |   -2.3e-03 | 1.5e-03 |  1.57 |
| B4   | P(r \<= bssd median) | 0.5000 | 0.4975 |    2.5e-03 | 1.5e-03 |  1.67 |
| B4   | mean of r            | 1.5218 | 1.5259 |   -4.1e-03 | 2.2e-03 |  1.84 |
| B5   | P(r \> 1)            | 0.7192 | 0.7181 |    1.1e-03 | 1.7e-03 |  0.66 |
| B5   | P(r \> 1.2)          | 0.6089 | 0.6072 |    1.8e-03 | 1.7e-03 |  1.04 |
| B5   | P(r \> 1.5)          | 0.4485 | 0.4464 |    2.0e-03 | 1.6e-03 |  1.27 |
| B5   | P(r \<= bssd median) | 0.5000 | 0.5015 |   -1.5e-03 | 1.7e-03 |  0.87 |
| B5   | mean of r            | 1.5231 | 1.5205 |    2.6e-03 | 2.9e-03 |  0.90 |
| B6   | P(r \> 1)            | 0.8152 | 0.8130 |    2.2e-03 | 1.5e-03 |  1.46 |
| B6   | P(r \> 1.2)          | 0.7580 | 0.7558 |    2.3e-03 | 1.6e-03 |  1.40 |
| B6   | P(r \> 1.5)          | 0.6679 | 0.6644 |    3.6e-03 | 1.7e-03 |  2.05 |
| B6   | P(r \<= bssd median) | 0.5000 | 0.5028 |   -2.8e-03 | 1.7e-03 |  1.61 |
| B6   | mean of r            | 2.3795 | 2.3646 |    1.5e-02 | 5.8e-03 |  2.57 |
| B7   | P(r \> 1)            | 0.8777 | 0.8752 |    2.5e-03 | 1.5e-03 |  1.68 |
| B7   | P(r \> 1.2)          | 0.7999 | 0.7971 |    2.9e-03 | 1.8e-03 |  1.57 |
| B7   | P(r \> 1.5)          | 0.6582 | 0.6548 |    3.5e-03 | 2.2e-03 |  1.58 |
| B7   | P(r \<= bssd median) | 0.5000 | 0.5017 |   -1.7e-03 | 2.3e-03 |  0.74 |
| B7   | mean of r            | 1.9489 | 1.9450 |    3.9e-03 | 4.8e-03 |  0.82 |

In all 35 items, the bssd value is within 2.6 standard errors of the
MCMC value, and 33 items are within two standard errors. The largest
absolute difference in a posterior probability is 0.0036, while the SEs
of the MCMC probabilities are between 0.0015 and 0.0023. Differences of
this size are expected from the Monte Carlo error of MCMC alone, so the
difference between the numerical integration of bssd and MCMC is
negligible in practice. The validation script judges an item PASS when
the difference is at most four standard errors; 35 of the 35 items are
PASS.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227-244.
<https://doi.org/10.1080/10543406.2024.2395532>

Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
detection of potential risk using inference on blinded safety data.
*Pharmaceutical Statistics*, 17(6), 823-834.
<https://doi.org/10.1002/pst.1898>
