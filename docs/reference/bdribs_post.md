# Posterior probability of a safety signal under the BDRIBS model

Computes the posterior distribution of the relative risk under the
Poisson-based BDRIBS model of Mukhopadhyay et al. (2018) as specified in
Duncan (2026). With the total exposure time as `exposure`, the model is
BDRIBS(1); with the total time at risk (time after the first occurrence
of the event removed), it is BDRIBS(2). Both can be computed from
subject-level data with
[`time_at_risk`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md).

## Usage

``` r
bdribs_post(
  y,
  exposure,
  k = 1,
  lambda_hp,
  pt_hp,
  df = 1,
  p_prior = c(0.5, 0.5),
  rr0 = 1,
  cutoff = 0.95,
  t_star = NULL,
  level = 0.95,
  details = TRUE
)
```

## Arguments

- y:

  Number of subjects with the event in the blinded data (non-negative
  integer).

- exposure:

  Total exposure time (BDRIBS(1)) or total time at risk (BDRIBS(2)) of
  the blinded population (positive).

- k:

  Randomization ratio of study drug to placebo (`k:1`).

- lambda_hp:

  Event rate in the historical placebo population.

- pt_hp:

  Total patient time at risk in the historical placebo population, or
  `Inf` for a known background rate `lambda_hp`.

- df:

  Discount factor for the historical data (`df >= 1` discounts the
  historical information).

- p_prior:

  Shape parameters `c(a, b)` of the Beta prior for \\p\\. The default
  `c(0.5, 0.5)` is the prior used by Duncan (2026).

- rr0:

  Critical value of the rate ratio (\\\theta_c\\ in Duncan, 2026).

- cutoff:

  Posterior probability threshold for a signal.

- t_star:

  Optional time at which the cumulative incidences are summarized.

- level:

  Probability of the equal-tailed credible interval of the rate ratio.

- details:

  Logical; if `FALSE`, only `prob`, `signal` and `incidence` are
  computed, and `rr` and `estimates` are `NA`. This is faster and is
  used in
  [`sim_oc`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md).

## Value

An object of class `bssd_post`, a list with elements

- `model`:

  `"BDRIBS"`.

- `prob`:

  Posterior probability \\P(r \> rr_0 \mid Y)\\.

- `signal`:

  Logical; `prob > cutoff`.

- `rr`:

  Posterior mean, median, and the limits of the equal-tailed credible
  interval of \\r\\.

- `estimates`:

  Posterior means of \\p\\, \\r\\, the placebo event rate \\\delta\\
  (`lambda_p`) and the study drug event rate \\r \delta\\ (`lambda_d`).

- `incidence`:

  Posterior means of the cumulative incidences at `t_star` for placebo
  and study drug and of their difference, or `NULL` when `t_star` is not
  given.

- `rr_label`, `rr0`, `cutoff`, `level`, `t_star`:

  Settings used.

- `data`, `prior`:

  Data summary and prior parameters.

- `integration`:

  Description of the numerical integration.

## Details

The model is \$\$Y \sim \mathrm{Poisson}\[E \delta (r k + 1) / (k +
1)\],\$\$ \$\$\delta \sim \mathrm{Gamma}(\lambda\_{HP} PT\_{HP} / DF,
PT\_{HP} / DF),\$\$ \$\$r = p / \[k (1 - p)\], \quad p \sim
\mathrm{Beta}(a, b),\$\$ where \\Y\\ is the number of subjects with the
event, \\E\\ is the exposure, \\k\\ is the randomization ratio (study
drug to placebo \\k:1\\), \\\delta\\ is the event rate in the placebo
group, and \\r\\ is the rate ratio of study drug to placebo. The Gamma
distribution is parameterized by shape and rate.

Because the Gamma prior is conjugate to the Poisson likelihood,
\\\delta\\ is integrated out analytically, and the marginal posterior
density of \\p\\ is proportional to \$\$p^{a - 1} (1 - p)^{b + \alpha -
1} \[E + \beta (k + 1) (1 - p)\]^{-(\alpha + Y)}\$\$ with \\\alpha =
\lambda\_{HP} PT\_{HP} / DF\\ and \\\beta = PT\_{HP} / DF\\. All
posterior quantities are one-dimensional integrals over
\\\mathrm{logit}(p)\\, which are evaluated by adaptive quadrature
([`integrate`](https://rdrr.io/r/stats/integrate.html)) with a relative
tolerance of `1e-10`.

A signal is declared when \\P(r \> rr_0 \mid Y) \> \\ `cutoff` (Duncan,
2026, Equation 3).

When `t_star` is given, the posterior means of the cumulative incidences
at `t_star` are computed under the exponential distribution implied by
the constant event rates, that is, \\1 - \exp(-\delta t^\*)\\ for
placebo and \\1 - \exp(-r \delta t^\*)\\ for study drug.

With `pt_hp = Inf`, the prior of \\\delta\\ is a point mass at
`lambda_hp`, so the inference is conditional on a known background rate.
This is the screening analysis of Mukhopadhyay et al. (2018) and
corresponds to `inf.type = 1` of the bdribs package; the marginal
posterior density of \\p\\ is then proportional to \$\$p^{a - 1} (1 -
p)^{b - Y - 1} \exp\[-E \lambda\_{HP} / ((k + 1) (1 - p))\].\$\$

The time unit of `exposure`, `lambda_hp`, `pt_hp` and `t_star` must be
the same.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
detection of potential risk using inference on blinded safety data.
*Pharmaceutical Statistics*, 17(6), 823–834.
[doi:10.1002/pst.1898](https://doi.org/10.1002/pst.1898)

## See also

[`time_at_risk`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md),
[`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md),
[`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)

## Examples

``` r
# 5 subjects with the event in 500 patient-years, randomization 2:1,
# historical placebo rate 18 / 4000 per patient-year
fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                   pt_hp = 4000)
fit
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.8199
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 2.451, median 2.123, 95% interval (0.195, 6.588)

# known background rate of 0.011 per patient-year and uniform prior for p
# (hypothetical example 2 of Mukhopadhyay et al., 2018)
bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011, pt_hp = Inf,
            p_prior = c(1, 1))
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.753
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.522, median 1.438, 95% interval (0.379, 3.143)
```
