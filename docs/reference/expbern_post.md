# Posterior probability of a safety signal under the EXP_BERN model

Computes the posterior distribution of the rate ratio under the
exponential Bernoulli (EXP_BERN) mixture model of Duncan (2026), which
models whether each subject had the event during the individual
follow-up time and does not assume a Poisson distribution for the number
of subjects with the event.

## Usage

``` r
expbern_post(
  y,
  t,
  k = 1,
  lambda_hp,
  pt_hp,
  df = 1,
  p_prior = c(0.5, 0.5),
  rr0 = 1,
  cutoff = 0.95,
  t_star = NULL,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Vector of event indicators (1 if the subject had at least one
  occurrence of the event during follow-up, 0 otherwise).

- t:

  Vector of follow-up times (positive), in the same time unit as
  `lambda_hp`, `pt_hp` and `t_star`.

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

  Logical; if `FALSE`, the quantiles of \\r\\ are not computed (faster;
  used in
  [`sim_oc`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)).

- control:

  List of settings for the numerical integration: `n_coarse` (cells per
  coordinate of the coarse grid, default 41), `n_lambda` (Gauss-Legendre
  nodes, default 48), `n_p` (Chebyshev points, default 49), `check`
  (default `TRUE`), `eps` (default `1e-8`) and `max_iter` (default 12).

## Value

An object of class `bssd_post`; see
[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md).
The element `model` is `"EXP_BERN"`, `estimates` contains the posterior
means of \\p\\, \\r\\, \\\lambda_p\\ and \\r \lambda_p\\ (`lambda_d`),
and `integration` contains the nodes used and the error estimate.

## Details

For subject \\i\\ with follow-up time \\t_i\\, the model is \$\$P(Y_i
= 1) = w_p \[1 - \exp(-\lambda_p t_i)\] + w_d \[1 - \exp(-r \lambda_p
t_i)\],\$\$ with \\w_p = 1 / (k + 1)\\ and \\w_d = k / (k + 1)\\, and
\$\$\lambda_p \sim \mathrm{Gamma}(\lambda\_{HP} PT\_{HP} / DF, PT\_{HP}
/ DF), \quad r = p / \[k (1 - p)\], \quad p \sim \mathrm{Beta}(a,
b).\$\$ The Gamma distribution is parameterized by shape and rate.

The two-dimensional posterior is integrated numerically on the
coordinates \\\log \lambda_p\\ and the prior distribution function of
\\p\\. The latter keeps the tails of the posterior of \\p\\ under
control: the likelihood does not vanish as \\r \to 0\\ or \\r \to
\infty\\, so these tails follow the prior. A coarse grid first locates
the region that holds all but `control$eps` of the posterior mass in
each marginal. On that region, Gauss-Legendre quadrature is used for
\\\log \lambda_p\\, and the marginal density of the second coordinate is
interpolated at Chebyshev points, which gives its distribution function
at any point (used for `prob` and the quantiles of \\r\\). When
`control$check` is `TRUE`, the computation is repeated with about 1.5
times as many nodes, the finer result is returned, and the absolute
difference in `prob` is reported as the integration error estimate.

Because the likelihood stays positive as \\r \to \infty\\, the posterior
means of \\r\\ and \\r \lambda_p\\ are infinite when the second shape
parameter of the Beta prior is 1 or less, as for the default
`c(0.5, 0.5)`. They are then returned as `Inf`, and the posterior median
of \\r\\ is the summary to use.

A signal is declared when \\P(r \> rr_0 \mid \mathrm{data}) \> \\
`cutoff`.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
[`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)

## Examples

``` r
set.seed(1)
n <- 100
t <- runif(n, 10, 52)
y <- rbinom(n, 1, 1 - exp(-0.004 * t))
fit <- expbern_post(y, t, k = 1, lambda_hp = 0.003, pt_hp = 2000, df = 2,
                    t_star = 52)
fit
#> Model: EXP_BERN 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.6589
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean Inf, median 1.573, 95% interval (0.02978, 10.1)
#> Cumulative incidence at t_star = 52: placebo 0.1739, study drug 0.2425, difference 0.06865
```
