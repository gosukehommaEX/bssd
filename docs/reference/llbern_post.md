# Posterior probability of a safety signal under the LL_BERN model

Computes the posterior distribution of the relative risk under the
log-logistic Bernoulli (LL_BERN) mixture model of Duncan (2026), which
allows the hazard of the event to change over time.

## Usage

``` r
llbern_post(
  y,
  t,
  k = 1,
  mu_hp,
  sigma_hp,
  pt_hp,
  df_p = 1,
  mu_d_prior = c(mu_hp, 1),
  sigma_d_prior = c(sigma_hp, 1),
  t_star,
  rr0 = 1,
  cutoff = 0.95,
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

  Vector of follow-up times (positive), in the time unit for which
  `mu_hp` and `sigma_hp` were estimated.

- k:

  Randomization ratio of study drug to placebo (`k:1`).

- mu_hp:

  Location parameter \\\mu\\ of the log-logistic distribution estimated
  from historical placebo data (positive).

- sigma_hp:

  Scale parameter \\\sigma\\ of the log-logistic distribution estimated
  from historical placebo data (positive).

- pt_hp:

  Total patient time of the historical placebo population.

- df_p:

  Discount factor for the historical placebo data.

- mu_d_prior:

  Shape and rate of the Gamma prior for \\\mu_d\\; the default is
  `c(mu_hp, 1)`.

- sigma_d_prior:

  Shape and rate of the Gamma prior for \\\sigma_d\\; the default is
  `c(sigma_hp, 1)`.

- t_star:

  Time at which the relative risk and the cumulative incidences are
  evaluated.

- rr0:

  Critical value of the relative risk.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible interval of the relative
  risk.

- details:

  Logical; if `FALSE`, the quantiles of \\RR(t^\*)\\ are not computed
  (faster; used in
  [`sim_oc`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)).

- control:

  List of settings for the numerical integration: `n_coarse` (cells per
  coordinate of the coarse grid, default 12), `n_outer` (Gauss-Legendre
  nodes per coordinate, default 16), `n_inner` (Chebyshev points,
  default 121), `check` (default `TRUE`), `eps` (default `1e-8`) and
  `max_iter` (default 12).

## Value

An object of class `bssd_post`; see
[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md).
The element `model` is `"LL_BERN"`, `rr` summarizes \\RR(t^\*)\\,
`estimates` contains the posterior means of \\\mu_p\\, \\\sigma_p\\,
\\\mu_d\\ and \\\sigma_d\\, and `integration` contains the nodes used
and the error estimate.

## Details

For subject \\i\\ with follow-up time \\t_i\\, the model is \$\$P(Y_i
= 1) = w_p F(t_i; \mu_p, \sigma_p) + w_d F(t_i; \mu_d, \sigma_d),\$\$
where \\F\\ is the log-logistic cumulative incidence
([`ll_cdf`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)),
\\w_p = 1 / (k + 1)\\ and \\w_d = k / (k + 1)\\. The priors are
\$\$\mu_p \sim \mathrm{Gamma}(PT\_{HP} \mu / DF_p, PT\_{HP} / DF_p),
\quad \sigma_p \sim \mathrm{Gamma}(PT\_{HP} \sigma / DF_p, PT\_{HP} /
DF_p),\$\$ \$\$\mu_d \sim \mathrm{Gamma}(\mu, 1), \quad \sigma_d \sim
\mathrm{Gamma}(\sigma, 1),\$\$ with shape and rate parameters, where
\\\mu\\ and \\\sigma\\ are the historical placebo values (`mu_hp` and
`sigma_hp`). The priors for the study drug parameters can be changed
with `mu_d_prior` and `sigma_d_prior`.

Duncan (2026) does not state which relative risk is used for the signal
criterion of this model. This function uses the ratio of the cumulative
incidences at the time `t_star`, \$\$RR(t^\*) = F(t^\*; \mu_d, \sigma_d)
/ F(t^\*; \mu_p, \sigma_p),\$\$ and declares a signal when \\P\[RR(t^\*)
\> rr_0 \mid \mathrm{data}\] \> \\ `cutoff`.

The four-dimensional posterior is integrated numerically on the
coordinates \\\log \mu_p\\, \\\log \sigma_p\\, \\\log \mu_d\\ and the
prior distribution function of \\\sigma_d\\; the last one keeps the
heavy left tail of the prior of \\\sigma_d\\ under control. A coarse
grid first locates the region that holds all but `control$eps` of the
posterior mass in each marginal. On that region, Gauss-Legendre
quadrature is used for three coordinates and Chebyshev interpolation for
\\\log \mu_d\\. Since \\RR(t^\*) \> rr_0\\ holds exactly when \\\mu_d\\
is below a value determined by the other three parameters, the posterior
probability and the quantiles of \\RR(t^\*)\\ are obtained from the
cumulative integral of the interpolant without a discontinuous
integrand. When `control$check` is `TRUE`, the computation is repeated
with about 1.25 times as many Gauss-Legendre nodes and 1.5 times as many
Chebyshev points, the finer result is returned, and twice the absolute
difference in `prob` between the two results is reported as the
integration error estimate.

When there are few events, the prior of \\\sigma_d\\ used by Duncan
(2026) puts much mass near zero, where the likelihood changes sharply in
\\\mu_d\\ near the logarithms of the follow-up times. The integration
error then decreases only in proportion to the inverse of the number of
Chebyshev points, and its sign changes irregularly. With 1.5 times as
many points, the error of the finer result is then about twice the
difference between the two results, which is the reason for the factor 2
in the error estimate. The default number of Chebyshev points is large
for the same reason, and the error estimate is not reliable when
`control$n_inner` is much smaller than the default. If the reported
error estimate is too large for the purpose, increase `control$n_inner`.

The posterior mean of \\RR(t^\*)\\ can be unstable because its posterior
distribution may have a heavy right tail; the posterior median is a more
stable summary.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`ll_cdf`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md),
[`ll_mu`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md),
[`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)

## Examples

``` r
set.seed(1)
n <- 60
t <- runif(n, 4, 104)
y <- rbinom(n, 1, ll_cdf(t, mu = 6.5, sigma = 0.8))
fit <- llbern_post(y, t, k = 1, mu_hp = 6.5, sigma_hp = 0.8, pt_hp = 200,
                   t_star = 104, control = list(check = FALSE))
fit
#> Model: LL_BERN 
#> Relative risk: ratio of cumulative incidences at t_star = 104 
#> P(RR(t_star) > 1 | data) = 0.2741
#> Signal (posterior probability > 0.95): no
#> Posterior of RR(t_star): mean 0.7703, median 0.1454, 95% interval (1.012e-119, 4.212)
#> Cumulative incidence at t_star = 104: placebo 0.09706, study drug 0.06457, difference -0.03249
```
