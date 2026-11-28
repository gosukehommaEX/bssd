# Posterior of the event rates under the exposure-time Poisson model

Computes the posterior distribution of the event rates of study drug and
control under the Bayesian exposure-time model of Schnell and Ball
(2016) from blinded subject-level data: the number of events and the
follow-up time of each subject.

## Usage

``` r
exptime_post(
  y,
  t,
  k = 1,
  prior_c,
  prior_d,
  rate0 = NULL,
  rr0 = 1,
  cutoff = 0.95,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Numbers of events of the subjects (non-negative integers).

- t:

  Follow-up times of the subjects (positive), in the time unit of the
  rates.

- k:

  Randomization ratio of study drug to control (`k:1`).

- prior_c:

  Shape and rate `c(alpha_c, beta_c)` of the Gamma prior of the control
  rate.

- prior_d:

  Shape and rate `c(alpha_d, beta_d)` of the Gamma prior of the study
  drug rate.

- rate0:

  Optional critical value of the study drug rate.

- rr0:

  Critical value of the rate ratio, used for the signal when `rate0` is
  `NULL`.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible intervals.

- details:

  Logical; if `FALSE`, the quantiles are not computed (they are `NA`).
  This is faster.

- control:

  List of settings of the numerical integration: `n_coarse` (cells per
  coordinate of the coarse grid, default 41), `n_u` and `n_v` (nodes of
  the two coordinates, defaults 48 and 49), `check` (default `TRUE`),
  `eps` (default 1e-8) and `max_iter` (default 12).

## Value

An object of class `bssd_post`, a list with elements

- `model`:

  `"EXPTIME"`.

- `prob`:

  Posterior probability of the tested quantity (`lambda_d > rate0`, or
  `r > rr0` when `rate0` is `NULL`).

- `signal`:

  Logical; `prob > cutoff`.

- `probs`:

  Both posterior probabilities: `drug` (\\\lambda_d \> rate_0\\, `NA`
  without `rate0`) and `ratio` (\\r \> rr_0\\).

- `rr`:

  Posterior mean, median and the limits of the equal-tailed credible
  interval of the rate ratio \\r\\.

- `rate`:

  The same summaries of \\\lambda_d\\.

- `estimates`:

  Posterior means of `lambda_c`, `lambda_d`, `r` and the difference
  \\\lambda_d - \lambda_c\\.

- `rr_label`, `rr0`, `rr_name`, `rate_name`, `description`, `cutoff`,
  `level`:

  Labels and settings used by `print`.

- `data`, `prior`:

  Data summary and prior parameters.

- `integration`:

  Description of the numerical integration: the first coordinate
  (`coordinate`: `"log"`, `"logit"` or `"normal"`), the numbers of
  nodes, the integration box in the first coordinate and \\z\\, and the
  error estimate.

## Details

For subject \\i\\ with follow-up time \\t_i\\ and \\y_i\\ events, the
model is \$\$Y_i \mid a_i \sim \mathrm{Poisson}(\lambda\_{a_i} t_i),
\quad a_i \sim \mathrm{Bernoulli}(k / (k + 1)),\$\$ \$\$\lambda_c \sim
\mathrm{Gamma}(\alpha_c, \beta_c), \quad \lambda_d \sim
\mathrm{Gamma}(\alpha_d, \beta_d),\$\$ where \\a_i = 1\\ means study
drug and the Gamma distributions are parameterized by shape and rate.
Summing over the unknown treatment assignment of each subject gives the
likelihood \$\$\prod_i \[\lambda_c^{y_i} e^{-\lambda_c t_i} / (k + 1) +
k \lambda_d^{y_i} e^{-\lambda_d t_i} / (k + 1)\],\$\$ so the posterior
is two-dimensional. Schnell and Ball (2016) sample it by Gibbs sampling;
here it is integrated numerically. A coarse grid locates the region that
holds all but `control$eps` of the posterior mass in the coordinates
\\\log \lambda_c\\ and \\z = \Phi^{-1}(F_d(\lambda_d))\\, where \\F_d\\
is the prior distribution function of \\\lambda_d\\. When \\\alpha_c \<
1\\, the prior density of \\\log \lambda_c\\ decreases slowly as
\\\lambda_c \to 0\\, where the blinded data do not bound the posterior,
and \\z_c = \Phi^{-1}(F_c(\lambda_c))\\, where \\F_c\\ is the prior
distribution function of \\\lambda_c\\, is used instead of \\\log
\lambda_c\\. The first coordinate is integrated by Gauss-Legendre
quadrature. For each of its nodes, the interval of \\z\\ that holds the
conditional posterior mass is located separately, and the density is
interpolated at Chebyshev points in \\z\\ or in \\F_d(\lambda_d)\\,
whichever represents it more accurately. This keeps the integration
accurate when the data conflict with an informative prior of
\\\lambda_d\\. When `control$check` is `TRUE`, the computation is
repeated with about 1.5 times as many nodes, the finer result is
returned, and the largest absolute difference in the posterior
probabilities is reported as the integration error estimate.

The prior of the control rate is expected to be informative (for
example, a posterior from historical control data); without it the
blinded data do not separate the two rates. With a vague prior of the
control rate, the posterior can have features that are narrow compared
with its spread, and the integration error estimate shows when more
nodes (`control$n_u`) are needed. When \\\alpha_c \le 1\\, the posterior
density of \\\lambda_c\\ near 0 is proportional to the prior density,
and the posterior mean of the rate ratio is infinite (`Inf`).

When `rate0` is given, a signal is declared when \\P(\lambda_d \> rate_0
\mid \mathrm{data}) \> \\ `cutoff`, which is the rule of Schnell and
Ball (2016); otherwise the rate ratio \\r = \lambda_d / \lambda_c\\ is
tested against `rr0`. Both probabilities are returned in `probs`.

With event indicators as `y` and the times at risk (follow-up up to the
first event) as `t`, the model is the exponential model for the time to
the first event.

## References

Schnell, P. M. and Ball, G. (2016). A Bayesian exposure-time method for
clinical trial safety monitoring with blinded data. *Therapeutic
Innovation & Regulatory Science*, 50(6), 833–838.
[doi:10.1177/2168479016656702](https://doi.org/10.1177/2168479016656702)

## See also

[`elicit_prior`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md)
for priors from quantiles,
[`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
for the exponential Bernoulli mixture model.

## Examples

``` r
set.seed(1)
n <- 300
t <- runif(n, 6, 18)
drug <- rbinom(n, 1, 2 / 3)
y <- rpois(n, ifelse(drug == 1, 0.016, 0.011) * t)
# priors of Schnell and Ball (2016): control rate from 339 events in 26064
# person-months, study drug rate with one hundredth of that information
fit <- exptime_post(y, t, k = 2, prior_c = c(339, 26064),
                    prior_d = c(339 / 100 + 4, 26064 / 100 + 167),
                    rate0 = 0.0156, cutoff = 0.62)
fit
#> Model: EXPTIME 
#> Tested quantity: event rate lambda_d of study drug
#> P(lambda_d > 0.0156 | data) = 0.5387
#> Signal (posterior probability > 0.62): no
#> Posterior of r: mean 1.234, median 1.221, 95% interval (0.8188, 1.721)
#> Posterior of lambda_d: mean 0.01599, median 0.01587, 95% interval (0.01092, 0.02176)
```
