# Posterior of the risk elevation from the total event count

Computes the posterior distribution of a metric of the difference in
event risk between study drug and control (rate ratio, risk difference
or odds ratio) from the blinded total event count, with the method of
Gould and Wang (2017), which uses external information about the control
event rate.

## Usage

``` r
canary_post(
  y,
  n,
  k = 1,
  prior_c,
  metric = c("ratio", "difference", "odds"),
  metric_prior = NULL,
  likelihood = c("binomial", "poisson"),
  crit = NULL,
  cutoff = 0.95,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Total number of subjects with the event (binomial) or of events
  (Poisson) in the blinded data.

- n:

  Number of subjects (binomial) or total exposure (Poisson).

- k:

  Randomization ratio of study drug to control (`k:1`).

- prior_c:

  Parameters of the prior of the control event probability (`c(a, b)` of
  a Beta distribution, binomial) or rate (shape and rate of a Gamma
  distribution, Poisson).

- metric:

  One of `"ratio"`, `"difference"` and `"odds"`; only `"ratio"` is
  available for the Poisson likelihood.

- metric_prior:

  `c(meanlog, sdlog)` of the lognormal prior of the ratio or odds ratio,
  or `NULL` for the default (see Details). Ignored for the difference.

- likelihood:

  `"binomial"` or `"poisson"`.

- crit:

  Critical value of the metric; by default 1 for the ratio and the odds
  ratio and 0 for the difference.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible interval of the metric.

- details:

  Logical; if `FALSE`, the quantiles are not computed (they are `NA`).

- control:

  List of settings of the numerical integration, as in
  [`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md)
  (`n_u` is not used for the Poisson likelihood).

## Value

An object of class `bssd_post` with model `"CANARY"`, `prob` (posterior
probability that the metric exceeds `crit`), `signal`, `rr` (posterior
mean, median and credible interval of the metric), `estimates`
(posterior means of the control and study drug probabilities or rates
and of the metric), labels and settings, `data`, `prior` and
`integration`.

## Details

**Binomial likelihood.** Of \\N\\ subjects, \\N_d = \mathrm{round}(k N /
(k + 1))\\ are on study drug and \\N_c = N - N_d\\ on control (the
expected arm sizes, as in Gould and Wang, 2017). The total count \\X\\
is the sum of independent \\\mathrm{Binomial}(N_d, p_d)\\ and
\\\mathrm{Binomial}(N_c, p_c)\\ counts, and its probability is the
convolution of the two. The control probability has the prior \\p_c \sim
\mathrm{Beta}(a, b)\\. The metric is

- `"ratio"`: \\R = p_d / p_c\\ with a lognormal prior truncated to \\R
  \le 1 / p_c\\, as in the code of Gould and Wang (2017);

- `"odds"`: the odds ratio \\Q\\ with a lognormal prior;

- `"difference"`: \\D = p_d - p_c\\ with a uniform prior on \\(-p_c, 1 -
  p_c)\\, which is a uniform prior on \\p_d\\ independent of \\p_c\\.

The two-dimensional posterior is integrated on \\\mathrm{logit}(p_c)\\
(Gauss-Legendre quadrature) and the logarithm of the metric, or \\p_d\\
for the difference (Chebyshev interpolation).

**Poisson likelihood.** With the total exposure \\E\\ (`n`) and the
control event rate \\\lambda_c\\ per unit of exposure, \\X \sim
\mathrm{Poisson}(E \lambda_c (1 + k R) / (k + 1))\\ with \\\lambda_c
\sim \mathrm{Gamma}(a, b)\\ (shape and rate) and a lognormal prior on
\\R\\. The control rate is integrated out analytically and the posterior
of \\\log R\\ is integrated numerically. Gould and Wang (2017) write the
model with the expected number of control events; their Gamma prior with
rate \\b'\\ on that number corresponds to \\b = b' E / (k + 1)\\. The
likelihood is that of the BDRIBS model
([`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md));
only the prior of the rate ratio differs.

The default `metric_prior` is the lognormal distribution with mean 1 and
probability 0.95 between 0.1 and 10 (Gould and Wang, 2017, Section 3.1);
`elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)` gives its
parameters. A signal is declared when \\P(M \> crit \mid X) \> \\
`cutoff`.

## References

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92–104.
[doi:10.1002/sim.7129](https://doi.org/10.1002/sim.7129)

## See also

[`signal_count`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
and
[`count_oc`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)
for the critical event counts and their operating characteristics.

## Examples

``` r
# 12 subjects with the event among 100, control prior Beta(1, 49)
canary_post(y = 12, n = 100, prior_c = c(1, 49))
#> Model: CANARY 
#> Tested quantity: ratio R = p_d / p_c of the event probabilities
#> P(R > 1 | data) = 0.8417
#> Signal (posterior probability > 0.95): no
#> Posterior of R: mean 3.156, median 2.389, 95% interval (0.3521, 10.53)
# the same data with the risk difference
canary_post(y = 12, n = 100, prior_c = c(1, 49), metric = "difference")
#> Model: CANARY 
#> Tested quantity: difference D = p_d - p_c of the event probabilities
#> P(D > 0 | data) = 0.9956
#> Signal (posterior probability > 0.95): yes
#> Posterior of D: mean 0.2108, median 0.2113, 95% interval (0.06536, 0.3504)
```
