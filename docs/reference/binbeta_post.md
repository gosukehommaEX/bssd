# Posterior of the event proportions under the blinded binomial model

Computes the posterior distribution of the proportions of subjects with
the event on study drug and control under the Bayesian hierarchical
model of Ball and Schnell (2016), from the blinded number of enrolled
subjects and the number of subjects with the event.

## Usage

``` r
binbeta_post(
  y,
  n,
  k = 1,
  prior_c,
  prior_d,
  p0 = NULL,
  rr0 = 1,
  cutoff = 0.95,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Number of subjects with the event (non-negative integer).

- n:

  Number of enrolled subjects (positive integer, at least `y`).

- k:

  Randomization ratio of study drug to control (`k:1`).

- prior_c:

  Shape parameters `c(a, b)` of the Beta prior of the control
  proportion.

- prior_d:

  Shape parameters `c(c, d)` of the Beta prior of the study drug
  proportion.

- p0:

  Optional critical value of the study drug proportion (the maximum
  acceptable proportion).

- rr0:

  Critical value of the ratio of the proportions, used for the signal
  when `p0` is `NULL`.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible intervals.

- details:

  Logical; if `FALSE`, the quantiles are not computed (they are `NA`).

- control:

  List of settings of the numerical integration, as in
  [`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md).

## Value

An object of class `bssd_post` as described in
[`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md),
with model `"BINBETA"`, the control and study drug proportions `theta_c`
and `theta_d`, and their ratio `rr` in place of the rates and the rate
ratio.

## Details

With \\N\\ enrolled subjects of whom \\Z\\ had the event, the model is
\$\$M_c \sim \mathrm{Binomial}(N, 1 / (k + 1)), \quad Y_c \mid M_c \sim
\mathrm{Binomial}(M_c, \theta_c), \quad Y_d \mid M_c \sim
\mathrm{Binomial}(N - M_c, \theta_d),\$\$ \$\$Z = Y_c + Y_d, \quad
\theta_c \sim \mathrm{Beta}(a, b), \quad \theta_d \sim \mathrm{Beta}(c,
d).\$\$ Under simple randomization each subject is on control with
probability \\1 / (k + 1)\\, so \\Z \sim \mathrm{Binomial}(N, \pi)\\
with \\\pi = (\theta_c + k \theta_d) / (k + 1)\\, and the posterior is
two-dimensional. It is integrated numerically on the coordinates
\\\mathrm{logit}(\theta_c)\\ and the prior distribution function of
\\\theta_d\\ as in
[`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md).
When \\a \< 1\\ or \\b \< 1\\, the prior density of
\\\mathrm{logit}(\theta_c)\\ decreases slowly in a tail, which the
blinded data do not bound, and \\z_c = \Phi^{-1}(F_c(\theta_c))\\, where
\\F_c\\ is the prior distribution function of \\\theta_c\\, is used
instead. For \\rr_0 \> 1\\, the Gauss-Legendre rule over the first
coordinate is split at \\\theta_c = 1 / rr_0\\, where the posterior
probability of the ratio has a kink because \\\theta_d \le 1\\. When \\a
\le 1\\, the posterior density of \\\theta_c\\ near 0 is proportional to
the prior density, and the posterior mean of the ratio is infinite
(`Inf`).

When `p0` is given, a signal is declared when \\P(\theta_d \> p_0 \mid
Z) \> \\ `cutoff` (the rule of Ball and Schnell, 2016); otherwise the
ratio \\\theta_d / \theta_c\\ is tested against `rr0`. With a single
study, the likelihood is that of the model of Goren et al. (2020), see
[`multibin_post`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md).

## References

Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
the FDA IND reporting final rule. In J. Lin et al. (eds.), *Statistical
Applications from Clinical Trials and Personalized Medicine to Finance
and Business Analytics*, ICSA Book Series in Statistics, 201–209.
Springer, Cham.
[doi:10.1007/978-3-319-42568-9_17](https://doi.org/10.1007/978-3-319-42568-9_17)

## See also

[`signal_count`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
for the largest acceptable number of subjects with the event,
[`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
for the model with fixed arm sizes of Gould and Wang (2017).

## Examples

``` r
# 40 of 300 subjects with the event, randomization 1:1, control prior from
# 30 of 300 historical control subjects, vague prior for study drug
fit <- binbeta_post(y = 40, n = 300, k = 1, prior_c = c(31, 271),
                    prior_d = c(1, 1), p0 = 0.15)
fit
#> Model: BINBETA 
#> Tested quantity: proportion theta_d of study drug
#> P(theta_d > 0.15 | data) = 0.6615
#> Signal (posterior probability > 0.95): no
#> Posterior of rr: mean 1.724, median 1.652, 95% interval (0.7135, 3.151)
#> Posterior of theta_d: mean 0.1689, median 0.1677, 95% interval (0.08793, 0.2566)
```
