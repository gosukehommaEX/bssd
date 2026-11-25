# Posterior of the rate ratio from blinded event counts of several studies

Computes the posterior distribution of a rate ratio common to several
blinded studies from the total numbers of subjects with the event in
each study, with the Bayesian model of Goren et al. (2020).

## Usage

``` r
multibin_post(
  y,
  n,
  k = 1,
  prior_c,
  lrr_prior = c(0, 5),
  rr0 = 1,
  cutoff = 0.95,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Numbers of subjects with the event, one per study.

- n:

  Numbers of enrolled subjects, one per study.

- k:

  Randomization ratios of study drug to control, one per study or one
  for all.

- prior_c:

  Shape parameters of the Beta prior of the control probabilities:
  `c(a, b)` for all studies, or a matrix with one row `c(a_s, b_s)` per
  study.

- lrr_prior:

  Mean and standard deviation `c(m, v)` of the normal prior of \\\log
  \delta\\; the default is that of Goren et al. (2020).

- rr0:

  Critical value of the rate ratio.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible interval.

- details:

  Logical; if `FALSE`, the quantiles are not computed (they are `NA`).

- control:

  List of settings of the numerical integration: `n_coarse` (default
  41), `n_v` (Chebyshev points, default 49), `n_theta` (Gauss-Legendre
  nodes per panel of the integral over each control probability, default
  12), `check` (default `TRUE`), `eps` (default 1e-8) and `max_iter`
  (default 12).

## Value

An object of class `bssd_post` with model `"MULTIBIN"`, `prob`,
`signal`, `rr` (posterior median and credible interval of \\\delta\\;
the mean is `NA`), `estimates` (posterior mean of \\\log \delta\\),
labels and settings, `data`, `prior` and `integration`.

## Details

Study \\s\\ has \\N_s\\ enrolled subjects, \\Y_s\\ subjects with the
event, and the probability \\r_s = k_s / (k_s + 1)\\ of study drug. The
model is \$\$Y_s \sim \mathrm{Binomial}(N_s, \theta_s (1 - r_s + r_s
\delta)), \quad \theta_s \sim \mathrm{Beta}(a_s, b_s), \quad \log \delta
\sim N(m, v^2),\$\$ restricted to \\\theta_s \le 1 / \delta\\, so that
the study drug probability \\\delta \theta_s\\ is at most 1. The control
probabilities \\\theta_s\\ may differ between studies, and \\\delta\\ is
common. Given \\\delta\\, each \\\theta_s\\ is integrated out by
composite Gauss-Legendre quadrature on \\\mathrm{logit}(\theta_s)\\ over
the region where the integrand is not negligible, so the posterior of
\\\delta\\ is one-dimensional; it is integrated at Chebyshev points of
the prior distribution function of \\\log \delta\\. Goren et al. (2020)
sample the same posterior by Hamiltonian Monte Carlo.

A signal is declared when \\P(\delta \> rr_0 \mid \mathrm{data}) \> \\
`cutoff`. With a single study, the likelihood is that of
[`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)
with \\\theta_d = \delta \theta_c\\.

The posterior mean of \\\delta\\ is not reported, because with the vague
prior recommended by Goren et al. (2020) (\\v = 5\\) it depends on the
far right tail of the prior; the posterior median and the credible
interval are reported instead, and `estimates` holds the posterior mean
of \\\log \delta\\.

## References

Goren, E., Lin, L.-A., Blaustein, R. O. and Ball, G. (2020). Bayesian
meta-analysis of safety outcomes using blinded clinical trial data.
*Therapeutic Innovation & Regulatory Science*, 54, 1557–1565.
[doi:10.1007/s43441-020-00183-7](https://doi.org/10.1007/s43441-020-00183-7)

## Examples

``` r
# four studies with randomization 1:1 and a control prior with mean 0.02
multibin_post(y = c(6, 10, 25, 47), n = c(200, 400, 800, 1600), k = 1,
              prior_c = c(2, 98))
#> Model: MULTIBIN 
#> Tested quantity: rate ratio delta common to the studies
#> P(delta > 1 | data) = 0.7692
#> Signal (posterior probability > 0.95): no
#> Posterior of delta: mean NA, median 1.675, 95% interval (0.1729, 5.239)
```
