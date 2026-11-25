# CANARY: the canary in the coal mine

``` r

library(bssd)
```

## Model

Gould and Wang (2017) use external information about the event rate in
patients not receiving the test product to judge whether the blinded
total event count \\X\\ is higher than expected. With \\N\\ subjects,
\\N_d = \mathrm{round}(k N / (k + 1))\\ on study drug and \\N_c = N -
N_d\\ on control (the expected arm sizes), \\X\\ is the sum of
\\\mathrm{Binomial}(N_d, p_d)\\ and \\\mathrm{Binomial}(N_c, p_c)\\
counts. The control probability has an informative Beta prior, and the
metric of the difference (the ratio \\R = p_d / p_c\\, the difference,
or the odds ratio) has a diffuse prior: for the ratio, a lognormal
distribution with mean 1 and probability 0.95 between 0.1 and 10,
truncated to \\R \le 1 / p_c\\. For large trials and rare events, a
Poisson likelihood with a Gamma prior of the control rate is used
instead.
[`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
integrates the two-dimensional posterior (one-dimensional for the
Poisson likelihood) numerically.

``` r

canary_post(y = 12, n = 100, k = 1, prior_c = c(1, 49))
#> Model: CANARY 
#> Tested quantity: ratio R = p_d / p_c of the event probabilities
#> P(R > 1 | data) = 0.8417
#> Signal (posterior probability > 0.95): no
#> Posterior of R: mean 3.156, median 2.389, 95% interval (0.3521, 10.53)
canary_post(y = 12, n = 100, k = 1, prior_c = c(1, 49),
            metric = "difference")
#> Model: CANARY 
#> Tested quantity: difference D = p_d - p_c of the event probabilities
#> P(D > 0 | data) = 0.9956
#> Signal (posterior probability > 0.95): yes
#> Posterior of D: mean 0.2108, median 0.2113, 95% interval (0.06536, 0.3504)
```

## Decision rules and their operating characteristics

Over the range of counts that matters for monitoring, the posterior
probability increases with the event count, so a rule “\\P(R \>
R\_\mathrm{crit} \mid X) \> \gamma\_\mathrm{crit}\\” is equivalent to
\\X \ge X\_\mathrm{crit}\\. (For counts close to \\N\\, both event
probabilities must be close to 1, and the probability of a large ratio
falls again;
[`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
therefore searches upward from 0.)
[`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
finds \\X\_\mathrm{crit}\\ for each sample size (the control charts of
Figure 2), and
[`count_oc()`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)
gives the exact probability of reaching it under true event rates
(Figure 3 and Table V).

``` r

sc <- signal_count(n = c(100, 200, 500), method = "canary",
                   prior_c = c(20, 980), crit = 2, cutoff = 0.8)
sc
#>     n y_crit      prob prob_below
#> 1 100      8 0.8516000  0.7357937
#> 2 200     12 0.8264072  0.7406849
#> 3 500     24 0.8159702  0.7685790
count_oc(sc$y_crit, sc$n, rate_c = 0.02, rate_d = 0.08)
#> [1] 0.1259040 0.2990768 0.6115370
```

## Reproduction of Gould and Wang (2017)

`inst/reproduce/reproduce_Gould2017.R` compares

- the 24 prior parameters of Table I, from
  [`elicit_prior()`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md):
  24 PASS;
- the 24 critical counts of Table III (binomial likelihood): 11 PASS and
  13 EXPLAINED;
- the 4 critical counts of Table V (Poisson likelihood): 0 PASS and 4
  EXPLAINED;
- the 36 probabilities of Table V, from
  [`count_oc()`](https://gosukehommaex.github.io/bssd/reference/count_oc.md):
  36 PASS.

The published critical counts were computed by the authors’ R code on
fixed grids. A transcription of the code, kept outside the package
because the code has no stated license, reproduces 28 of the 28
published counts
(`inst/validation/python/ref_gould2017_author_grid.csv`). The authors’
code itself, run by a script outside the package, gives 28 of the 28
published counts (`inst/validation/results/gould2017_author_code.csv`).
The differences from the exact counts of bssd come from the grids (100
values of \\R\\ in \\(0, 10)\\, with the threshold 2 evaluated at the
grid point 1.98, for the binomial model; a coarse grid of the expected
control count for the Poisson model). The binomial counts of Table III
were computed with the lognormal prior between 0.1 and 8, the default of
the authors’ code, which is used here as well.

| quantity | published | bssd | judgement |
|:---|---:|---:|:---|
| prior (1, 49), N = 100, R_crit = 1, gamma = 0.8 | 12 | 12 | PASS |
| prior (1, 49), N = 100, R_crit = 2, gamma = 0.5 | 11 | 11 | PASS |
| prior (1, 49), N = 100, R_crit = 2, gamma = 0.8 | 18 | 17 | EXPLAINED |
| prior (1, 49), N = 1000, R_crit = 1, gamma = 0.8 | 90 | 87 | EXPLAINED |
| prior (1, 49), N = 1000, R_crit = 2, gamma = 0.5 | 90 | 87 | EXPLAINED |
| prior (1, 49), N = 1000, R_crit = 2, gamma = 0.8 | 144 | 139 | EXPLAINED |
| prior (20, 980), N = 100, R_crit = 1, gamma = 0.8 | 6 | 6 | PASS |
| prior (20, 980), N = 100, R_crit = 2, gamma = 0.5 | 6 | 6 | PASS |
| prior (20, 980), N = 100, R_crit = 2, gamma = 0.8 | 8 | 8 | PASS |
| prior (20, 980), N = 1000, R_crit = 1, gamma = 0.8 | 31 | 30 | EXPLAINED |
| prior (20, 980), N = 1000, R_crit = 2, gamma = 0.5 | 36 | 36 | PASS |
| prior (20, 980), N = 1000, R_crit = 2, gamma = 0.8 | 44 | 44 | PASS |
| prior (5, 45), N = 100, R_crit = 1, gamma = 0.8 | 21 | 21 | PASS |
| prior (5, 45), N = 100, R_crit = 2, gamma = 0.5 | 23 | 23 | PASS |
| prior (5, 45), N = 100, R_crit = 2, gamma = 0.8 | 31 | 30 | EXPLAINED |
| prior (5, 45), N = 1000, R_crit = 1, gamma = 0.8 | 182 | 178 | EXPLAINED |
| prior (5, 45), N = 1000, R_crit = 2, gamma = 0.5 | 206 | 204 | EXPLAINED |
| prior (5, 45), N = 1000, R_crit = 2, gamma = 0.8 | 275 | 272 | EXPLAINED |
| prior (100, 900), N = 100, R_crit = 1, gamma = 0.8 | 16 | 16 | PASS |
| prior (100, 900), N = 100, R_crit = 2, gamma = 0.5 | 18 | 18 | PASS |
| prior (100, 900), N = 100, R_crit = 2, gamma = 0.8 | 22 | 21 | EXPLAINED |
| prior (100, 900), N = 1000, R_crit = 1, gamma = 0.8 | 118 | 116 | EXPLAINED |
| prior (100, 900), N = 1000, R_crit = 2, gamma = 0.5 | 157 | 155 | EXPLAINED |
| prior (100, 900), N = 1000, R_crit = 2, gamma = 0.8 | 173 | 171 | EXPLAINED |

## Numerical integration

For the binomial likelihood, the posterior is integrated over
\\\mathrm{logit}(p_c)\\ by Gauss-Legendre quadrature, and for each of
its nodes over an interval of \\\log R\\ (or \\p_d\\) located from the
conditional posterior; for the ratio the interval ends at the truncation
point \\R = 1 / p_c\\, where the density of \\\log R\\ can be largest
when many subjects had the event. Against independent Python
computations on dense grids split exactly at the thresholds, the largest
absolute difference of the 27 compared quantities (five data sets: the
ratio, the difference, the odds ratio, the Poisson likelihood, and the
ratio when all 100 subjects had the event) is 5.8e-07. For the Poisson
likelihood, MCMC in JAGS gives 0.2214 for \\P(R \> 2 \mid X)\\ against
0.2213 from bssd (PASS).

## References

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92–104.
<https://doi.org/10.1002/sim.7129>
