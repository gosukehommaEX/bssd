# BINBETA and MULTIBIN: blinded binomial models

``` r

library(bssd)
```

## One study: Ball and Schnell (2016)

With \\N\\ enrolled subjects of whom \\Z\\ had the event, the number on
control is \\M_c \sim \mathrm{Binomial}(N, 1 / (k + 1))\\, the numbers
with the event are \\Y_c \sim \mathrm{Binomial}(M_c, \theta_c)\\ and
\\Y_d \sim \mathrm{Binomial}(N - M_c, \theta_d)\\, \\Z = Y_c + Y_d\\,
and the proportions have Beta priors. Summing over \\M_c\\ gives \\Z
\sim \mathrm{Binomial}(N, (\theta_c + k \theta_d) / (k + 1))\\, so the
posterior of \\(\theta_c, \theta_d)\\ is two-dimensional.
[`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)
integrates it numerically; Ball and Schnell (2016) raise a signal when
the posterior probability that \\\theta_d\\ exceeds the maximum
acceptable proportion is larger than a threshold.

``` r

fit <- binbeta_post(y = 40, n = 300, k = 1, prior_c = c(31, 271),
                    prior_d = c(1, 1), p0 = 0.15, cutoff = 0.8)
fit
#> Model: BINBETA 
#> Tested quantity: proportion theta_d of study drug
#> P(theta_d > 0.15 | data) = 0.6615
#> Signal (posterior probability > 0.8): no
#> Posterior of rr: mean 1.724, median 1.652, 95% interval (0.7135, 3.151)
#> Posterior of theta_d: mean 0.1689, median 0.1677, 95% interval (0.08793, 0.2566)
```

Before the trial, the largest acceptable number of subjects with the
event can be computed for each number of enrolled subjects (Figure 1 of
Ball and Schnell, 2016);
[`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
returns the smallest number with a signal, which is one more than the
largest acceptable number.

``` r

signal_count(n = c(50, 100, 200, 300), method = "binbeta",
             prior_c = c(31, 271), prior_d = c(1, 1), p0 = 0.15,
             cutoff = 0.8)
#>     n y_crit      prob prob_below
#> 1  50      8 0.8091985  0.6895111
#> 2 100     16 0.8599635  0.7912753
#> 3 200     30 0.8455281  0.7965627
#> 4 300     43 0.8087599  0.7649654
```

The binomial expansion of the likelihood makes the posterior a finite
mixture of products of Beta distributions, which gives exact values. The
largest absolute difference of the 21 compared quantities from these
exact values (`inst/validation/validate_methods_python.R`, relative
tolerance 1e-5) is 8.7e-07.

## Several studies: Goren et al. (2020)

For several blinded studies, Goren et al. (2020) assume \\Y_s \sim
\mathrm{Binomial}(N_s, \theta_s (1 - r_s + r_s \delta))\\ with
study-specific control probabilities \\\theta_s \sim \mathrm{Beta}(a_s,
b_s)\\, a common rate ratio with \\\log \delta \sim N(m, v^2)\\, and
\\\theta_s \le 1 / \delta\\. Given \\\delta\\, each \\\theta_s\\ is
integrated out separately, so the posterior of \\\delta\\ is
one-dimensional and
[`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)
integrates it without MCMC.

The real data of Goren et al. (2020, Table 3) are eight studies; the
homogeneous control prior is from 20 events among 1687 control patients
of a previous trial.

``` r

y <- c(8, 4, 5, 6, 20, 2, 7, 39)
n <- c(354, 305, 307, 582, 244, 68, 307, 1612)
k <- c(2, 2, 1, 1, 2, 1, 2, 1)
fit <- multibin_post(y, n, k = k, prior_c = c(20, 1667), rr0 = 2)
fit
#> Model: MULTIBIN 
#> Tested quantity: rate ratio delta common to the studies
#> P(delta > 2 | data) = 0.9828
#> Signal (posterior probability > 0.95): yes
#> Posterior of delta: mean NA, median 2.951, 95% interval (2.061, 4.093)
```

`inst/reproduce/reproduce_Goren2020.R` compares the posterior
probabilities of Table 4 (homogeneous and heterogeneous priors, all
studies and without study 5), which the authors computed by Hamiltonian
Monte Carlo:

| quantity | published | bssd | tolerance | judgement |
|:---|---:|---:|---:|:---|
| P(delta \> 1), full analysis, homogeneous priors | 1.000 | 1.000 | 0.0020 | PASS |
| P(delta \> 2), full analysis, homogeneous priors | 0.984 | 0.983 | 0.0035 | PASS |
| P(delta \> 3), full analysis, homogeneous priors | 0.462 | 0.462 | 0.0125 | PASS |
| P(delta \> 1), full analysis, heterogeneous priors | 1.000 | 1.000 | 0.0020 | PASS |
| P(delta \> 2), full analysis, heterogeneous priors | 0.836 | 0.839 | 0.0095 | PASS |
| P(delta \> 3), full analysis, heterogeneous priors | 0.158 | 0.158 | 0.0095 | PASS |
| P(delta \> 1), sensitivity analysis, homogeneous priors | 0.998 | 0.998 | 0.0020 | PASS |
| P(delta \> 2), sensitivity analysis, homogeneous priors | 0.561 | 0.563 | 0.0125 | PASS |
| P(delta \> 3), sensitivity analysis, homogeneous priors | 0.034 | 0.035 | 0.0035 | PASS |
| P(delta \> 1), sensitivity analysis, heterogeneous priors | 0.994 | 0.995 | 0.0035 | PASS |
| P(delta \> 2), sensitivity analysis, heterogeneous priors | 0.486 | 0.485 | 0.0125 | PASS |
| P(delta \> 3), sensitivity analysis, heterogeneous priors | 0.031 | 0.031 | 0.0035 | PASS |

The simulation of Table 2
(`inst/reproduce/reproduce_Goren2020_simulation.R`) has not been run for
this installation.

Against nested adaptive quadrature in Python, the largest absolute
difference of the 10 compared quantities is 3.7e-07.

## Comparison with MCMC

| model | quantity | MCMC | bssd | judgement |
|:---|:---|---:|---:|:---|
| BINBETA | P(theta_d \> 0.15 \| data) | 0.6608 | 0.6615 | PASS |
| BINBETA | P(rr \> 1 \| data) | 0.8965 | 0.8969 | PASS |
| BINBETA | MCMC proportion below the bssd posterior median | 0.5000 | 0.5015 | PASS |
| MULTIBIN | P(delta \> 1 \| data) | 0.7656 | 0.7692 | PASS |
| MULTIBIN | posterior mean of log(delta) | 0.3546 | 0.3694 | PASS |
| MULTIBIN | MCMC proportion below the bssd posterior median | 0.5000 | 0.4994 | PASS |

The JAGS models are in `inst/validation/validate_jags_methods.R`; for
BINBETA, the number on control and the numbers with the event in each
arm are latent variables.

## References

Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
the FDA IND reporting final rule. In J. Lin et al. (eds.), *Statistical
Applications from Clinical Trials and Personalized Medicine to Finance
and Business Analytics*, ICSA Book Series in Statistics, 201–209.
Springer, Cham. <https://doi.org/10.1007/978-3-319-42568-9_17>

Goren, E., Lin, L.-A., Blaustein, R. O. and Ball, G. (2020). Bayesian
meta-analysis of safety outcomes using blinded clinical trial data.
*Therapeutic Innovation & Regulatory Science*, 54, 1557–1565.
<https://doi.org/10.1007/s43441-020-00183-7>
