# Smallest event count that gives a signal

For each number of subjects (or exposure), finds the smallest blinded
event count for which the posterior probability exceeds `cutoff`.
Plotted against the number of subjects, these counts are the boundaries
of Ball and Schnell (2016, Figure 1) and the control charts of Gould and
Wang (2017, Figure 2 and Table III).

## Usage

``` r
signal_count(
  n,
  method = c("canary", "binbeta"),
  ...,
  cutoff = 0.95,
  y_max = NULL
)
```

## Arguments

- n:

  Numbers of subjects (or exposures for the Poisson likelihood of
  [`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)).

- method:

  `"canary"`
  ([`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md))
  or `"binbeta"`
  ([`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)).

- ...:

  Further arguments of the posterior function (priors, `k`, critical
  values); `y`, `n`, `cutoff` and `details` are set here.

- cutoff:

  Posterior probability threshold for a signal.

- y_max:

  Largest count to consider; by default `n` (rounded down). It must be
  given for the Poisson likelihood.

## Value

A data frame with columns `n`, `y_crit` (the smallest count with a
signal, `NA` if none of the evaluated counts up to `y_max` gives one),
`prob` (the posterior probability at `y_crit`, or at `y_max` if `y_crit`
is `NA`) and `prob_below` (at `y_crit - 1`, `NA` if `y_crit` is 0).

## Details

The posterior probability increases with the event count up to the first
count with a signal, but not necessarily beyond it: when the event
probabilities are bounded (binomial likelihood), counts near `n` imply
event probabilities near 1 in both arms, and the probability of a large
ratio falls again. The counts 1, 2, 4, 8, ... (up to `y_max`) are
therefore evaluated until the first one with a signal, and the smallest
count with a signal is found by bisection between that count and the one
before it. A signal at count \\y\\ means `prob > cutoff`, the rule of
all functions of the package; Gould and Wang (2017) use \\\ge\\, which
differs only when the probability equals `cutoff` exactly. The largest
acceptable count of Ball and Schnell (2016) is the result minus 1.

## References

Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
the FDA IND reporting final rule. In J. Lin et al. (eds.), *Statistical
Applications from Clinical Trials and Personalized Medicine to Finance
and Business Analytics*, ICSA Book Series in Statistics, 201–209.
Springer, Cham.
[doi:10.1007/978-3-319-42568-9_17](https://doi.org/10.1007/978-3-319-42568-9_17)

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92–104.
[doi:10.1002/sim.7129](https://doi.org/10.1002/sim.7129)

## See also

[`count_oc`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)

## Examples

``` r
# Gould and Wang (2017), Table III: control prior Beta(20, 980), signal if
# P(R > 1) > 0.8
signal_count(n = 100, method = "canary", prior_c = c(20, 980),
             cutoff = 0.8)
#>     n y_crit      prob prob_below
#> 1 100      6 0.8421326  0.7167684
```
