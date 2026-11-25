# Probability that the blinded event count reaches a critical value

Computes the exact probability that the blinded total event count is at
least `y_crit` under true event probabilities or rates of study drug and
control: the operating characteristics of a rule based on
[`signal_count`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
(Gould and Wang, 2017, Figure 3 and Table V).

## Usage

``` r
count_oc(
  y_crit,
  n,
  rate_c,
  rate_d,
  k = 1,
  likelihood = c("binomial", "poisson"),
  allocation = c("fixed", "random")
)
```

## Arguments

- y_crit:

  Critical counts (`NA` gives probability 0).

- n:

  Numbers of subjects, or total exposures for the Poisson likelihood.

- rate_c, rate_d:

  True event probabilities (binomial) or rates (Poisson) of control and
  study drug.

- k:

  Randomization ratio of study drug to control (`k:1`).

- likelihood:

  `"binomial"` or `"poisson"`.

- allocation:

  `"fixed"` or `"random"` arm sizes (binomial likelihood only).

## Value

Numeric vector of probabilities; the arguments `y_crit`, `n`, `rate_c`
and `rate_d` are recycled to a common length.

## Details

With `likelihood = "binomial"` and `allocation = "fixed"`, the count is
the sum of \\\mathrm{Binomial}(N_d, p_d)\\ and \\\mathrm{Binomial}(N_c,
p_c)\\ counts with \\N_d = \mathrm{round}(k N / (k + 1))\\ (the model of
[`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md));
with `allocation = "random"`, it is \\\mathrm{Binomial}(N, (p_c + k p_d)
/ (k + 1))\\ (the model of
[`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)).
With `likelihood = "poisson"`, `n` is the total exposure and the count
is \\\mathrm{Poisson}(n (\lambda_c + k \lambda_d) / (k + 1))\\.

## References

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92–104.
[doi:10.1002/sim.7129](https://doi.org/10.1002/sim.7129)

## See also

[`signal_count`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)

## Examples

``` r
# probability of 6 or more events among 100 subjects when the true
# probabilities are 0.02 (control) and 0.02, 0.04 or 0.08 (study drug)
count_oc(6, 100, rate_c = 0.02, rate_d = c(0.02, 0.04, 0.08))
#> [1] 0.01548364 0.08048695 0.38392900
```
