# Location parameter of the log-logistic distribution for a given incidence

Returns the location parameter \\\mu\\ for which the log-logistic
cumulative incidence (see
[`ll_cdf`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)) at
time `t` equals `ip`: \$\$\mu = \log t - \sigma \log\[ip / (1 -
ip)\].\$\$

## Usage

``` r
ll_mu(ip, t, sigma)
```

## Arguments

- ip:

  Cumulative incidence at time `t`, in `(0, 1)`.

- t:

  Time (positive).

- sigma:

  Scale parameter (positive).

## Value

Numeric vector of location parameters. The arguments are recycled to a
common length.

## See also

[`ll_cdf`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)

## Examples

``` r
# 3.5% incidence at week 104 with sigma = 0.8 gives mu close to 7.3
ll_mu(0.035, t = 104, sigma = 0.8)
#> [1] 7.297815
```
