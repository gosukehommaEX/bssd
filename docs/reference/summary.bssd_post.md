# Summarize a posterior for safety signal detection as a data frame

Summarize a posterior for safety signal detection as a data frame

## Usage

``` r
# S3 method for class 'bssd_post'
summary(object, ...)
```

## Arguments

- object:

  An object of class `bssd_post`.

- ...:

  Not used.

## Value

A data frame with one row and the columns `model`, `prob`, `signal`,
`rr_mean`, `rr_median`, `rr_lower`, `rr_upper`, `ip_placebo`, `ip_drug`,
`ip_difference` and `error`. The incidence columns are `NA` when
`t_star` was not given.

## See also

[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
[`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md),
[`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)

## Examples

``` r
fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                   pt_hp = 4000, t_star = 1)
summary(fit)
#>    model     prob signal  rr_mean rr_median  rr_lower rr_upper  ip_placebo
#> 1 BDRIBS 0.819888  FALSE 2.450603  2.123375 0.1950283 6.588037 0.004655827
#>      ip_drug ip_difference error
#> 1 0.01058191   0.005926085    NA
```
