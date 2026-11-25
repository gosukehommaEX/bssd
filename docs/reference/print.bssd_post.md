# Print a posterior summary for safety signal detection

Summaries that are not available (`NA`), such as the posterior mean of
[`multibin_post`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)
or the quantiles with `details = FALSE`, are not shown.

## Usage

``` r
# S3 method for class 'bssd_post'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `bssd_post`.

- digits:

  Number of significant digits.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
[`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md),
[`multibin_post`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md),
[`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md),
[`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md),
[`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md),
[`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md),
[`mapbern_post`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)

## Examples

``` r
fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                   pt_hp = 4000)
print(fit)
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.8199
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 2.451, median 2.123, 95% interval (0.195, 6.588)
```
