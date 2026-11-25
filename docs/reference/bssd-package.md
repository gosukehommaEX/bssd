# bssd: Bayesian Safety Signal Detection for Blinded Clinical Trial Data

The package computes the posterior probability of a safety signal from
blinded clinical trial data. It implements the models examined by Duncan
(2026): the BDRIBS model of Mukhopadhyay et al. (2018) with the total
exposure time or the total time at risk
([`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)),
the exponential Bernoulli mixture model
([`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)),
and the log-logistic Bernoulli mixture model
([`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)).
It also implements the exposure-time Poisson model of Schnell and Ball
(2016)
([`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md)),
the binomial model of Ball and Schnell (2016)
([`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)),
the method of Gould and Wang (2017)
([`canary_post`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)),
the meta-analytic predictive prior and mixture likelihood of Brock et
al. (2023)
([`map_beta`](https://gosukehommaex.github.io/bssd/reference/map_beta.md),
[`mapbern_post`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)),
the model of Goren et al. (2020) for several studies
([`multibin_post`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)),
and the tools of Waterhouse et al. (2022) for BDRIBS
([`pooled_ratio`](https://gosukehommaex.github.io/bssd/reference/pooled_ratio.md),
[`bdribs_grid`](https://gosukehommaex.github.io/bssd/reference/bdribs_grid.md)).
All posterior quantities are obtained by numerical integration, so no
Markov chain Monte Carlo sampling is needed.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
detection of potential risk using inference on blinded safety data.
*Pharmaceutical Statistics*, 17(6), 823–834.
[doi:10.1002/pst.1898](https://doi.org/10.1002/pst.1898)

## See also

Useful links:

- <https://github.com/gosukehommaEX/bssd>

- <https://gosukehommaex.github.io/bssd/>

- Report bugs at <https://github.com/gosukehommaEX/bssd/issues>

## Author

**Maintainer**: Gosuke Homma <gosuke.homma@astellas.com>
([ORCID](https://orcid.org/0000-0002-6854-5627))

Authors:

- Gosuke Homma <gosuke.homma@astellas.com>
  ([ORCID](https://orcid.org/0000-0002-6854-5627))
