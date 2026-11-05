# bssd: Bayesian Safety Signal Detection for Blinded Clinical Trial Data

The package computes the posterior probability of a safety signal from
blinded clinical trial data with the models examined by Duncan (2026):
the BDRIBS model of Mukhopadhyay et al. (2018) with the total exposure
time or the total time at risk
([`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)),
the exponential Bernoulli mixture model
([`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)),
and the log-logistic Bernoulli mixture model
([`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)).
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
