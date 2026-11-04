# Validation

The scripts in this folder compare the functions of bssd with independent
computations. Run them from the package root after `devtools::load_all()` (or
`library(bssd)`), for example `source("inst/validation/validate_jags.R")`.

| Script | Comparison | Requirements | Run time |
|---|---|---|---|
| `validate_bdribs_package.R` | `bdribs_post()` against the bdribs package (MCMC by the authors of the BDRIBS method), with a known and with a Gamma-distributed background rate | JAGS 4.x, rjags, coda, bdribs | a few minutes |
| `validate_jags.R` | `expbern_post()` and `llbern_post()` against MCMC with JAGS models written from the model specifications of Duncan (2026) | JAGS 4.x, rjags, coda | 10 to 20 minutes |
| `validate_convergence.R` | default integration grid against a finer grid; the difference must not exceed the reported integration error | none | a few minutes |
| `summarize.R` | collects the judgements into `results/SUMMARY.md` | none | seconds |

The rjags 4-x package on CRAN does not work with JAGS 5, so JAGS 4.x is
needed.

Each script writes `results/validate_<name>.csv` with one row per item: the
reference value, the bssd value, the tolerance and the judgement (PASS,
EXPLAINED, FAIL or INFO, as in `inst/reproduce`). Reasons for EXPLAINED items
are registered in `explained.csv`. For the MCMC comparisons the tolerance is
four Monte Carlo standard errors computed from the effective sample size.
Each script ends with the self-check of `inst/reproduce/judge.R`, which
shifts the reference values of the PASS items by five times their tolerance
and stops unless all of them become FAIL.

The folder `python/` contains the Python scripts that computed the expected
values of the tests in `tests/testthat/`.
