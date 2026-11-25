# Validation

The scripts in this folder compare the functions of bssd with independent
computations. Run them from the package root after `devtools::load_all()` (or
`library(bssd)`), for example `source("inst/validation/validate_jags.R")`.

| Script | Comparison | Requirements | Run time |
|---|---|---|---|
| `validate_bdribs_package.R` | `bdribs_post()` against the bdribs package (MCMC by the authors of the BDRIBS method), with a known and with a Gamma-distributed background rate | JAGS 4.x, rjags, coda, bdribs | a few minutes |
| `validate_jags.R` | `expbern_post()` and `llbern_post()` against MCMC with JAGS models written from the model specifications of Duncan (2026) | JAGS 4.x, rjags, coda | 10 to 20 minutes |
| `validate_convergence.R` | default integration grid against a finer grid (for LL_BERN, 24 Gauss-Legendre nodes and 721 Chebyshev points) and, for one LL_BERN data set, against an independent Python reference (`python/ref_llbern_breakpoints.py`); the difference must not exceed the reported integration error | none | a few minutes |
| `validate_sim_oc.R` | simulation results of `sim_oc()` (`inst/reproduce/results/Duncan2026_simulation.rds`) against an independent Python simulation (`python/ref_sim_oc.py`): power of BDRIBS(1) and BDRIBS(2) in Figures 3, 4, 5 and 9, power and risk difference of EXP_BERN in scenario 5.2; tolerance four standard errors of the difference of two Monte Carlo estimates | none (reads `python/ref_sim_oc.csv`) | seconds |
| `validate_methods_python.R` | the functions added in 0.2.0 (`exptime_post()`, `binbeta_post()`, `canary_post()`, `mapbern_post()`, `map_beta()`, `multibin_post()`) against independent Python computations (`python/ref_*.py`, collected in `python/ref_methods.csv`): composite Gauss-Legendre rules split exactly at the thresholds, an exact finite mixture of Beta products for `binbeta_post()` and nested adaptive quadrature for `multibin_post()`; tolerance 2e-6 for probabilities and 1e-5 relative for means and quantiles | none (reads `python/ref_methods.csv`) | about one minute |
| `validate_jags_methods.R` | `exptime_post()`, `binbeta_post()`, `mapbern_post()`, `multibin_post()` and `canary_post()` (Poisson likelihood) against MCMC with JAGS models written from the model specifications of the articles | JAGS 4.x, rjags, coda | 5 to 10 minutes |
| `validate_gould2017_author_grid.R` | the critical counts of Tables III and V of Gould and Wang (2017) computed by a transcription of the authors' grid computations (`python/ref_gould2017_author_grid.csv`) against the published counts, which explains the differences from the exact counts of bssd in `inst/reproduce/reproduce_Gould2017.R` | none | seconds |
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
and stops unless all of them become FAIL. The self-check also stops when no
item is PASS (for example, when all tolerances are NA by mistake), because the
judgement is then vacuous.

The folder `python/` contains the Python scripts that computed the expected
values of the tests in `tests/testthat/` and the reference values used by
`validate_convergence.R`, `validate_sim_oc.R` and
`validate_methods_python.R`. The transcription of the code of Gould and Wang
(2017) is not in the package, because the authors' code has no stated
license; only its counts are included.

For the same reason, the authors' own code of Schnell and Ball (2016) and of
Gould and Wang (2017) is run by scripts outside the package
(`Rpackage/tools/external/` of the project). Their results are copied to
`results/schnell2016_author_code.csv` (posterior summaries of Table 2 from
the authors' sampler with the published seed and with long chains, from a
sampler with the subject-level likelihood, and from bssd) and
`results/gould2017_author_code.csv` (critical counts of Tables III and V
from the authors' code). They support the explanations in
`inst/reproduce/explained.csv`.
