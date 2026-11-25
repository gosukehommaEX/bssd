# Python reference values

These scripts computed the expected values written in `tests/testthat/`
before the R functions were run, and the reference values used by
`inst/validation/validate_convergence.R`,
`inst/validation/validate_sim_oc.R` and
`inst/validation/validate_methods_python.R`. They use methods that differ from
the package (NumPy 2.5 and SciPy 1.18 were used).

| Script | Method | Values used in |
|---|---|---|
| `ref_bdribs.py` | adaptive quadrature over p with `scipy.integrate.quad`; incidences of case A also by two-dimensional integration over (delta, p) | `test-bdribs_post.R` (cases A, B, C); output `ref_bdribs.json` |
| `ref_bdribs_fixed.py` | adaptive quadrature over p with a known background rate | `test-bdribs_post.R` (cases F1, F2, F3); output `ref_bdribs_fixed.json` |
| `ref_expbern2.py` | dense fixed Gauss-Legendre grid over (log lambda, logit p), split at logit(p0), 200 and 400 nodes per coordinate | `test-expbern_post.R` (settings S1, S2, S3); output `ref_expbern.json`. S1 and S2 were computed in a first run and S3 in a second run, which skips settings already in the JSON file |
| `ref_llbern.py` | all four parameters on their prior distribution-function scale, Gauss-Legendre nodes, mu_d split at the threshold of RR(t_star) | `test-llbern_post.R`: probabilities, posterior means and incidences with 48 nodes (the quantile part of this script was not used) |
| `ref_ll_q.py` | same as `ref_llbern.py` | `test-llbern_post.R`: median and upper limit of RR (32 nodes, run as `python3 ref_ll_q.py median` and `python3 ref_ll_q.py upper`) and the probability for k = 2 (`python3 ref_ll_q.py S2`, 48 nodes) |
| `ref_llbern_breakpoints.py` | log(mu_p) and log(sigma_p) with Gauss-Legendre nodes; sigma_d on its prior distribution-function scale with panels graded towards zero; mu_d with composite Gauss-Legendre panels whose breakpoints follow the logarithms of the follow-up times and the threshold of RR(t_star); no box refinement and no Chebyshev interpolation | `validate_convergence.R`: data set L1 (`data_llbern_L1.csv`, written by R with `sim_blinded()`); output `ref_llbern_breakpoints.csv`, which also contains the values for the data of `test-llbern_post.R` next to the reference values of `ref_llbern.py` and `ref_ll_q.py`; about 15 minutes |
| `ref_sim_oc.py` | trials simulated with NumPy from the design of `sim_blinded()`; BDRIBS with the placebo rate integrated out analytically and the Simpson rule over logit(p); EXP_BERN with the Simpson rule over (log(lambda_p), logit(p)) on the region located by a coarse grid; `python3 ref_sim_oc.py check` compares these functions with the reference values of `test-bdribs_post.R` and `test-expbern_post.R` | `validate_sim_oc.R`: Figures 3, 4, 5 and 9 (BDRIBS, 2000 trials) and scenario 5.2 (EXP_BERN, 1000 trials) of `reproduce_Duncan2026_simulation.R`; output `ref_sim_oc.csv`; about 30 minutes |

| `ref_exptime.py` | composite Gauss-Legendre rules on (log lambda_c, log lambda_d), split exactly at the threshold of the rate and of the rate ratio | `test-exptime_post.R` and `validate_methods_python.R` (data sets E1 and E2); output `ref_exptime.json` |
| `ref_binbeta.py` | exact finite mixture of products of Beta distributions from the binomial expansion of the likelihood | `test-binbeta_post.R` and `validate_methods_python.R` (cases B1 and B2); output `ref_binbeta.json` |
| `ref_canary.py` | binomial likelihood: composite Gauss-Legendre rules on (logit p_c, log M) or (logit p_c, logit p_d), split exactly at the threshold; Poisson likelihood: control rate integrated out analytically and `scipy.integrate.quad` over log R | `test-canary_post.R` and `validate_methods_python.R` (cases C1 to C4); output `ref_canary.json` |
| `ref_mapbern.py` | composite Gauss-Legendre rules on (logit pi_c, logit pi_d) split exactly at the threshold; for `map_beta()`, a dense grid over (logit a, logit b), and `scipy.integrate.quad` over pi_j for subject-level historical data | `test-mapbern_post.R`, `test-map_beta.R` and `validate_methods_python.R` (data set M1 and the MAP priors); output `ref_mapbern.json` |
| `ref_multibin.py` | nested adaptive quadrature with `scipy.integrate.quad` over the control probability of each study and over log delta, split at the threshold | `test-multibin_post.R` and `validate_methods_python.R` (cases G1 and G2); output `ref_multibin.json` |
| `ref_misc.py` | `scipy.optimize.fsolve` for `elicit_prior()`; the closed form for `pooled_ratio()`; exact sums of binomial and Poisson probabilities for `count_oc()`; the probabilities of `ref_canary.py` for `signal_count()` | `test-elicit_prior.R`, `test-pooled_ratio.R`, `test-count_oc.R` and `test-signal_count.R`; output `ref_misc.json` |
| `ref_methods_table.py` | collects the JSON files of the scripts above | `validate_methods_python.R`; outputs `ref_methods.csv` and `data_methods.csv` (the data sets of `ref_common.py`) |

`ref_common.py` contains the integration rules and the data sets (fixed
seeds) shared by the scripts for the functions added in bssd 0.2.0.
`ref_gould2017_author_grid.csv` contains the critical counts of Gould and
Wang (2017) computed by a transcription of the authors' code, which is kept
outside the package because the code has no stated license; it is read by
`inst/validation/validate_gould2017_author_grid.R`.

Run each script from this folder, for example `python3 ref_bdribs.py`.
