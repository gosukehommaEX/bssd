# Python reference values for the tests

These scripts computed the expected values written in `tests/testthat/`
before the R functions were run. They use methods that differ from the
package (NumPy 2.5 and SciPy 1.18 were used).

| Script | Method | Values used in |
|---|---|---|
| `ref_bdribs.py` | adaptive quadrature over p with `scipy.integrate.quad`; incidences of case A also by two-dimensional integration over (delta, p) | `test-bdribs_post.R` (cases A, B, C); output `ref_bdribs.json` |
| `ref_bdribs_fixed.py` | adaptive quadrature over p with a known background rate | `test-bdribs_post.R` (cases F1, F2, F3); output `ref_bdribs_fixed.json` |
| `ref_expbern2.py` | dense fixed Gauss-Legendre grid over (log lambda, logit p), split at logit(p0), 200 and 400 nodes per coordinate | `test-expbern_post.R` (settings S1, S2, S3); output `ref_expbern.json`. S1 and S2 were computed in a first run and S3 in a second run, which skips settings already in the JSON file |
| `ref_llbern.py` | all four parameters on their prior distribution-function scale, Gauss-Legendre nodes, mu_d split at the threshold of RR(t_star) | `test-llbern_post.R`: probabilities, posterior means and incidences with 48 nodes (the quantile part of this script was not used) |
| `ref_ll_q.py` | same as `ref_llbern.py` | `test-llbern_post.R`: median and upper limit of RR (32 nodes, run as `python3 ref_ll_q.py median` and `python3 ref_ll_q.py upper`) and the probability for k = 2 (`python3 ref_ll_q.py S2`, 48 nodes) |

Run each script from this folder, for example `python3 ref_bdribs.py`.
