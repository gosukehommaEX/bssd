# Reproduction of published results

The scripts in this folder reproduce numerical results of the articles whose
methods are implemented in bssd, with the functions of the package.

| Script | Article | Run time |
|---|---|---|
| `reproduce_Mukhopadhyay2018.R` | Mukhopadhyay, Waterhouse and Hartford (2018), Pharmaceutical Statistics 17(6):823-834: Table 1, hypothetical examples 1 and 2, operating characteristics in Section 5 | about one minute |
| `reproduce_Duncan2026.R` | Duncan (2026), Journal of Biopharmaceutical Statistics 36(2):227-244: Table 1, the descriptions of the scenarios in Section 4, Figure 1 | a few seconds |
| `reproduce_Duncan2026_simulation.R` | Duncan (2026): simulations for Figures 2 to 10 and the false positive rates in Section 4.2 (assumptions in the header of the script) | 2 to 3 hours with 1000 simulated trials per scenario on one core; the scenarios can run in parallel (`n_cores`) after `devtools::install()` |
| `reproduce_Duncan2026_figures.R` | Duncan (2026): figures and judgements from the simulation results | seconds |
| `reproduce_Duncan2026_prior.R` | Duncan (2026): simulations for Figures 4 to 6 with other priors of BDRIBS and EXP_BERN than in the main simulation, which explain the failed statements D_Q_F5a and D_Q_F6a (settings in the header of the script) | about 30 minutes with 2000 simulated trials per scenario on one core |
| `summarize.R` | collects the judgements into `results/SUMMARY.md` | seconds |

Run the scripts from the package root after `devtools::load_all()` (or
`library(bssd)`), in the order of the table, for example
`source("inst/reproduce/reproduce_Duncan2026.R")`.

Each script writes `results/reproduce_<Article>.csv` with one row per item and
the judgement:

- PASS: the reproduced value agrees with the published value within the
  tolerance given in the script.
- EXPLAINED: the values differ, and the reason is registered in
  `explained.csv`.
- FAIL: the values differ without a registered reason, or the reproduced
  value is NA.
- INFO: no published value to compare with, or a value shown for comparison
  only (the simulation results of Duncan, 2026, whose settings are not fully
  reported; the published values are read from the figures in
  `Duncan2026_digitized.csv`).

Each script ends with a self-check that shifts the published values of the
PASS items by five times their tolerance and stops unless all of them become
FAIL, so that the tolerances are not too loose to detect a difference. The
self-check also stops when no item is PASS (for example, when all tolerances
are NA by mistake), because the judgement is then vacuous.
`judge.R` contains the shared functions.
