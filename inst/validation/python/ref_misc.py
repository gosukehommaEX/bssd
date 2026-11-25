"""Reference values for elicit_prior(), pooled_ratio(), count_oc() and
signal_count().

elicit_prior(): the two equations F(q1) = p1 and F(q2) = p2 are solved with
scipy.optimize.fsolve from several starting values (two-dimensional
Newton-type search; the package uses nested one-dimensional root finding).
count_oc(): exact sums of binomial and Poisson probabilities.
signal_count(): the posterior probability of canary_post() (ratio,
binomial) by the dense rule of ref_canary.py for consecutive counts, also
for a case in which the probability falls again for counts near n.
Output: ref_misc.json.
"""
import numpy as np
from scipy import stats
from scipy.optimize import fsolve, brentq
import ref_common as rc
import ref_canary as rcan


def solve(fun, starts):
    for s in starts:
        x, info, ier, msg = fsolve(fun, s, full_output=True, xtol=1e-14)
        if np.max(np.abs(fun(x))) < 1e-11:
            return x
    raise RuntimeError("no solution")


def elicit_gamma(q, p):
    f = lambda th: [stats.gamma.cdf(q[0], th[0], scale=1 / th[1]) - p[0],
                    stats.gamma.cdf(q[1], th[0], scale=1 / th[1]) - p[1]]
    x = solve(lambda x: np.array(f(np.exp(x))), [np.log([1, 1]),
              np.log([0.5, 0.5]), np.log([2, 1 / q[0]]),
              np.log([0.3, 0.2])])
    return np.exp(x).tolist()


def elicit_beta(q, p):
    f = lambda th: [stats.beta.cdf(q[0], th[0], th[1]) - p[0],
                    stats.beta.cdf(q[1], th[0], th[1]) - p[1]]
    x = solve(lambda x: np.array(f(np.exp(x))),
              [np.log([1, 10]), np.log([2, 50]), np.log([5, 200]),
               np.log([0.5, 0.5])])
    return np.exp(x).tolist()


def count_oc_fixed(yc, n, pc, pd, k):
    n1 = int(round(k / (k + 1) * n))
    n0 = n - n1
    x = np.arange(n1 + 1)
    tot = 0.0
    pmf0 = stats.binom.pmf(np.arange(n0 + 1), n0, pc)
    pmf1 = stats.binom.pmf(x, n1, pd)
    conv = np.convolve(pmf1, pmf0)
    return float(conv[yc:].sum())


if __name__ == "__main__":
    res = {}
    res["gamma_schnell"] = elicit_gamma([0.5, 2], [0.5, 0.75])
    res["gamma_gould_r2_z08"] = elicit_gamma([1, 2], [0.5, 0.8])
    res["beta_example"] = elicit_beta([0.02, 0.04], [0.5, 0.95])
    res["lnorm_mean"] = list(rcan.lnorm_default(1, 0.1, 10, 0.95))
    res["lnorm_mean_8"] = list(rcan.lnorm_default(1, 0.1, 8, 0.95))
    z = stats.norm.ppf([0.1, 0.9])
    s = (np.log(4) - np.log(0.5)) / (z[1] - z[0])
    res["lnorm_q"] = [np.log(0.5) - s * z[0], s]
    res["pooled_ratio"] = (250 * 1 / 2 + 80 * 2 / 3) / (250 / 2 + 80 / 3)
    res["count_oc_fixed"] = [count_oc_fixed(6, 100, 0.02, pd, 1)
                             for pd in (0.02, 0.04, 0.08)]
    res["count_oc_random"] = [float(stats.binom.sf(5, 100,
                                                   (0.02 + 2 * pd) / 3))
                              for pd in (0.02, 0.04, 0.08)]
    res["count_oc_pois"] = [float(stats.poisson.sf(214, 20000 *
                                                   (0.007 + R * 0.007) / 2))
                            for R in (1, 2, 4)]
    # signal_count: canary_post ratio, prior Beta(20, 980), cutoff 0.8,
    # N = 100: probabilities at the counts 5 and 6
    mu, sd = rcan.lnorm_default()
    probs = {}
    for y in (5, 6):
        cs = {"id": "S%d" % y, "y": y, "n": 100, "k": 1,
              "prior_c": [20, 980],
              "metric": "ratio", "crit": 1.0}
        probs[str(y)] = rcan.run_binom(cs)["prob"]
    res["signal_count_probs"] = probs
    # signal_count: prior Beta(1, 49), R > 2 with cutoff 0.8, N = 100: the
    # probability passes 0.8 between the counts 16 and 17 (it falls again
    # for counts near N, where the posterior of p_c lies outside the range
    # of this rule, so those counts are not computed here)
    probs = {}
    for y in (16, 17):
        cs = {"id": "T%d" % y, "y": y, "n": 100, "k": 1, "prior_c": [1, 49],
              "metric": "ratio", "crit": 2.0}
        probs[str(y)] = rcan.run_binom(cs)["prob"]
    res["signal_count_nonmonotone_probs"] = probs
    rc.dump(res, "ref_misc.json")
    for k, v in res.items():
        print(k, v)
