"""Reference values for binbeta_post() (Ball and Schnell, 2016).

Exact computation: with pi = w_c theta_c + w_d theta_d and
1 - pi = w_c (1 - theta_c) + w_d (1 - theta_d), the binomial expansion of
pi^Z (1 - pi)^(N - Z) makes the posterior a finite mixture of products of
Beta distributions, theta_c ~ Beta(a + i, b + j) and
theta_d ~ Beta(c + Z - i, d + N - Z - j), with weights proportional to
C(Z, i) C(N - Z, j) w_c^(i + j) w_d^(N - i - j) B(a + i, b + j)
B(c + Z - i, d + N - Z - j). P(theta_d > p0) is a weighted sum of Beta tail
probabilities; P(theta_d > c theta_c) uses scipy.integrate.quad for each
component. Output: ref_binbeta.json.
"""
import numpy as np
from scipy import special, stats, integrate
from scipy.optimize import brentq
import ref_common as rc

CASES = [
    {"id": "B1", "y": 40, "n": 300, "k": 1, "prior_c": [31, 271],
     "prior_d": [1, 1], "p0": 0.15, "rr0": 1.0},
    {"id": "B2", "y": 12, "n": 80, "k": 2, "prior_c": [5, 45],
     "prior_d": [0.5, 0.5], "p0": None, "rr0": 1.5},
    # informative priors of both proportions that conflict with the data:
    # the posterior of theta_d lies in the upper tail of its prior (B3, B4)
    # or in the lower tail (B5); P(theta_d > theta_c) = 0.5 in B3 by
    # symmetry
    {"id": "B3", "y": 12, "n": 100, "k": 1, "prior_c": [1, 49],
     "prior_d": [1, 49], "p0": 0.05, "rr0": 1.0},
    {"id": "B4", "y": 100, "n": 100, "k": 1, "prior_c": [1, 99],
     "prior_d": [1, 99], "p0": 0.3, "rr0": 2.0},
    {"id": "B5", "y": 5, "n": 1000, "k": 1, "prior_c": [5, 495],
     "prior_d": [30, 270], "p0": 0.04, "rr0": 8.0},
]


def components(cs):
    Z, N, k = cs["y"], cs["n"], cs["k"]
    a, b = cs["prior_c"]
    c, d = cs["prior_d"]
    wc = 1 / (k + 1)
    i = np.arange(Z + 1)[:, None]
    j = np.arange(N - Z + 1)[None, :]
    lw = (special.gammaln(Z + 1) - special.gammaln(i + 1) -
          special.gammaln(Z - i + 1) + special.gammaln(N - Z + 1) -
          special.gammaln(j + 1) - special.gammaln(N - Z - j + 1) +
          (i + j) * np.log(wc) + (N - i - j) * np.log(1 - wc) +
          special.betaln(a + i, b + j) +
          special.betaln(c + Z - i, d + N - Z - j))
    w = np.exp(lw - lw.max())
    w = w / w.sum()
    keep = w > 1e-18
    I = np.broadcast_to(i, w.shape)[keep]
    J = np.broadcast_to(j, w.shape)[keep]
    return (w[keep], a + I, b + J, c + Z - I, d + N - Z - J)


def run_case(cs):
    w, a1, b1, a2, b2 = components(cs)
    out = {"id": cs["id"], "n_components": int(w.size)}
    if cs["p0"] is not None:
        out["prob_rate"] = float(np.sum(w * stats.beta.sf(cs["p0"], a2, b2)))

    def p_ratio_le(c):
        # one adaptive integral over theta_c of the mixture density times
        # the mixture distribution function of theta_d at c * theta_c
        f = lambda u: float(np.sum(w * stats.beta.pdf(u, a1, b1) *
                                   stats.beta.cdf(min(c * u, 1.0), a2, b2)))
        mode = float(np.sum(w * (a1 - 1) / (a1 + b1 - 2)))
        pts = sorted({mode, min(1 / c, 1.0) - 1e-12})
        pts = [x for x in pts if 0 < x < 1]
        return integrate.quad(f, 0, 1, points=pts, limit=500, epsabs=1e-15,
                              epsrel=1e-12)[0]
    # components with negligible weight are dropped for the ratio
    sel = w > 1e-12
    w, a1, b1, a2, b2 = [v[sel] for v in (w, a1, b1, a2, b2)]
    w = w / w.sum()
    out["prob_ratio"] = 1 - p_ratio_le(cs["rr0"])
    out["mean_theta_c"] = float(np.sum(w * a1 / (a1 + b1)))
    out["mean_theta_d"] = float(np.sum(w * a2 / (a2 + b2)))
    # E(1 / theta_c) is infinite in components with a1 <= 1
    if np.all(a1 > 1):
        out["mean_rr"] = float(np.sum(w * a2 / (a2 + b2) * (a1 + b1 - 1) /
                                      (a1 - 1)))
    cdf_d = lambda x: float(np.sum(w * stats.beta.cdf(x, a2, b2)))
    out["q_theta_d"] = [brentq(lambda x: cdf_d(x) - p, 1e-12, 1 - 1e-12,
                               xtol=1e-13) for p in [0.025, 0.5, 0.975]]
    out["q_rr"] = [brentq(lambda c: p_ratio_le(c) - p, 0.01, 100,
                          xtol=1e-11) for p in [0.025, 0.5, 0.975]]
    return out


if __name__ == "__main__":
    res = {"cases": CASES, "results": [run_case(c) for c in CASES]}
    rc.dump(res, "ref_binbeta.json")
    for r in res["results"]:
        print(r)
