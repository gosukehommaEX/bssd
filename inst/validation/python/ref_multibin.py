"""Reference values for multibin_post() (Goren et al., 2020).

The posterior density of w = log(delta) is computed with nested adaptive
quadrature: for each w, the control probability of each study is
integrated out by scipy.integrate.quad on (0, min(1, exp(-w))); the outer
integrals over w use scipy.integrate.quad split at the mode and at the
threshold, over 12 prior standard deviations or more. Output:
ref_multibin.json.
"""
import numpy as np
from scipy import special, stats, integrate
from scipy.optimize import brentq
import ref_common as rc

CASES = [
    {"id": "G1", "y": [6, 10, 25, 47], "n": [200, 400, 800, 1600],
     "k": [1, 1, 1, 1], "prior_c": [(2, 98)] * 4, "lrr": [0, 5],
     "rr0": 1.0},
    {"id": "G2", "y": [3, 9], "n": [150, 300], "k": [1, 2],
     "prior_c": [(2, 98), (4, 196)], "lrr": [0, 2], "rr0": 1.5},
    # informative priors of log(delta) that conflict with the data: the
    # posterior lies far in the upper (G3) or lower (G4) tail of the prior
    {"id": "G3", "y": [60, 120], "n": [1000, 2000], "k": [1, 1],
     "prior_c": [(20, 1980), (40, 3960)], "lrr": [0, 0.1], "rr0": 1.5},
    {"id": "G4", "y": [2, 3], "n": [1000, 2000], "k": [1, 1],
     "prior_c": [(200, 19800), (400, 39600)], "lrr": [0, 0.1],
     "rr0": 0.8},
]


def log_inner(w, y, n, r, a, b):
    c = 1 - r + r * np.exp(w)
    up = min(1.0, np.exp(-w))

    def g0(th):
        return ((a + y - 1) * np.log(th) + (b - 1) * np.log1p(-th) +
                (n - y) * np.log1p(-min(c * th, 1.0)))
    # mode on a grid for scaling and as a break point
    grid = np.linspace(up * 1e-6, up * (1 - 1e-9), 4001)
    vals = np.array([g0(x) for x in grid])
    m = vals.max()
    tm = grid[np.argmax(vals)]
    f = lambda th: np.exp(g0(th) - m)
    pts = [tm]
    v = (integrate.quad(f, 0, tm, limit=400, epsabs=0, epsrel=1e-12)[0] +
         integrate.quad(f, tm, up, limit=400, epsabs=0, epsrel=1e-12)[0])
    return np.log(v) + m - special.betaln(a, b) + y * np.log(c)


def make_logpost(cs):
    r = [k / (k + 1) for k in cs["k"]]
    m, s = cs["lrr"]

    def lp(w):
        tot = stats.norm.logpdf(w, m, s)
        for y, n, rr, (a, b) in zip(cs["y"], cs["n"], r, cs["prior_c"]):
            tot += log_inner(w, y, n, rr, a, b)
        return tot
    return lp


def run_case(cs):
    lp = make_logpost(cs)
    ws = np.linspace(-6, 6, 241)
    lv = np.array([lp(w) for w in ws])
    wm = ws[np.argmax(lv)]
    c0 = lv.max()
    f = lambda w: np.exp(lp(w) - c0)
    q = lambda lo, hi: integrate.quad(f, lo, hi, limit=400, epsabs=0,
                                      epsrel=1e-10)[0]
    m, s = cs["lrr"]
    lo, hi = wm - 6, wm + 6
    # outer limits: 12 prior standard deviations, and at least the central
    # range around the mode (for an informative prior that conflicts with
    # the data)
    a_lo, b_hi = min(m - 12 * s, lo), max(m + 12 * s, hi)
    tails = (integrate.quad(f, a_lo, lo, limit=200, epsabs=0,
                            epsrel=1e-8)[0],
             integrate.quad(f, hi, b_hi, limit=200, epsabs=0,
                            epsrel=1e-8)[0])
    Z = q(lo, wm) + q(wm, hi) + tails[0] + tails[1]

    def cdf(w):
        if w <= lo:
            return integrate.quad(f, a_lo, w, limit=200)[0] / Z
        if w <= wm:
            return (tails[0] + q(lo, w)) / Z
        if w <= hi:
            return (tails[0] + q(lo, wm) + q(wm, w)) / Z
        return 1 - integrate.quad(f, w, b_hi, limit=200)[0] / Z
    out = {"id": cs["id"], "prob": 1 - cdf(np.log(cs["rr0"])),
           "tail_mass": [t / Z for t in tails]}
    out["q"] = [np.exp(brentq(lambda w: cdf(w) - p, lo, hi, xtol=1e-12))
                for p in [0.025, 0.5, 0.975]]
    wf = lambda w: w * f(w)
    out["mean_log_delta"] = (
        integrate.quad(wf, a_lo, lo, limit=200, epsabs=0,
                       epsrel=1e-8)[0] +
        integrate.quad(wf, lo, wm, limit=400)[0] +
        integrate.quad(wf, wm, hi, limit=400)[0] +
        integrate.quad(wf, hi, b_hi, limit=200, epsabs=0,
                       epsrel=1e-8)[0]) / Z
    return out


if __name__ == "__main__":
    res = {"cases": CASES, "results": [run_case(c) for c in CASES]}
    rc.dump(res, "ref_multibin.json")
    for r in res["results"]:
        print(r)
