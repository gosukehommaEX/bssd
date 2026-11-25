"""Reference values for exptime_post() (Schnell and Ball, 2016).

Coordinates (x, z) = (log lambda_c, log lambda_d), composite Gauss-Legendre
rules, split exactly at z = log(rate0) for P(lambda_d > rate0) and at
z = x + log(c) for P(lambda_d / lambda_c > c). Output: ref_exptime.json.
"""
import numpy as np
from scipy import special, stats
import ref_common as rc

CASES = [
    {"id": "E1", "data": "E1", "k": 2, "prior_c": [30, 1500],
     "prior_d": [2, 70], "rate0": 0.03, "rr0": 1.2},
    {"id": "E2", "data": "E2", "k": 1, "prior_c": [50, 300],
     "prior_d": [1, 5], "rate0": None, "rr0": 1.0},
]


def make_logpost(y, t, k, pc, pd):
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    wc = 1 / (k + 1)

    def kern(lr, lam_log):
        # y * log(lambda) - lambda * t for log(lambda) = lam_log (matrix)
        return y * lam_log[..., None] - np.exp(lam_log)[..., None] * t

    def logpost(x, Z):
        X = np.broadcast_to(x[:, None], Z.shape)
        a = np.log(wc) + kern(None, X)
        b = np.log(1 - wc) + kern(None, Z)
        ll = np.logaddexp(a, b).sum(axis=-1)
        lpx = pc[0] * np.log(pc[1]) - special.gammaln(pc[0]) + pc[0] * X - \
            pc[1] * np.exp(X)
        lpz = pd[0] * np.log(pd[1]) - special.gammaln(pd[0]) + pd[0] * Z - \
            pd[1] * np.exp(Z)
        return ll + lpx + lpz
    return logpost


def run_case(cs, data):
    d = data[cs["data"]]
    lp = make_logpost(d["y"], d["t"], cs["k"], cs["prior_c"], cs["prior_d"])
    # range of x from the prior of lambda_c (informative) and of z from a
    # coarse scan
    xl = np.log(stats.gamma.ppf(1e-13, cs["prior_c"][0],
                                scale=1 / cs["prior_c"][1]))
    xh = np.log(stats.gamma.ppf(1 - 1e-13, cs["prior_c"][0],
                                scale=1 / cs["prior_c"][1]))
    zl = np.log(stats.gamma.ppf(1e-16, cs["prior_d"][0],
                                scale=1 / cs["prior_d"][1])) - 2
    zh = np.log(stats.gamma.ppf(1 - 1e-16, cs["prior_d"][0],
                                scale=1 / cs["prior_d"][1])) + 1
    xs = np.linspace(xl, xh, 121)
    zs = np.linspace(zl, zh, 241)
    L = lp(xs, np.broadcast_to(zs, (xs.size, zs.size)))
    keep = L > L.max() - 40
    zc = zs[np.any(keep, axis=0)]
    core = (zc.min() - 0.5, zc.max() + 0.5)
    S = rc.Split2D(lp, np.linspace(xl, xh, 41), zl, zh, core, n_core=100,
                   n_tail=30)
    out = {"id": cs["id"]}
    lev = [0.025, 0.5, 0.975]
    if cs["rate0"] is not None:
        out["prob_rate"] = 1 - S.below(np.log(cs["rate0"]))
    out["prob_ratio"] = 1 - S.below(S.x + np.log(cs["rr0"]))
    out["mean_lambda_c"] = S.mean(lambda x, z: np.exp(x))
    out["mean_lambda_d"] = S.mean(lambda x, z: np.exp(z))
    out["mean_r"] = S.mean(lambda x, z: np.exp(z - x))
    out["q_lambda_d"] = [S.quantile(lambda c: np.log(c), p, 1e-8, 10)
                         for p in lev]
    out["q_r"] = [S.quantile(lambda c: S.x + np.log(c), p, 1e-6, 1e4)
                  for p in lev]
    return out


if __name__ == "__main__":
    data = rc.make_data()
    res = {"cases": CASES, "results": [run_case(c, data) for c in CASES]}
    rc.dump(res, "ref_exptime.json")
    for r in res["results"]:
        print(r)
