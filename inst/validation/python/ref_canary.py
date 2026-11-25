"""Reference values for canary_post() (Gould and Wang, 2017).

Binomial likelihood: the convolution of two binomial counts with fixed arm
sizes, on (x, z) = (logit p_c, log M) for the ratio and the odds ratio and
(logit p_c, logit p_d) for the difference, with composite Gauss-Legendre
rules split exactly at the threshold. Poisson likelihood: the Gamma prior
of the control rate is integrated out analytically and log R is integrated
by scipy.integrate.quad. Output: ref_canary.json.
"""
import numpy as np
from scipy import special, stats, integrate
from scipy.optimize import brentq
import ref_common as rc

CASES = [
    {"id": "C1", "y": 12, "n": 100, "k": 1, "prior_c": [1, 49],
     "metric": "ratio", "crit": 1.0},
    {"id": "C2", "y": 12, "n": 100, "k": 1, "prior_c": [1, 49],
     "metric": "difference", "crit": 0.02},
    {"id": "C3", "y": 20, "n": 120, "k": 2, "prior_c": [5, 45],
     "metric": "odds", "crit": 1.2},
    {"id": "C4", "y": 250, "n": 20000, "k": 1, "prior_c": [75.71, 8500],
     "metric": "ratio", "likelihood": "poisson", "crit": 2.0},
]


def lnorm_default(mean=1.0, lo=0.1, hi=10.0, p=0.95):
    def f(s):
        mu = np.log(mean) - s * s / 2
        return (stats.norm.cdf((np.log(hi) - mu) / s) -
                stats.norm.cdf((np.log(lo) - mu) / s) - p)
    s = brentq(f, 1e-3, 10, xtol=1e-15)
    return np.log(mean) - s * s / 2, s


def softplus(x):
    return np.logaddexp(0, x)


def make_logpost(cs, mu, sd):
    Y, N, k = cs["y"], cs["n"], cs["k"]
    a, b = cs["prior_c"]
    n1 = int(round(k / (k + 1) * N))
    n0 = N - n1
    xs = np.arange(max(0, Y - n0), min(Y, n1) + 1)
    metric = cs["metric"]

    def logpost(x, Z):
        X = np.broadcast_to(x[:, None], Z.shape)
        lpc, l1pc = -softplus(-X), -softplus(X)
        with np.errstate(all="ignore"):
            if metric == "ratio":
                lpt = lpc + Z
                ok = lpt < 0
                lpt = np.where(ok, lpt, -1.0)
                l1pt = np.log1p(-np.exp(lpt))
            elif metric == "odds":
                lpt, l1pt = -softplus(-(X + Z)), -softplus(X + Z)
                ok = np.ones(Z.shape, bool)
            else:
                lpt, l1pt = -softplus(-Z), -softplus(Z)
                ok = np.ones(Z.shape, bool)
            terms = [special.gammaln(n1 + 1) - special.gammaln(x + 1) -
                     special.gammaln(n1 - x + 1) + x * lpt + (n1 - x) * l1pt +
                     special.gammaln(n0 + 1) - special.gammaln(Y - x + 1) -
                     special.gammaln(n0 - Y + x + 1) + (Y - x) * lpc +
                     (n0 - Y + x) * l1pc for x in xs]
            ll = special.logsumexp(np.stack(terms), axis=0)
            lp = ll + a * lpc + b * l1pc
            if metric == "ratio":
                lp = lp + stats.norm.logpdf(Z, mu, sd) - \
                    stats.norm.logcdf((-lpc - mu) / sd)
                lp = np.where(ok, lp, -np.inf)
            elif metric == "odds":
                lp = lp + stats.norm.logpdf(Z, mu, sd)
            else:
                # uniform prior of p_d on (0, 1), in logit coordinates
                lp = lp + lpt + l1pt
        return lp
    return logpost


def run_binom(cs):
    mu, sd = lnorm_default()
    lp = make_logpost(cs, mu, sd)
    a, b = cs["prior_c"]
    xl = special.logit(stats.beta.ppf(1e-13, a, b))
    xh = special.logit(stats.beta.ppf(1 - 1e-13, a, b))
    if cs["metric"] == "difference":
        zl, zh = -40.0, 40.0
    else:
        zl, zh = mu - 11 * sd, mu + 11 * sd
    xs = np.linspace(xl, xh, 121)
    zs = np.linspace(zl, zh, 401)
    L = lp(xs, np.broadcast_to(zs, (xs.size, zs.size)))
    keep = L > L.max() - 40
    zk = zs[np.any(keep, axis=0)]
    S = rc.Split2D(lp, np.linspace(xl, xh, 41), zl, zh,
                   (zk.min() - 0.5, zk.max() + 0.5), n_core=120, n_tail=40)
    sig = lambda v: 1 / (1 + np.exp(-v))
    out = {"id": cs["id"]}
    lev = [0.025, 0.5, 0.975]
    if cs["metric"] == "difference":
        split = lambda d: special.logit(np.clip(sig(S.x) + d, 1e-300,
                                                1 - 1e-16))
        out["prob"] = 1 - S.below(split(cs["crit"]))
        out["mean_metric"] = S.mean(lambda x, z: sig(z) - sig(x))
        out["q"] = [S.quantile(split, p, -0.999, 0.999) for p in lev]
    else:
        out["prob"] = 1 - S.below(np.log(cs["crit"]))
        out["mean_metric"] = S.mean(lambda x, z: np.exp(z))
        out["q"] = [S.quantile(lambda c: np.log(c), p, 1e-6, 1e4)
                    for p in lev]
    out["mean_p_c"] = S.mean(lambda x, z: sig(x))
    return out


def run_pois(cs):
    mu, sd = lnorm_default()
    Y, E, k = cs["y"], cs["n"], cs["k"]
    a, b = cs["prior_c"]
    m_of = lambda v: E * (1 + k * np.exp(v)) / (k + 1)
    f0 = lambda v: (Y * np.log(m_of(v)) - (a + Y) * np.log(b + m_of(v)) +
                    stats.norm.logpdf(v, mu, sd))
    vs = np.linspace(mu - 12 * sd, mu + 12 * sd, 20001)
    vm = vs[np.argmax(f0(vs))]
    c0 = f0(vm)
    f = lambda v: np.exp(f0(v) - c0)
    q = lambda g, lo, hi: integrate.quad(g, lo, hi, points=None, limit=500,
                                         epsabs=0, epsrel=1e-12)[0]
    lo, hi = mu - 14 * sd, mu + 14 * sd
    Zn = q(f, lo, vm) + q(f, vm, hi)
    cdf = lambda v: (q(f, lo, min(v, vm)) +
                     (q(f, vm, v) if v > vm else 0)) / Zn
    out = {"id": cs["id"], "prob": 1 - cdf(np.log(cs["crit"]))}
    out["mean_metric"] = (q(lambda v: np.exp(v) * f(v), lo, vm) +
                          q(lambda v: np.exp(v) * f(v), vm, hi)) / Zn
    out["q"] = [np.exp(brentq(lambda v: cdf(v) - p, lo + 1, hi - 1,
                              xtol=1e-13)) for p in [0.025, 0.5, 0.975]]
    lam = lambda v: (a + Y) / (b + m_of(v))
    out["mean_lambda_c"] = (q(lambda v: lam(v) * f(v), lo, vm) +
                            q(lambda v: lam(v) * f(v), vm, hi)) / Zn
    return out


if __name__ == "__main__":
    res = {"lnorm_default": list(lnorm_default()), "cases": CASES,
           "results": [run_pois(c) if c.get("likelihood") == "poisson" else
                       run_binom(c) for c in CASES]}
    rc.dump(res, "ref_canary.json")
    print(res["lnorm_default"])
    for r in res["results"]:
        print(r)
