"""Reference values for expbern_post with a dense fixed Gauss-Legendre grid.

Natural parameters (lambda, p) are integrated on (log lambda, logit p) over a
fixed wide box with many nodes; the box is split at logit(p0) so that
P(r > rr0) is a smooth integral. No box refinement and no Chebyshev
interpolation are used (unlike the package). Quantiles of r come from the
marginal distribution function of logit p, computed on piecewise GL panels.
"""
import json
import numpy as np
from scipy import optimize, special, stats

rng = np.random.default_rng(20261104)
N = 60
t = np.round(rng.uniform(5, 52, N), 1)
grp = rng.random(N) < 0.5
haz = np.where(grp, 2.5 * 0.004, 0.004)
tev = rng.exponential(1 / haz)
y = (tev <= t).astype(int)


def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * x + (a + b) / 2, w * (b - a) / 2


def loglik_grid(lam, r, wp):
    wd = 1 - wp
    out = np.empty((len(lam), len(r)))
    for i, l in enumerate(lam):
        h1 = l * t
        h2 = np.outer(r, h1)
        m1 = wp * (-np.expm1(-h1))[None, :] + wd * (-np.expm1(-h2))
        l0 = np.logaddexp(np.log(wp) - h1[None, :], np.log(wd) - h2)
        out[i] = np.where(y[None, :] == 1, np.log(m1), l0).sum(1)
    return out


def ref(k, lam_hp, pt_hp, df, a, b, rr0, t_star, n_l=400, n_z=400):
    al, be = lam_hp * pt_hp / df, pt_hp / df
    wp = 1 / (k + 1)
    u_lo, u_hi = np.log(stats.gamma.ppf(1e-13, al, scale=1 / be)), np.log(stats.gamma.ppf(1 - 1e-13, al, scale=1 / be))
    z_lo, z_hi = -60.0, 60.0
    z0 = np.log(k * rr0)
    u, wu = gl(n_l, u_lo, u_hi)
    lam = np.exp(u)
    lpu = stats.gamma.logpdf(lam, al, scale=1 / be) + u

    def panel(za, zb, n):
        # composite rule: breakpoints at -10 and 10 inside (za, zb)
        br = [za] + [v for v in (-10.0, 10.0) if za < v < zb] + [zb]
        zz, ww = zip(*[gl(n, br[j], br[j + 1]) for j in range(len(br) - 1)])
        z, wz = np.concatenate(zz), np.concatenate(ww)
        # stable log density of logit(p) under the Beta(a, b) prior
        lpz = a * special.log_expit(z) + b * special.log_expit(-z) - special.betaln(a, b)
        L = loglik_grid(lam, np.exp(z) / k, wp) + lpu[:, None] + lpz[None, :]
        return z, wz, L

    zA, wA, LA = panel(z_lo, z0, n_z)
    zB, wB, LB = panel(z0, z_hi, n_z)
    M = max(LA.max(), LB.max())
    WA = np.exp(LA - M) * wu[:, None] * wA[None, :]
    WB = np.exp(LB - M) * wu[:, None] * wB[None, :]
    Z = WA.sum() + WB.sum()
    prob = WB.sum() / Z
    W = np.hstack([WA, WB]) / Z
    z = np.concatenate([zA, zB])
    r = np.exp(z) / k
    out = dict(prob=prob,
               p_mean=(W * special.expit(z)[None, :]).sum(),
               r_mean=(W * r[None, :]).sum(),
               lambda_p=(W * lam[:, None]).sum(),
               lambda_d=(W * np.outer(lam, r)).sum(),
               ip_placebo=(W * (-np.expm1(-lam * t_star))[:, None]).sum(),
               ip_drug=(W * (-np.expm1(-np.outer(lam, r) * t_star))).sum())
    out["ip_difference"] = out["ip_drug"] - out["ip_placebo"]

    def cdf(zq):
        zc, wc, Lc = panel(z_lo, zq, 200)
        return (np.exp(Lc - M) * wu[:, None] * wc[None, :]).sum() / Z

    for nm, lev in (("r_lower", 0.025), ("r_median", 0.5), ("r_upper", 0.975)):
        zq = optimize.brentq(lambda v: cdf(v) - lev, z_lo + 1, z_hi - 1, xtol=1e-12)
        out[nm] = float(np.exp(zq) / k)
    return out


settings = {
    "S1": dict(k=1.0, lam_hp=0.004, pt_hp=2000.0, df=2.0, a=0.5, b=0.5, rr0=1.0, t_star=52.0),
    "S2": dict(k=2.0, lam_hp=0.005, pt_hp=3000.0, df=1.0, a=1.0, b=1.0, rr0=1.5, t_star=26.0),
    "S3": dict(k=1.0, lam_hp=0.004, pt_hp=1000.0, df=1.0, a=1.0, b=3.0, rr0=1.0, t_star=52.0),
}
if __name__ == "__main__":
    print("events", int(y.sum()), "of", N)
    import os
    res = json.load(open("ref_expbern.json")) if os.path.exists("ref_expbern.json") else dict(t=t.tolist(), y=y.tolist())
    for nm, s in settings.items():
        if nm in res:
            continue
        for n in (200, 400):
            o = ref(**s, n_l=n, n_z=n)
            print(nm, n, json.dumps({k_: round(v, 10) for k_, v in o.items()}), flush=True)
        res[nm] = o
    json.dump(res, open("ref_expbern.json", "w"), indent=1)
