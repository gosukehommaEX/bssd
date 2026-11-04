"""Reference values for bdribs_post with a known background rate (pt_hp = Inf).

Adaptive quadrature over p with scipy.integrate.quad; the posterior density of
p is Beta(a, b) times the Poisson likelihood with delta fixed.
"""
import json
import numpy as np
from scipy import integrate, optimize

Q = dict(limit=500, epsabs=0.0, epsrel=1e-12)


def ref(y, E, k, d0, a, b, rr0, t_star, level):
    lk = lambda p: (a - 1) * np.log(p) + (b - 1) * np.log1p(-p) - y * np.log1p(-p) - E * d0 / ((k + 1) * (1 - p))
    grid = np.linspace(1e-9, 1 - 1e-9, 400001)
    pm = grid[np.argmax(lk(grid) - (a - 1) * np.log(grid) - (b - 1) * np.log1p(-grid))]
    m = lk(pm)
    f = lambda p: np.exp(lk(p) - m)
    pc = k * rr0 / (k * rr0 + 1)

    def I(g, lo=0.0, hi=1.0):
        br = [lo] + [q for q in sorted({pc, pm}) if lo < q < hi] + [hi]
        return sum(integrate.quad(lambda p: g(p) * f(p), br[j], br[j + 1], **Q)[0] for j in range(len(br) - 1))

    Z = I(lambda p: 1.0)
    r = lambda p: p / (k * (1 - p))
    out = dict(prob=I(lambda p: 1.0, pc, 1.0) / Z, p_mean=I(lambda p: p) / Z, r_mean=I(r) / Z)
    out["lambda_d"] = out["r_mean"] * d0
    cdf = lambda x: I(lambda p: 1.0, 0.0, x) / Z
    for nm, lev in (("r_lower", (1 - level) / 2), ("r_median", 0.5), ("r_upper", (1 + level) / 2)):
        pq = optimize.brentq(lambda x: cdf(x) - lev, 1e-14, 1 - 1e-14, xtol=1e-15, rtol=1e-13)
        out[nm] = r(pq)
    out["ip_placebo"] = 1 - np.exp(-d0 * t_star)
    out["ip_drug"] = 1 - I(lambda p: np.exp(-r(p) * d0 * t_star)) / Z
    out["ip_difference"] = out["ip_drug"] - out["ip_placebo"]
    return out


if __name__ == "__main__":
    cases = {
        "F1": dict(y=9, E=540.0, k=2.0, d0=0.011, a=1.0, b=1.0, rr0=1.0, t_star=0.5, level=0.90),
        "F2": dict(y=9, E=540.0, k=2.0, d0=0.011, a=1.0, b=1.0, rr0=1.3, t_star=0.5, level=0.90),
        "F3": dict(y=3, E=900.0, k=1.0, d0=0.004, a=0.5, b=0.5, rr0=1.2, t_star=1.0, level=0.95),
    }
    res = {nm: ref(**c) for nm, c in cases.items()}
    for nm, v in res.items():
        print(nm, json.dumps(v))
    json.dump(res, open("ref_bdribs_fixed.json", "w"), indent=1)
