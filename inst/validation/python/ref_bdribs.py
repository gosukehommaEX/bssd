"""Independent reference values for bdribs_post (integration over p, not logit p)."""
import json
import numpy as np
from scipy import integrate, optimize, special, stats

QK = dict(limit=500, epsabs=0.0, epsrel=1e-12)


def bdribs_ref(y, E, k, lam_hp, pt_hp, df, a, b, rr0, t_star):
    al, be = lam_hp * pt_hp / df, pt_hp / df
    A = al + y

    def lker(p):
        return (a - 1) * np.log(p) + (b + al - 1) * np.log1p(-p) - A * np.log(E + be * (k + 1) * (1 - p))

    # scale by the maximum of the kernel without the endpoint singularities
    res = optimize.minimize_scalar(lambda p: -(lker(p) - (a - 1) * np.log(p) - (b - 1) * np.log1p(-p)),
                                   bounds=(1e-12, 1 - 1e-12), method="bounded")
    m = lker(res.x)
    f = lambda p: np.exp(lker(p) - m)
    pc = k * rr0 / (k * rr0 + 1)
    pts = sorted(set([pc, res.x]))

    def I(g, lo=0.0, hi=1.0):
        br = [lo] + [q for q in pts if lo < q < hi] + [hi]
        return sum(integrate.quad(lambda p: g(p) * f(p), br[j], br[j + 1], **QK)[0] for j in range(len(br) - 1))

    Z = I(lambda p: 1.0)
    prob = I(lambda p: 1.0, pc, 1.0) / Z
    cdf = lambda p: I(lambda q: 1.0, 0.0, p) / Z
    qs = []
    for lev in (0.025, 0.5, 0.975):
        pq = optimize.brentq(lambda p: cdf(p) - lev, 1e-15, 1 - 1e-15, xtol=1e-15, rtol=1e-13)
        qs.append(pq / (k * (1 - pq)))
    B = lambda p: be + E / ((k + 1) * (1 - p))
    r = lambda p: p / (k * (1 - p))
    out = dict(prob=prob, r_median=qs[1], r_lower=qs[0], r_upper=qs[2],
               p_mean=I(lambda p: p) / Z, r_mean=I(r) / Z,
               lambda_p=I(lambda p: A / B(p)) / Z, lambda_d=I(lambda p: r(p) * A / B(p)) / Z)
    if t_star is not None:
        sp = I(lambda p: (B(p) / (B(p) + t_star)) ** A) / Z
        sd = I(lambda p: (B(p) / (B(p) + r(p) * t_star)) ** A) / Z
        out.update(ip_placebo=1 - sp, ip_drug=1 - sd, ip_difference=sp - sd)
    return out


def bdribs_2d_incidence(y, E, k, lam_hp, pt_hp, df, a, b, t_star):
    """Brute two-dimensional check over (delta, p) for the incidences."""
    al, be = lam_hp * pt_hp / df, pt_hp / df
    dmax = stats.gamma.ppf(1 - 1e-13, al + y, scale=1 / be) * 4

    def joint(d, p):
        lam = E * d / ((k + 1) * (1 - p))
        return stats.poisson.pmf(y, lam) * stats.gamma.pdf(d, al, scale=1 / be) * stats.beta.pdf(p, a, b)

    opts = dict(epsabs=0, epsrel=1e-10)
    Z = integrate.dblquad(joint, 0, 1, 0, dmax, **opts)[0]
    sp = integrate.dblquad(lambda d, p: np.exp(-d * t_star) * joint(d, p), 0, 1, 0, dmax, **opts)[0] / Z
    sd = integrate.dblquad(lambda d, p: np.exp(-p / (k * (1 - p)) * d * t_star) * joint(d, p), 0, 1, 0, dmax, **opts)[0] / Z
    return dict(ip_placebo=1 - sp, ip_drug=1 - sd)


cases = {
    "A": dict(y=5, E=500.0, k=2.0, lam_hp=18 / 4000, pt_hp=4000.0, df=1.0, a=0.5, b=0.5, rr0=1.0, t_star=1.0),
    "B": dict(y=30, E=8000.0, k=1.0, lam_hp=0.002, pt_hp=10000.0, df=2.0, a=0.5, b=0.5, rr0=1.2, t_star=52.0),
    "C": dict(y=0, E=300.0, k=1.0, lam_hp=0.01, pt_hp=500.0, df=1.0, a=1.0, b=1.0, rr0=1.0, t_star=10.0),
}
out = {}
for nm, c in cases.items():
    out[nm] = bdribs_ref(**c)
    print(nm, json.dumps(out[nm], indent=None))
chk = bdribs_2d_incidence(5, 500.0, 2.0, 18 / 4000, 4000.0, 1.0, 0.5, 0.5, 1.0)
print("A 2D check:", chk)
json.dump(out, open("ref_bdribs.json", "w"), indent=1)
