"""Reference values for expbern_post when the data conflict with an
informative prior of p (cases X1 to X3).

The posterior is integrated on the coordinates u = log(lambda_p) and
z = logit(p) = log(k r), not on the prior distribution function of p used by
bssd, and no Chebyshev interpolation is used:

- u: composite Gauss-Legendre rules on the range where the maximum of the
  log posterior density over z is within `drop` of the overall maximum
  (located on a grid of 400 points between the prior quantiles 1e-14 and
  1 - 1e-14 of lambda_p and a grid of z of spacing 0.1 over [-40, 40]);
- z: for every node of u, the range where the log density is within `drop`
  of the maximum of that row is located on a grid of spacing 0.05 over
  [-40, 40], and z is integrated by composite Gauss-Legendre rules on that
  range, split at the threshold (probabilities and quantiles).

Each case is computed with two numbers of nodes; the output keeps the finer
result and the largest change between the two.

Data (the same as in inst/validation/validate_methods_python.R and
tests/testthat/test-expbern_post.R): t = 4, 5, ..., 52 repeated 8 times
(392 subjects) and y_i = 1 when (37 i) mod 100 is below 15 (data set X1,
59 events, about 2.7 times the events expected from the historical placebo
rate) or below 3 (data set X3, 11 events). Priors: lambda_hp =
-log(0.9) / 52, pt_hp = 392 / 2 * 52, df = 2; X1 and X3 use Beta(10, 10) for
p (k = 1 and 2), X2 uses Beta(30, 30) and rr0 = 1.5.

Run from this folder: python3 ref_expbern_rows.py (about 15 minutes).
Output: ref_expbern_rows.json.
"""
import json
import numpy as np
from scipy import optimize, special, stats

LAM_HP = -np.log(0.9) / 52


def make_data(thr):
    t = np.tile(np.arange(4, 53), 8).astype(float)
    i = np.arange(1, len(t) + 1)
    return ((i * 37) % 100 < thr).astype(int), t


def leg(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * x + (a + b) / 2, w * (b - a) / 2


def composite(a, b, panels, n):
    e = np.linspace(a, b, panels + 1)
    xs, ws = zip(*[leg(n, e[j], e[j + 1]) for j in range(panels)])
    return np.concatenate(xs), np.concatenate(ws)


def loglik(t, y, wp, lam, r):
    """log-likelihood on the product grid lam x r."""
    wd = 1 - wp
    out = np.empty((len(lam), len(r)))
    y1 = y == 1
    for j, la in enumerate(lam):
        h1 = la * t
        h2 = np.outer(r * la, t)
        l1 = np.log(wp * (-np.expm1(-h1))[None, :] + wd * (-np.expm1(-h2)))
        l0 = np.logaddexp(np.log(wp) - h1[None, :], np.log(wd) - h2)
        out[j] = np.where(y1[None, :], l1, l0).sum(1)
    return out


class Posterior:
    def __init__(self, y, t, k, lambda_hp, pt_hp, df, p_prior, n_out, n_in,
                 drop=60.0):
        self.y, self.t, self.k = y, t, k
        self.G = stats.gamma(lambda_hp * pt_hp / df, scale=df / pt_hp)
        self.a, self.b = p_prior
        self.wp = 1 / (k + 1)
        self.n_in = n_in
        u_lo = np.log(self.G.ppf(1e-14))
        u_hi = np.log(self.G.ppf(1 - 1e-14))
        uc = np.linspace(u_lo, u_hi, 400)
        Lc = self.logjoint(uc, np.arange(-40, 40.0001, 0.1))
        rm = Lc.max(1)
        keep = np.where(rm > rm.max() - drop)[0]
        du = uc[1] - uc[0]
        self.u, self.wu = composite(max(uc[keep[0]] - du, u_lo),
                                    min(uc[keep[-1]] + du, u_hi), *n_out)
        zf = np.arange(-40, 40.0001, 0.05)
        Lf = self.logjoint(self.u, zf)
        rmf = Lf.max(1)
        self.M = rmf.max()
        self.zl = np.empty(len(self.u))
        self.zh = np.empty(len(self.u))
        for i in range(len(self.u)):
            ii = np.where(Lf[i] > rmf[i] - drop)[0]
            self.zl[i] = max(zf[ii[0]] - 0.05, -40)
            self.zh[i] = min(zf[ii[-1]] + 0.05, 40)
        self.Z = self.below(np.inf)

    def logjoint(self, u, z):
        u = np.atleast_1d(u)
        z = np.atleast_1d(z)
        ll = loglik(self.t, self.y, self.wp, np.exp(u), np.exp(z) / self.k)
        lpu = self.G.logpdf(np.exp(u)) + u
        lpz = (self.a * special.log_expit(z) + self.b * special.log_expit(-z)
               - special.betaln(self.a, self.b))
        return ll + lpu[:, None] + lpz[None, :]

    def below(self, zq):
        """integral of the joint density over z < zq (scaled by exp(-M))."""
        tot = 0.0
        for i in range(len(self.u)):
            hi = min(zq, self.zh[i])
            if hi <= self.zl[i]:
                continue
            z, wz = composite(self.zl[i], hi, *self.n_in)
            L = self.logjoint(self.u[i], z)[0]
            tot += self.wu[i] * (np.exp(L - self.M) * wz).sum()
        return tot

    def prob(self, rr0):
        return 1 - self.below(np.log(self.k * rr0)) / self.Z

    def quantile(self, p):
        f = lambda zq: self.below(zq) / self.Z - p
        zq = optimize.brentq(f, self.zl.min(), self.zh.max(), xtol=1e-11)
        return float(np.exp(zq) / self.k)

    def moments(self, t_star):
        acc = np.zeros(7)
        for i in range(len(self.u)):
            z, wz = composite(self.zl[i], self.zh[i], *self.n_in)
            L = self.logjoint(self.u[i], z)[0]
            w = self.wu[i] * np.exp(L - self.M) * wz
            lam = np.exp(self.u[i])
            r = np.exp(z) / self.k
            fp = -np.expm1(-lam * t_star)
            fd = -np.expm1(-r * lam * t_star)
            acc += np.array([w.sum(), (w * special.expit(z)).sum(),
                             (w * r).sum(), (w * lam).sum(),
                             (w * r * lam).sum(), (w * fp).sum(),
                             (w * fd).sum()])
        m = acc[1:] / acc[0]
        return dict(mean_p=m[0], mean_r=m[1], mean_lambda_p=m[2],
                    mean_lambda_d=m[3], ip_placebo=m[4], ip_drug=m[5],
                    ip_difference=m[5] - m[4])


CASES = (
    ("X1", 15, dict(k=1.0, p_prior=(10.0, 10.0)), 1.0),
    ("X2", 15, dict(k=1.0, p_prior=(30.0, 30.0)), 1.5),
    ("X3", 3, dict(k=2.0, p_prior=(10.0, 10.0)), 1.0),
)


def compute(thr, s, rr0, n_out, n_in):
    y, t = make_data(thr)
    post = Posterior(y, t, lambda_hp=LAM_HP, pt_hp=len(t) / 2 * 52, df=2.0,
                     n_out=n_out, n_in=n_in, **s)
    out = dict(prob=post.prob(rr0),
               q=[post.quantile(p) for p in (0.025, 0.5, 0.975)])
    out.update(post.moments(52.0))
    return out


if __name__ == "__main__":
    res = dict(results=[], data={})
    for name, thr in (("X1", 15), ("X3", 3)):
        y, t = make_data(thr)
        res["data"][name] = dict(y=y.tolist(), t=t.tolist())
    for cid, thr, s, rr0 in CASES:
        coarse = compute(thr, s, rr0, (16, 32), (8, 32))
        fine = compute(thr, s, rr0, (24, 48), (12, 48))
        change = max(abs(fine["prob"] - coarse["prob"]),
                     max(abs(x / c - 1)
                         for x, c in zip(fine["q"], coarse["q"])),
                     max(abs(fine[key] / coarse[key] - 1)
                         for key in fine if key.startswith(("mean", "ip_"))))
        row = dict(id=cid, dataset="X1" if thr == 15 else "X3", rr0=rr0,
                   change=change, **fine)
        print(json.dumps(row), flush=True)
        res["results"].append(row)
    with open("ref_expbern_rows.json", "w") as fh:
        json.dump(res, fh, indent=1)
