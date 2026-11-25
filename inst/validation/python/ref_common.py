"""Shared tools and data sets of the Python reference computations for the
models added in bssd 0.2.0 (ref_exptime.py, ref_mapbern.py, ref_canary.py,
ref_binbeta.py, ref_multibin.py, ref_misc.py).

The reference integrals do not use the algorithms of the package: they use
composite Gauss-Legendre rules on log or logit coordinates of both
parameters, split exactly at the threshold of the tested quantity (the split
point may differ from row to row), with fixed, dense panels; or closed-form
finite sums; or adaptive quadrature (scipy.integrate.quad).

Requires numpy and scipy. Run with python3 -B from this folder.
"""
import json
import numpy as np
from scipy import special
from scipy.optimize import brentq

GL_X, GL_W = np.polynomial.legendre.leggauss(8)


def unit_rule(panels, q=8):
    """Composite Gauss-Legendre rule on (0, 1)."""
    x, w = np.polynomial.legendre.leggauss(q)
    e = np.linspace(0.0, 1.0, panels + 1)
    h = np.diff(e) / 2
    m = (e[1:] + e[:-1]) / 2
    nodes = (h[:, None] * x[None, :] + m[:, None]).ravel()
    weights = (h[:, None] * w[None, :]).ravel()
    return nodes, weights


def rule(edges, q=8):
    """Composite Gauss-Legendre rule with the given panel edges."""
    x, w = np.polynomial.legendre.leggauss(q)
    e = np.asarray(edges, dtype=float)
    h = np.diff(e) / 2
    m = (e[1:] + e[:-1]) / 2
    nodes = (h[:, None] * x[None, :] + m[:, None]).ravel()
    weights = (h[:, None] * w[None, :]).ravel()
    return nodes, weights


def graded_edges(lo, hi, core_lo, core_hi, n_core, n_tail):
    """Panel edges dense on [core_lo, core_hi] and coarser outside."""
    core_lo = max(core_lo, lo)
    core_hi = min(core_hi, hi)
    parts = [np.linspace(core_lo, core_hi, n_core + 1)]
    if core_lo > lo:
        parts.insert(0, np.linspace(lo, core_lo, n_tail + 1)[:-1])
    if core_hi < hi:
        parts.append(np.linspace(core_hi, hi, n_tail + 1)[1:])
    return np.concatenate(parts)


class Split2D:
    """Posterior on (x, z) with a z-threshold that may depend on x.

    logpost(x, Z): x vector (nx), Z matrix (nx, nz) -> matrix of log
    densities. x_edges: panel edges in x. z_lo, z_hi: range of z; z_core:
    (lo, hi) of the dense part in z. z_hi_of: optional function of the x
    nodes giving an upper limit of z for each x (the density is zero above
    it), so that a density that ends at a limit depending on x is
    integrated exactly up to that limit.
    """

    def __init__(self, logpost, x_edges, z_lo, z_hi, z_core, n_core=120,
                 n_tail=40, q=8, z_hi_of=None):
        self.logpost = logpost
        self.x, self.wx = rule(x_edges, q)
        self.z_lo, self.z_hi = z_lo, z_hi
        self.z_core = z_core
        self.n_core, self.n_tail, self.q = n_core, n_tail, q
        if z_hi_of is None:
            self.z_top = None
            z, wz = rule(graded_edges(z_lo, z_hi, z_core[0], z_core[1],
                                      n_core, n_tail), q)
            Z = np.broadcast_to(z, (self.x.size, z.size))
            L = logpost(self.x, Z)
            self.lmax = np.max(L[np.isfinite(L)])
            D = np.exp(L - self.lmax) * self.wx[:, None] * wz[None, :]
            D[~np.isfinite(D)] = 0.0
            self.norm = D.sum()
            self.D = D / self.norm
            self.z = z
            return
        # one grid of z per x node, up to the limit of the row
        self.z_top = np.minimum(z_hi_of(self.x), z_hi)
        rows = []
        for i in range(self.x.size):
            zr, wr = self.row_rule(self.z_top[i])
            if zr.size == 0:
                rows.append((zr, wr, np.zeros(0)))
                continue
            L = self.logpost(self.x[i:i + 1], zr[None, :])[0]
            rows.append((zr, wr, L))
        self.lmax = max(np.max(L[np.isfinite(L)]) for _, _, L in rows
                        if L.size and np.any(np.isfinite(L)))
        dens = []
        for i, (zr, wr, L) in enumerate(rows):
            d = np.exp(L - self.lmax) * wr * self.wx[i]
            d[~np.isfinite(d)] = 0.0
            dens.append(d)
        self.norm = sum(d.sum() for d in dens)
        self.rows = [(zr, d / self.norm) for (zr, _, _), d in zip(rows, dens)]

    def row_rule(self, top, split=None):
        """Graded rule on (z_lo, min(top, split))."""
        hi = top if split is None else min(top, split)
        if hi <= self.z_lo:
            return np.zeros(0), np.zeros(0)
        e = graded_edges(self.z_lo, hi, min(self.z_core[0], hi),
                         min(self.z_core[1], hi), self.n_core, self.n_tail)
        e = e[np.concatenate(([True], np.diff(e) > 0))]
        if e.size < 2:
            return np.zeros(0), np.zeros(0)
        return rule(e, self.q)

    def below(self, s):
        """P(z <= s(x)); s is a vector over the x nodes."""
        s = np.clip(np.broadcast_to(s, self.x.shape), self.z_lo, self.z_hi)
        if self.z_top is not None:
            total = 0.0
            for i, si in enumerate(s):
                z, wz = self.row_rule(self.z_top[i], si)
                if z.size == 0:
                    continue
                L = self.logpost(self.x[i:i + 1], z[None, :])[0]
                v = np.exp(L - self.lmax) * wz
                v[~np.isfinite(v)] = 0.0
                total += self.wx[i] * v.sum()
            return total / self.norm
        total = 0.0
        for i, si in enumerate(s):
            if si <= self.z_lo:
                continue
            core_hi = min(self.z_core[1], si)
            core_lo = min(self.z_core[0], si)
            e = graded_edges(self.z_lo, si, core_lo, core_hi,
                             self.n_core, self.n_tail)
            e = e[np.concatenate(([True], np.diff(e) > 0))]
            if e.size < 2:
                continue
            z, wz = rule(e, self.q)
            L = self.logpost(self.x[i:i + 1], z[None, :])[0]
            v = np.exp(L - self.lmax) * wz
            v[~np.isfinite(v)] = 0.0
            total += self.wx[i] * v.sum()
        return total / self.norm

    def mean(self, f):
        """Posterior mean of f(x, Z)."""
        if self.z_top is not None:
            tot = 0.0
            with np.errstate(all="ignore"):
                for i, (zr, d) in enumerate(self.rows):
                    if zr.size:
                        tot += np.sum(d * f(np.full(zr.size, self.x[i]), zr))
            return float(tot)
        Z = np.broadcast_to(self.z, (self.x.size, self.z.size))
        with np.errstate(all="ignore"):
            F = f(self.x[:, None], Z)
        return float(np.sum(self.D * F))

    def quantile(self, s_of_c, p, lo, hi):
        """Value c with P(z <= s_of_c(c)) = p, by bisection on [lo, hi]."""
        return brentq(lambda c: self.below(s_of_c(c)) - p, lo, hi,
                      xtol=1e-11, rtol=1e-11)


def logaddexp(a, b):
    return np.logaddexp(a, b)


def make_data():
    """Data sets of the reference cases (fixed seeds)."""
    out = {}
    rs = np.random.RandomState(20261124)
    n = 60
    t = np.round(rs.uniform(3, 18, n), 1)
    arm = rs.binomial(1, 2 / 3, n)
    y = rs.poisson(np.where(arm == 1, 0.03, 0.02) * t)
    out["E1"] = {"y": y.tolist(), "t": t.tolist()}
    rs = np.random.RandomState(20261125)
    n = 80
    t = np.round(rs.uniform(1, 2, n), 2)
    arm = rs.binomial(1, 0.5, n)
    lam = np.where(arm == 1, 0.25, 0.15)
    tev = rs.exponential(1 / lam)
    y = (tev < t).astype(int)
    tr = np.where(y == 1, np.maximum(np.round(tev, 3), 0.001), t)
    out["E2"] = {"y": y.tolist(), "t": tr.tolist()}
    rs = np.random.RandomState(20261126)
    n = 60
    t = rs.choice([0.5, 1.0], n)
    arm = rs.binomial(1, 0.75, n)
    pi = np.where(arm == 1, 0.25, 0.12)
    y = rs.binomial(1, 1 - (1 - pi) ** t)
    out["M1"] = {"y": y.tolist(), "t": t.tolist()}
    # E3: data far in the upper tail of the informative prior of the study
    # drug rate of Schnell and Ball (2016)
    rs = np.random.RandomState(20261128)
    n = 120
    t = np.round(rs.uniform(6, 18, n), 1)
    arm = rs.binomial(1, 2 / 3, n)
    y = rs.poisson(np.where(arm == 1, 0.08, 0.013) * t)
    out["E3"] = {"y": y.tolist(), "t": t.tolist()}
    return out


def dump(obj, path):
    with open(path, "w", newline="\n") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
        f.write("\n")


def lbeta(a, b):
    return special.betaln(a, b)
