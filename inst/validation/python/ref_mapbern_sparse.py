"""Reference values for mapbern_post() with MAP priors from map_beta() for
historical control arms with few events (cases MS1 to MS7). The MAP prior
of Brock et al. (2023) is pi_j ~ Beta(a b N, (1 - a) b N), N = sum N_j,
a ~ Beta(1/3, 1/3), b ~ Beta(5, 5); the trial has complete follow-up
(t = 1), randomization 3:1, the prior Beta(1/3, 1/3) of pi_d and Y events
among n subjects, so that the blinded likelihood is q^Y (1 - q)^(n - Y)
with q = (pi_c + 3 pi_d) / 4.

The method differs from the package: the posterior of (logit a, logit b)
is integrated by composite Gauss-Legendre (GL) rules on the box where its
log density is within 45 of the maximum (located on a grid of step 0.1 and
0.05 over [-110, 20] x [-14, 14]), instead of the mixture of map_beta().
For each node, pi_c ~ Beta(s1, s2) is integrated on
z = qnorm(I(pi_c; s1, s2)) by composite GL rules on [-9, 9]. pi_d is
integrated on s = pi_d^(1/3) (pi_d <= 1/2) and t = (1 - pi_d)^(1/3)
(pi_d > 1/2), where the prior density times the Jacobian is smooth; for
P(pi_d > x) the s or t interval starts exactly at the threshold x
(x = rr0 * pi_c or x = rate0). Each case is computed twice, the second time
with about 1.5 times as many nodes in each coordinate; the second result
and the largest change are written to ref_mapbern_sparse.json.

Run as python3 ref_mapbern_sparse.py [case ...]; cases already in the JSON
file are skipped. About 1.5 hours for all cases (10 to 16 minutes per
case, measured with two cases running at the same time on two cores).
"""
import numpy as np
from scipy import special

S = 0.5 ** (1 / 3)
LB13 = special.betaln(1 / 3, 1 / 3)


def leg(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * x + (a + b) / 2, w * (b - a) / 2


def comp(edges, n):
    xs, ws = zip(*[leg(n, edges[j], edges[j + 1]) for j in range(len(edges) - 1)])
    return np.concatenate(xs), np.concatenate(ws)


def ab_logpost(LA, LB, events, sizes):
    N = sizes.sum()
    A, B = special.expit(LA), special.expit(LB)
    s1 = A * B * N
    s2 = (1 - A) * B * N
    lp = ((1 / 3) * special.log_expit(LA) + (1 / 3) * special.log_expit(-LA)
          + 5 * special.log_expit(LB) + 5 * special.log_expit(-LB))
    ll = -len(sizes) * special.betaln(s1, s2)
    for xj, Nj in zip(events, sizes):
        ll = ll + special.betaln(xj + s1, Nj - xj + s2)
    return lp + ll


def ab_box(events, sizes, drop=45.0):
    """box of (logit a, logit b) where the log posterior density is within
    'drop' of its maximum, located on a uniform grid over [-110, 20] x
    [-14, 14] (step 0.1 and 0.05) and widened by two grid steps."""
    ga = np.linspace(-110, 20, 1301)
    gb = np.linspace(-14, 14, 561)
    L = ab_logpost(ga[:, None], gb[None, :], events, sizes)
    mask = L > L.max() - drop
    ia = np.where(mask.any(1))[0]
    ib = np.where(mask.any(0))[0]
    return (ga[max(ia.min() - 2, 0)], ga[min(ia.max() + 2, len(ga) - 1)],
            gb[max(ib.min() - 2, 0)], gb[min(ib.max() + 2, len(gb) - 1)])


def ab_posterior(events, sizes, na, nb):
    N = sizes.sum()
    a_lo, a_hi, b_lo, b_hi = ab_box(events, sizes)
    ua, wa = comp(np.linspace(a_lo, a_hi, na[0] + 1), na[1])
    ub, wb = comp(np.linspace(b_lo, b_hi, nb[0] + 1), nb[1])
    LA, LB = np.meshgrid(ua, ub, indexing="ij")
    A, B = special.expit(LA), special.expit(LB)
    s1 = A * B * N
    s2 = (1 - A) * B * N
    lp = ((1 / 3) * special.log_expit(LA) + (1 / 3) * special.log_expit(-LA)
          + 5 * special.log_expit(LB) + 5 * special.log_expit(-LB))
    ll = -len(sizes) * special.betaln(s1, s2)
    for xj, Nj in zip(events, sizes):
        ll = ll + special.betaln(xj + s1, Nj - xj + s2)
    lw = (lp + ll + np.log(np.outer(wa, wb))).ravel()
    w = np.exp(lw - lw.max())
    w /= w.sum()
    keep = w > 1e-15
    return w[keep] / w[keep].sum(), s1.ravel()[keep], s2.ravel()[keep]


def pc_of_z(z, s1, s2):
    """pi_c = I^{-1}(Phi(z); s1, s2) from the smaller tail; arrays (K, nz)."""
    Z = np.broadcast_to(z, (len(s1), len(z)))
    A1 = np.broadcast_to(s1[:, None], Z.shape)
    A2 = np.broadcast_to(s2[:, None], Z.shape)
    out = np.empty(Z.shape)
    lo = Z <= 0
    out[lo] = special.betaincinv(A1[lo], A2[lo], special.ndtr(Z[lo]))
    out[~lo] = 1 - special.betaincinv(A2[~lo], A1[~lo], special.ndtr(-Z[~lo]))
    bad = np.isnan(out)
    if bad.any():
        # betaincinv of SciPy returns NaN for a few arguments: bisection on
        # betainc instead (on the smaller tail)
        a1, a2, zz = A1[bad], A2[bad], Z[bad]
        target = special.ndtr(zz)
        lo_x, hi_x = np.zeros(len(zz)), np.ones(len(zz))
        for _ in range(200):
            mid = (lo_x + hi_x) / 2
            f = special.betainc(a1, a2, mid) - target
            lo_x = np.where(f < 0, mid, lo_x)
            hi_x = np.where(f < 0, hi_x, mid)
        out[bad] = (lo_x + hi_x) / 2
    return out


def reference(events, sizes, Y, n, k=3.0, rr0s=(1.0, 2.0), rate0=0.05,
              na=(52, 10), nb=(14, 10), nz=(12, 12), ns=(24, 16), npart=(4, 16),
              chunk=150):
    wc = 1 / (k + 1)
    w, s1, s2 = ab_posterior(np.asarray(events), np.asarray(sizes), na, nb)
    z, wz = comp(np.linspace(-9, 9, nz[0] + 1), nz[1])
    wz = wz * np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi)   # prior density of z is N(0, 1)
    sg, ws = comp(np.linspace(0, S, ns[0] + 1), ns[1])
    pd_lo, pd_hi = sg ** 3, 1 - sg ** 3
    dens = 3 * (1 - sg ** 3) ** (-2 / 3) / np.exp(LB13) * ws     # weight per node (both halves)
    q_hat = min(max(Y / n, 1e-12), 1 - 1e-12) if n > 0 else 0.5
    lmax = Y * np.log(q_hat) + (n - Y) * np.log1p(-q_hat)

    def lik(pc, pdv):
        q = wc * pc + (1 - wc) * pdv
        with np.errstate(divide="ignore"):
            return np.exp(Y * np.log(q) + (n - Y) * np.log1p(-q) - lmax)

    xg, xw = comp(np.linspace(0, 1, npart[0] + 1), npart[1])   # on [0, 1]

    def above(PC, x):
        """integral of p(pi_d) L over pi_d > x (x has the shape of PC)."""
        x = np.clip(x, 0, 1)
        lowx = x <= 0.5
        s0 = np.where(lowx, np.cbrt(x), S)                       # lower half from s0 to S
        t1 = np.where(lowx, S, np.cbrt(np.maximum(1 - x, 0)))    # upper half from 0 to t1
        sv = s0[..., None] + (S - s0)[..., None] * xg
        swt = (S - s0)[..., None] * xw * 3 * (1 - sv ** 3) ** (-2 / 3) / np.exp(LB13)
        tv = t1[..., None] * xg
        twt = t1[..., None] * xw * 3 * (1 - tv ** 3) ** (-2 / 3) / np.exp(LB13)
        P = PC[..., None]
        return (lik(P, sv ** 3) * swt).sum(-1) + (lik(P, 1 - tv ** 3) * twt).sum(-1)

    tot = 0.0
    acc = {f"prob_ratio_{r:g}": 0.0 for r in rr0s}
    acc.update(prob_rate=0.0, mean_pc=0.0, mean_pd=0.0)
    for c0 in range(0, len(w), chunk):
        wk, a1, a2 = w[c0:c0 + chunk], s1[c0:c0 + chunk], s2[c0:c0 + chunk]
        PC = pc_of_z(z, a1, a2)
        P = PC[..., None]
        L1 = lik(P, pd_lo) * dens
        L2 = lik(P, pd_hi) * dens
        g = L1.sum(-1) + L2.sum(-1)
        gd = (L1 * pd_lo).sum(-1) + (L2 * pd_hi).sum(-1)
        W = wk[:, None] * wz[None, :]
        tot += (W * g).sum()
        acc["mean_pc"] += (W * g * PC).sum()
        acc["mean_pd"] += (W * gd).sum()
        for r in rr0s:
            acc[f"prob_ratio_{r:g}"] += (W * above(PC, r * PC)).sum()
        acc["prob_rate"] += (W * above(PC, np.full(PC.shape, rate0))).sum()
    return {key: val / tot for key, val in acc.items()}, len(w)


CASES = {
    "MS1": dict(events=(1, 1, 0), sizes=(40, 40, 40), Y=1, n=50),
    "MS2": dict(events=(0, 0, 0), sizes=(40, 40, 40), Y=1, n=50),
    "MS3": dict(events=(2, 1, 1), sizes=(40, 40, 40), Y=2, n=100),
    "MS4": dict(events=(0, 0, 0), sizes=(40, 40, 40), Y=10, n=50),
    "MS5": dict(events=(2, 5, 7, 5, 10, 14, 10, 19, 28),
                sizes=(50, 50, 50, 100, 100, 100, 200, 200, 200), Y=95, n=1000),
    "MS6": dict(events=(1, 2, 3), sizes=(100, 100, 100), Y=3, n=200),
    "MS7": dict(events=(9, 11, 10), sizes=(100, 100, 100), Y=30, n=200),
}


if __name__ == "__main__":
    import json
    import os
    import sys
    import time
    path = "ref_mapbern_sparse.json"
    out = json.load(open(path)) if os.path.exists(path) else {
        "note": "mapbern_post() with MAP priors from few historical events; "
                "k = 3, prior_d = Beta(1/3, 1/3), t = 1, rr0 = 1 and 2, rate0 = 0.05",
        "results": []}
    done = {r["id"] for r in out["results"]}
    for nm in (sys.argv[1:] or list(CASES)):
        if nm in done:
            continue
        c = CASES[nm]
        t0 = time.time()
        r1, k1 = reference(c["events"], np.array(c["sizes"]), c["Y"], c["n"])
        r2, k2 = reference(c["events"], np.array(c["sizes"]), c["Y"], c["n"],
                           na=(78, 12), nb=(20, 12), nz=(18, 14), ns=(36, 18),
                           npart=(6, 18))
        res = dict(id=nm, events=list(c["events"]), sizes=list(c["sizes"]),
                   Y=c["Y"], n=c["n"], **r2,
                   change=max(abs(r2[key] - r1[key]) for key in r1),
                   components=[k1, k2], seconds=round(time.time() - t0))
        print(json.dumps(res), flush=True)
        # re-read, so that runs of different cases at the same time add up
        out = json.load(open(path)) if os.path.exists(path) else out
        out["results"] = [r for r in out["results"] if r["id"] != nm] + [res]
        out["results"].sort(key=lambda r: r["id"])
        with open(path, "w", newline="\n") as f:
            json.dump(out, f, indent=1, sort_keys=True)
            f.write("\n")
