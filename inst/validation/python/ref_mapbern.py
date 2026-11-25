"""Reference values for mapbern_post() and map_beta() (Brock et al., 2023).

mapbern_post(): coordinates (x, z) = (logit pi_c, logit pi_d), composite
Gauss-Legendre rules, split exactly at z = logit(min(c * pi_c, 1)).
map_beta() with summary data: the posterior of (logit a, logit b) on a dense
composite Gauss-Legendre grid (Beta-binomial likelihoods), and the MAP prior
as the resulting mixture; with subject-level data, pi_j is integrated out by
scipy.integrate.quad on (0, 1). Output: ref_mapbern.json.
"""
import numpy as np
from scipy import special, stats, integrate
import ref_common as rc

MIX1 = {"w": [0.7, 0.3], "shape1": [20, 2], "shape2": [180, 18]}
CASES = [
    {"id": "M1", "data": "M1", "k": 3, "prior_c": MIX1,
     "prior_d": [1 / 3, 1 / 3], "rr0": 1.0},
    {"id": "M2", "data": "M1", "k": 3, "prior_c": MIX1,
     "prior_d": [1 / 3, 1 / 3], "rr0": 1.5},
    # an informative prior of pi_d (mean 0.02) that conflicts with the data
    {"id": "M3", "data": "M1", "k": 3, "prior_c": MIX1,
     "prior_d": [2, 98], "rr0": 0.25},
]
MAP_CASES = [
    {"id": "P1", "events": [9, 11, 10], "n": [100, 100, 100]},
    {"id": "P2", "events": [3, 9, 1, 6], "n": [40, 60, 30, 80]},
]


def softplus(x):
    return np.logaddexp(0, x)


def make_logpost(y, t, k, mix, pd):
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    wc = 1 / (k + 1)
    w = np.asarray(mix["w"], float)
    s1 = np.asarray(mix["shape1"], float)
    s2 = np.asarray(mix["shape2"], float)

    def kern(h):
        ht = h[..., None] * t
        with np.errstate(divide="ignore"):
            return np.where(y > 0.5, np.log(-np.expm1(-ht)), -ht)

    def logpost(x, Z):
        X = np.broadcast_to(x[:, None], Z.shape)
        # hazard -log(1 - plogis(x)) = log(1 + exp(x))
        a = np.log(wc) + kern(softplus(X))
        b = np.log(1 - wc) + kern(softplus(Z))
        ll = np.logaddexp(a, b).sum(axis=-1)
        l1x, l0x = -softplus(-X), -softplus(X)
        comp = (s1 * l1x[..., None] + s2 * l0x[..., None] +
                np.log(w) - special.betaln(s1, s2))
        lpx = special.logsumexp(comp, axis=-1)
        lpz = pd[0] * (-softplus(-Z)) + pd[1] * (-softplus(Z)) - \
            special.betaln(pd[0], pd[1])
        return ll + lpx + lpz
    return logpost


def run_case(cs, data):
    d = data[cs["data"]]
    lp = make_logpost(d["y"], d["t"], cs["k"], cs["prior_c"], cs["prior_d"])
    xs = np.linspace(-12, 4, 161)
    zs = np.linspace(-20, 20, 401)
    L = lp(xs, np.broadcast_to(zs, (xs.size, zs.size)))
    keep = L > L.max() - 40
    xk = xs[np.any(keep, axis=1)]
    zk = zs[np.any(keep, axis=0)]
    S = rc.Split2D(lp, np.linspace(xk.min() - 0.5, xk.max() + 0.5, 41),
                   -150, 150, (zk.min() - 0.5, zk.max() + 0.5), n_core=120,
                   n_tail=60)
    sig = lambda v: 1 / (1 + np.exp(-v))
    split = lambda c: special.logit(np.minimum(c * sig(S.x), 1 - 1e-16))
    out = {"id": cs["id"]}
    out["prob_ratio"] = 1 - S.below(split(cs["rr0"]))
    out["mean_pi_c"] = S.mean(lambda x, z: sig(x))
    out["mean_pi_d"] = S.mean(lambda x, z: sig(z))
    out["q_pi_d"] = [S.quantile(lambda c: special.logit(c), p, 1e-12,
                                1 - 1e-12) for p in [0.025, 0.5, 0.975]]
    out["q_rr"] = [S.quantile(split, p, 1e-6, 1e3)
                   for p in [0.025, 0.5, 0.975]]
    return out


def map_summary(cs, a_prior=(1 / 3, 1 / 3), b_prior=(5, 5), n_side=160):
    ev = np.asarray(cs["events"], float)
    n = np.asarray(cs["n"], float)
    N = n.sum()

    def logpost(xa, xb):
        a = 1 / (1 + np.exp(-xa))
        b = 1 / (1 + np.exp(-xb))
        s1 = a * b * N
        s2 = (1 - a) * b * N
        out = (a_prior[0] * -softplus(-xa) + a_prior[1] * -softplus(xa) +
               b_prior[0] * -softplus(-xb) + b_prior[1] * -softplus(xb))
        for e, m in zip(ev, n):
            out = out + special.betaln(s1 + e, s2 + m - e) - \
                special.betaln(s1, s2)
        return out, s1, s2
    # locate the posterior on a coarse grid
    xa = np.linspace(-12, 12, 241)
    xb = np.linspace(-15, 15, 301)
    XA, XB = np.meshgrid(xa, xb, indexing="ij")
    L, _, _ = logpost(XA, XB)
    keep = L > L.max() - 46
    ra = xa[np.any(keep, axis=1)]
    rb = xb[np.any(keep, axis=0)]
    na, wa = rc.rule(np.linspace(ra.min() - 0.2, ra.max() + 0.2, n_side // 8
                                 + 1))
    nb, wb = rc.rule(np.linspace(rb.min() - 0.2, rb.max() + 0.2, n_side // 8
                                 + 1))
    XA, XB = np.meshgrid(na, nb, indexing="ij")
    L, S1, S2 = logpost(XA, XB)
    W = np.exp(L - L.max()) * wa[:, None] * wb[None, :]
    W = (W / W.sum()).ravel()
    S1 = S1.ravel()
    S2 = S2.ravel()
    mean = np.sum(W * S1 / (S1 + S2))
    m2 = np.sum(W * S1 * (S1 + 1) / ((S1 + S2) * (S1 + S2 + 1)))
    pts = [0.05, 0.08, 0.1, 0.12, 0.15]
    cdf = [float(np.sum(W * stats.beta.cdf(p, S1, S2))) for p in pts]
    A = (1 / (1 + np.exp(-XA))).ravel()
    B = (1 / (1 + np.exp(-XB))).ravel()
    return {"id": cs["id"], "mean": mean, "sd": np.sqrt(m2 - mean ** 2),
            "cdf_at": pts, "cdf": cdf, "mean_a": float(np.sum(W * A)),
            "mean_b": float(np.sum(W * B))}


def map_subject(y, t, study, a_prior=(1 / 3, 1 / 3), b_prior=(5, 5)):
    """MAP prior moments from subject-level data with quad inner integrals
    on a moderate (logit a, logit b) grid (independent of the package)."""
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    study = np.asarray(study)
    arms = [np.where(study == s)[0] for s in np.unique(study)]
    N = float(len(y))

    def arm_lik(idx):
        yy, tt = y[idx], t[idx]

        def f(p):
            with np.errstate(divide="ignore"):
                return np.sum(np.where(yy > 0.5,
                                       np.log(-np.expm1(tt * np.log1p(-p))),
                                       tt * np.log1p(-p)))
        return f
    fl = [arm_lik(i) for i in arms]
    # profile of each arm likelihood for scaling
    pg = np.linspace(1e-4, 0.6, 3000)
    lmax = [max(f(p) for p in pg) for f in fl]

    def log_marg(s1, s2):
        tot = 0.0
        for f, m in zip(fl, lmax):
            g = lambda p: np.exp(f(p) - m + (s1 - 1) * np.log(p) +
                                 (s2 - 1) * np.log1p(-p) -
                                 special.betaln(s1, s2))
            mode = (s1 - 1) / (s1 + s2 - 2) if s1 > 1 and s2 > 1 else 0.1
            v = integrate.quad(g, 0, 1, points=[mode], limit=400,
                               epsabs=0, epsrel=1e-11)[0]
            tot += np.log(v) + m
        return tot
    def lpost(u, v):
        a = 1 / (1 + np.exp(-u))
        b = 1 / (1 + np.exp(-v))
        return (log_marg(a * b * N, (1 - a) * b * N) +
                a_prior[0] * np.log(a) + a_prior[1] * np.log1p(-a) +
                b_prior[0] * np.log(b) + b_prior[1] * np.log1p(-b))
    # coarse scan to locate the posterior of (logit a, logit b)
    ua = np.linspace(-6, 0, 31)
    ub = np.linspace(-8, 8, 33)
    C = np.array([[lpost(u, v) for v in ub] for u in ua])
    keep = C > C.max() - 30
    ra = ua[np.any(keep, axis=1)]
    rb = ub[np.any(keep, axis=0)]
    xan, wa = rc.rule(np.linspace(ra.min() - 0.4, ra.max() + 0.4, 9), q=6)
    xbn, wb = rc.rule(np.linspace(rb.min() - 1, rb.max() + 1, 9), q=6)
    res = np.array([[lpost(u, v) for v in xbn] for u in xan])
    W = np.exp(res - res.max()) * wa[:, None] * wb[None, :]
    W = W / W.sum()
    A = 1 / (1 + np.exp(-xan))[:, None] * np.ones_like(W)
    B = np.ones_like(W) / (1 + np.exp(-xbn))[None, :]
    S1 = A * B * N
    S2 = (1 - A) * B * N
    mean = np.sum(W * S1 / (S1 + S2))
    m2 = np.sum(W * S1 * (S1 + 1) / ((S1 + S2) * (S1 + S2 + 1)))
    edge = float(W[[0, -1], :].sum() + W[:, [0, -1]].sum())
    return {"mean": mean, "sd": float(np.sqrt(m2 - mean ** 2)),
            "edge_mass": edge}


if __name__ == "__main__":
    data = rc.make_data()
    res = {"cases": CASES, "results": [run_case(c, data) for c in CASES],
           "map_cases": MAP_CASES,
           "map_results": [map_summary(c) for c in MAP_CASES]}
    # subject-level data: arms with mixed follow-up
    rs = np.random.RandomState(20261127)
    sizes = [40, 50, 60]
    rates = [0.10, 0.08, 0.12]
    y, t, st = [], [], []
    for j, (m, p) in enumerate(zip(sizes, rates)):
        tt = rs.choice([0.5, 1.0, 1.5], m)
        y += rs.binomial(1, 1 - (1 - p) ** tt).tolist()
        t += tt.tolist()
        st += [j + 1] * m
    res["map_subject_data"] = {"y": y, "t": t, "study": st}
    res["map_subject"] = map_subject(y, t, st)
    rc.dump(res, "ref_mapbern.json")
    for r in res["results"] + res["map_results"]:
        print(r)
    print(res["map_subject"])
