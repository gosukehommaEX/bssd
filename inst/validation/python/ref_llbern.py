"""Independent reference values for llbern_post.

All four parameters are integrated on their prior distribution-function scale
(the prior becomes uniform), with Gauss-Legendre nodes and no box refinement.
The indicator of RR(t*) > q is handled by splitting the mu_d coordinate at the
threshold for each (mu_p, sigma_p, sigma_d) node, with separate Gauss-Legendre
nodes on each piece.
"""
import json
import sys
import time
import numpy as np
from scipy import optimize, special, stats

rng = np.random.default_rng(7)
N = 60
t = np.round(rng.uniform(4, 104, N), 1)
grp = rng.random(N) < 0.5
mu_i = np.where(grp, 5.3, 6.0)
u = rng.random(N)
tev = np.exp(mu_i + 0.8 * np.log(u / (1 - u)))
y = (tev <= t).astype(int)


def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2


def make(k, mu_hp, sg_hp, pt_hp, df_p, t_star, md_pr=None, sd_pr=None):
    md_pr = (mu_hp, 1.0) if md_pr is None else md_pr
    sd_pr = (sg_hp, 1.0) if sd_pr is None else sd_pr
    G = [stats.gamma(pt_hp * mu_hp / df_p, scale=df_p / pt_hp), stats.gamma(pt_hp * sg_hp / df_p, scale=df_p / pt_hp),
         stats.gamma(md_pr[0], scale=1 / md_pr[1]), stats.gamma(sd_pr[0], scale=1 / sd_pr[1])]
    wp, wd = 1 / (k + 1), k / (k + 1)
    lt = np.log(t)
    lts = np.log(t_star)

    def engine(n, q_list, want_means=False):
        x, w = gl01(n)
        mp, sp, sd = G[0].ppf(x), G[1].ppf(x), np.maximum(G[3].ppf(x), 1e-300)
        MP, SP, SD = [a.ravel() for a in np.meshgrid(mp, sp, sd, indexing="ij")]
        WT = np.einsum("a,b,c->abc", w, w, w).ravel()
        zp = (lt[None, :] - MP[:, None]) / SP[:, None]
        Ap = np.where(y[None, :] == 1, wp * special.expit(zp), wp * special.expit(-zp))  # (T, N)
        Fps = special.expit((lts - MP) / SP)

        def piece_integrals(ua, ub, extra=None):
            # integrate likelihood over u_md in [ua, ub] for each triple with n GL nodes
            xx, ww = gl01(n)
            U = ua[:, None] + (ub - ua)[:, None] * xx[None, :]  # (T, n)
            Wn = (ub - ua)[:, None] * ww[None, :]
            MD = G[2].ppf(U)
            res = np.zeros(len(ua))
            mres = {}
            vals = np.empty(U.shape)
            for j in range(n):
                zd = (lt[None, :] - MD[:, j][:, None]) / SD[:, None]
                Bd = np.where(y[None, :] == 1, wd * special.expit(zd), wd * special.expit(-zd))
                vals[:, j] = np.log(Ap + Bd).sum(1)
            return MD, Wn, vals

        out = {}
        # full integral (no split) for normalization and means
        MD, Wn, vals = piece_integrals(np.zeros(len(MP)), np.ones(len(MP)))
        M = vals.max()
        L = np.exp(vals - M) * Wn * WT[:, None]
        Z = L.sum()
        if want_means:
            Fd = special.expit((lts - MD) / SD[:, None])
            Fp = Fps[:, None]
            out["mu_p"] = (L * MP[:, None]).sum() / Z
            out["sigma_p"] = (L * SP[:, None]).sum() / Z
            out["mu_d"] = (L * MD).sum() / Z
            out["sigma_d"] = (L * SD[:, None]).sum() / Z
            out["ip_placebo"] = (L * Fp).sum() / Z
            out["ip_drug"] = (L * Fd).sum() / Z
            out["ip_difference"] = out["ip_drug"] - out["ip_placebo"]
            out["rr_mean"] = (L * Fd / Fp).sum() / Z
        probs = []
        for q in q_list:
            arg = q * Fps
            ok = arg < 1
            mstar = np.where(ok, lts - SD * special.logit(np.minimum(arg, 1 - 1e-16)), -np.inf)
            ub = np.where(ok & (mstar > 0), G[2].cdf(np.maximum(mstar, 0)), 0.0)
            ub = np.clip(ub, 0, 1)
            MD1, Wn1, v1 = piece_integrals(np.zeros(len(MP)), ub)
            P1 = (np.exp(v1 - M) * Wn1 * WT[:, None]).sum()
            probs.append(P1 / Z)
        out["probs"] = probs
        return out

    return engine


if __name__ == "__main__":
    print("events", int(y.sum()), "of", N)
    res = dict(t=t.tolist(), y=y.tolist())
    eng = make(1.0, 6.0, 0.8, 200.0, 1.0, 104.0)
    for n in (48,):
        t0 = time.time()
        o = eng(n, [1.0, 2.0], want_means=True)
        print(n, {kk: (round(v, 7) if not isinstance(v, list) else [round(x, 7) for x in v]) for kk, v in o.items()},
              round(time.time() - t0, 1), "s", flush=True)
    res["S1"] = o
    n_q = 48
    qs = {}
    for nm, lev in (("rr_lower", 0.025), ("rr_median", 0.5), ("rr_upper", 0.975)):
        f = lambda lq: eng(n_q, [np.exp(lq)])["probs"][0] - (1 - lev)
        qs[nm] = float(np.exp(optimize.brentq(f, -690.0, 18.0, xtol=1e-7)))
        print(nm, qs[nm], flush=True)
    res["S1"].update(qs)
    eng2 = make(2.0, 6.0, 0.8, 200.0, 2.0, 52.0)
    for n in (40, 48):
        o2 = eng2(n, [1.5])
        print("S2", n, o2["probs"], flush=True)
    res["S2"] = dict(prob=o2["probs"][0])
    json.dump(res, open("ref_llbern.json", "w"), indent=1)
