"""Independent simulation of the operating characteristics computed by sim_oc
in reproduce_Duncan2026_simulation.R, for validate_sim_oc.R.

The trials are simulated with NumPy from the design of sim_blinded (uniform
accrual, round(n * k / (k + 1)) subjects on study drug chosen at random,
exponential or log-logistic times to the event, follow-up time
min(follow_up, tau - entry) at the analysis time tau, subjects not yet
entered excluded). The posterior quantities are computed without the code
of the package:

- BDRIBS: the Gamma prior of the placebo rate is integrated out analytically,
  and the density of z = logit(p),
  p^a (1 - p)^(b + alpha) [E + beta (k + 1) (1 - p)]^(-(alpha + Y)),
  is integrated with the composite Simpson rule on a uniform grid of
  z in [-80, 80] with step 0.005 that has logit(p0) as a node
  (p0 = k rr0 / (1 + k rr0));
- EXP_BERN: the posterior of (log(lambda_p), logit(p)) is evaluated on a
  coarse uniform grid to find the region that holds the posterior mass, and
  then integrated with the composite Simpson rule on a uniform grid of
  61 x 801 points on that region (split at logit(p0)).

The check mode (python3 ref_sim_oc.py check) compares these functions with
the reference values of tests/testthat/test-bdribs_post.R (cases A and B)
and test-expbern_post.R (the first two settings); the differences were below
3e-12 (BDRIBS) and 5e-6 (EXP_BERN).

Scenarios (time in weeks; the settings of reproduce_Duncan2026_simulation.R):
Figures 3 and 4 (BDRIBS(1) and BDRIBS(2)), Figures 5 and 9 (BDRIBS(2)), 2000
simulated trials each, and scenario 5.2 of Figure 8 (EXP_BERN, looks 1 to 3),
1000 simulated trials.

Output: ref_sim_oc.csv. Run from this folder: python3 ref_sim_oc.py
(about 30 minutes).
"""
import csv
import os
import sys
import numpy as np
from scipy import special, stats


def log_expit(z):
    return -np.logaddexp(0.0, -z)


def simpson_weights(n, h):
    """Composite Simpson weights for n (odd) equally spaced points."""
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return w * h / 3.0


def split_grid(lo, hi, z0, n_half):
    """Simpson nodes and weights on [lo, z0] and [z0, hi] (z0 shared)."""
    za = np.linspace(lo, z0, n_half)
    zb = np.linspace(z0, hi, n_half)
    wa = simpson_weights(n_half, (z0 - lo) / (n_half - 1))
    wb = simpson_weights(n_half, (hi - z0) / (n_half - 1))
    z = np.concatenate([za, zb[1:]])
    w = np.concatenate([wa[:-1], [wa[-1] + wb[0]], wb[1:]])
    # the shared node z0 is counted with its weight from [z0, hi] only
    w_up = np.concatenate([np.zeros(n_half - 1), [wb[0]], wb[1:]])
    return z, w, w_up


def bdribs_prob(y, e, k, lambda_hp, pt_hp, df, a=0.5, b=0.5, rr0=1.0,
                step=0.005, width=80.0, chunk=200):
    """P(r > rr0 | Y, E) for arrays of Y and E."""
    alpha = lambda_hp * pt_hp / df
    beta = pt_hp / df
    z0 = special.logit(k * rr0 / (1 + k * rr0))
    n_half = int(round(width / step)) + 1
    z, w, w_up = split_grid(z0 - width, z0 + width, z0, n_half)
    lp, lq = log_expit(z), log_expit(-z)
    y = np.asarray(y, float)
    e = np.asarray(e, float)
    out = np.empty(len(y))
    for s in range(0, len(y), chunk):
        ys, es = y[s:s + chunk, None], e[s:s + chunk, None]
        ld = (a * lp + (b + alpha) * lq)[None, :] - (alpha + ys) * np.log(
            es + beta * (k + 1) * np.exp(lq)[None, :])
        f = np.exp(ld - ld.max(1, keepdims=True))
        out[s:s + chunk] = (f * w_up).sum(1) / (f * w).sum(1)
    return out


def expbern_fit(y, t, k, lambda_hp, pt_hp, df, a=0.5, b=0.5, rr0=1.0,
                t_star=52.0, n_u=61, n_z_half=401):
    """P(r > rr0 | data) and the posterior mean of F_d(t*) - F_p(t*)."""
    alpha = lambda_hp * pt_hp / df
    beta = pt_hp / df
    wp, wd = 1 / (k + 1), k / (k + 1)
    z0 = special.logit(k * rr0 / (1 + k * rr0))
    y = np.asarray(y, bool)
    t1, t0 = np.asarray(t, float)[y], np.asarray(t, float)[~y]
    g = stats.gamma(alpha, scale=1 / beta)

    def logpost(u, z):
        lam = np.exp(u)
        r = np.exp(z) / k
        lam_d = r * lam
        ll = np.zeros(lam.shape)
        if len(t1):
            ll += np.log(-(wp * np.expm1(-lam[..., None] * t1) +
                           wd * np.expm1(-lam_d[..., None] * t1))).sum(-1)
        if len(t0):
            ll += np.log(wp * np.exp(-lam[..., None] * t0) +
                         wd * np.exp(-lam_d[..., None] * t0)).sum(-1)
        # Gamma prior of lambda on the log scale, Beta prior of p on logit
        return ll + alpha * u - beta * lam + a * log_expit(z) + b * log_expit(-z)

    # coarse grid to locate the posterior mass
    u_lo, u_hi = np.log(g.ppf(1e-13)), np.log(g.ppf(1 - 1e-13))
    uc = np.linspace(u_lo - 1.0, u_hi + 1.0, 61)
    zc = np.linspace(z0 - 80, z0 + 80, 481)
    U, Z = np.meshgrid(uc, zc, indexing="ij")
    L = logpost(U, Z)
    m = np.exp(L - L.max())
    mu, mz = m.sum(1), m.sum(0)
    iu = np.where(mu > 1e-14 * mu.max())[0]
    iz = np.where(mz > 1e-14 * mz.max())[0]
    hu, hz = uc[1] - uc[0], zc[1] - zc[0]
    u_a, u_b = uc[iu[0]] - 2 * hu, uc[iu[-1]] + 2 * hu
    z_a, z_b = zc[iz[0]] - 2 * hz, zc[iz[-1]] + 2 * hz
    if iz[0] == 0 or iz[-1] == len(zc) - 1 or iu[0] == 0 or iu[-1] == len(uc) - 1:
        raise RuntimeError("posterior mass at the edge of the coarse grid")
    # fine grid (Simpson), split at z0
    uf = np.linspace(u_a, u_b, n_u)
    wu = simpson_weights(n_u, (u_b - u_a) / (n_u - 1))
    z_a, z_b = min(z_a, z0 - hz), max(z_b, z0 + hz)
    zf, wz, wz_up = split_grid(z_a, z_b, z0, n_z_half)
    U, Z = np.meshgrid(uf, zf, indexing="ij")
    L = logpost(U, Z)
    f = np.exp(L - L.max()) * wu[:, None]
    lam = np.exp(U)
    lam_d = np.exp(Z) / k * lam
    rd = np.exp(-lam * t_star) - np.exp(-lam_d * t_star)
    tot = (f * wz[None, :]).sum()
    prob = (f * wz_up[None, :]).sum() / tot
    mean_rd = (f * rd * wz[None, :]).sum() / tot
    return prob, mean_rd


def simulate(rng, n, k, accrual, follow_up, taus, dist, par_p, par_d):
    n_d = int(round(n * k / (k + 1)))
    drug = np.zeros(n, bool)
    drug[rng.permutation(n)[:n_d]] = True
    entry = rng.uniform(0, accrual, n)
    if dist == "exponential":
        tev = rng.exponential(1 / np.where(drug, par_d, par_p))
    else:
        mu = np.where(drug, par_d[0], par_p[0])
        sg = np.where(drug, par_d[1], par_p[1])
        u = rng.uniform(size=n)
        tev = np.exp(mu + sg * np.log(u / (1 - u)))
    looks = []
    for tau in taus:
        fu = np.minimum(follow_up, tau - entry)
        keep = fu > 0
        ev = tev <= fu
        tar = np.where(ev, tev, fu)
        looks.append(dict(fu=fu[keep], ev=ev[keep], tar=tar[keep]))
    return looks


def ll_cdf(t, mu, sigma):
    return special.expit((np.log(t) - mu) / sigma)


def scenarios():
    rate_of = lambda ip, t: -np.log(1 - ip) / t
    sc = []
    for i, ip_d in enumerate([0.10, 0.15, 0.20, 0.25], 1):
        sc.append(dict(figure="3", point=i, n=350, follow_up=52,
                       taus=[4 * 52 / 3], dist="exponential",
                       par_p=rate_of(0.10, 52), par_d=rate_of(ip_d, 52),
                       ip_hp=0.10, models=["BDRIBS1", "BDRIBS2"], cutoff=0.925))
    for i, ip_d in enumerate([0.30, 0.40, 0.45, 0.50], 1):
        sc.append(dict(figure="4", point=i, n=200, follow_up=52,
                       taus=[26 + 52], dist="exponential",
                       par_p=rate_of(0.30, 52), par_d=rate_of(ip_d, 52),
                       ip_hp=0.30, models=["BDRIBS1", "BDRIBS2"], cutoff=0.95))
    for i, mu_d in enumerate([9.85, 8.35, 7.4], 1):
        sc.append(dict(figure="5", point=i, n=375, follow_up=52,
                       taus=[3 * (52 / 3)], dist="loglogistic",
                       par_p=(9.85, 2.0), par_d=(mu_d, 2.0),
                       ip_hp=ll_cdf(52, 9.85, 2.0), models=["BDRIBS2"],
                       cutoff=0.95))
    for i, mu_d in enumerate([7.3, 6.5, 6.1, 5.85], 1):
        sc.append(dict(figure="9", point=i, n=150, follow_up=104,
                       taus=[3 * 39], dist="loglogistic",
                       par_p=(7.3, 0.8), par_d=(mu_d, 0.8),
                       ip_hp=ll_cdf(104, 7.3, 0.8), models=["BDRIBS2"],
                       cutoff=0.95))
    f52 = ll_cdf(52, 6.17, 1.6)
    sc.append(dict(figure="8", point=1, n=400, follow_up=52,
                   taus=[52 / 3, 2 * 52 / 3, 3 * 52 / 3], dist="exponential",
                   par_p=rate_of(f52, 52), par_d=rate_of(f52, 52),
                   ip_hp=f52, models=["EXP_BERN"], cutoff=0.95))
    for s in sc:
        s["lambda_hp"] = rate_of(s["ip_hp"], s["follow_up"])
        s["pt_hp"] = s["n"] / 2 * s["follow_up"]
        s["df"] = 2.0
    return sc


def run(n_bdribs, n_expbern, seed0=20271106):
    rows = []
    for j, s in enumerate(scenarios()):
        rng = np.random.default_rng(seed0 + j)
        n_sim = n_expbern if "EXP_BERN" in s["models"] else n_bdribs
        sims = [simulate(rng, s["n"], 1.0, 26.0, s["follow_up"], s["taus"],
                         s["dist"], s["par_p"], s["par_d"]) for _ in range(n_sim)]
        for look in range(len(s["taus"])):
            data = [sm[look] for sm in sims]
            y = np.array([d["ev"].sum() for d in data])
            for m in s["models"]:
                if m in ("BDRIBS1", "BDRIBS2"):
                    e = np.array([(d["fu"] if m == "BDRIBS1" else d["tar"]).sum()
                                  for d in data])
                    prob = bdribs_prob(y, e, 1.0, s["lambda_hp"], s["pt_hp"],
                                       s["df"])
                    rd = None
                else:
                    res = np.array([expbern_fit(d["ev"], d["fu"], 1.0,
                                                s["lambda_hp"], s["pt_hp"],
                                                s["df"], t_star=s["follow_up"])
                                    for d in data])
                    prob, rd = res[:, 0], res[:, 1]
                p = np.mean(prob > s["cutoff"])
                rows.append([s["figure"], s["point"], look + 1, m, "power",
                             s["cutoff"], p, np.sqrt(p * (1 - p) / n_sim), n_sim])
                if rd is not None:
                    rows.append([s["figure"], s["point"], look + 1, m, "rd", "",
                                 rd.mean(), rd.std(ddof=1) / np.sqrt(n_sim),
                                 n_sim])
                print(rows[-1], flush=True)
    return rows


TEST_T = [38.3, 50.3, 31.8, 45.1, 38.6, 8.5, 19.3, 42.3, 38.7, 40.7, 25.5,
          46.3, 10.7, 5, 22.9, 46.9, 40, 30.7, 40.9, 9.1, 16.2, 50.6, 31.6,
          33.1, 40.6, 20.6, 13.8, 25, 42.5, 30.1, 16.4, 15, 27.2, 47.9, 23.8,
          9.7, 7.9, 6.8, 39.6, 14, 43.4, 29.5, 24.7, 38, 39.2, 25, 36, 25.8,
          27.4, 48.6, 44.2, 6, 6.6, 33.1, 43.9, 32.4, 49.9, 44.5, 26.3, 47.7]
TEST_Y = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
          0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
          0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def check():
    """Compare with the reference values of the tests of the package."""
    out = [
        ("BDRIBS case A", bdribs_prob([5], [500], 2.0, 18 / 4000, 4000, 1.0)[0],
         0.8198879688642029),
        ("BDRIBS case B", bdribs_prob([30], [8000], 1.0, 0.002, 10000, 2.0,
                                      rr0=1.2)[0], 0.8989448489448683),
    ]
    p, d = expbern_fit(TEST_Y, TEST_T, 1.0, 0.004, 2000, 2.0, t_star=52)
    out += [("EXP_BERN setting 1, prob", p, 0.6313681141356096),
            ("EXP_BERN setting 1, difference", d, 0.07354565606388641)]
    p, d = expbern_fit(TEST_Y, TEST_T, 2.0, 0.005, 3000, 1.0, a=1, b=1,
                       rr0=1.5, t_star=26)
    out += [("EXP_BERN setting 2, prob", p, 0.26752577330727506),
            ("EXP_BERN setting 2, difference", d, 0.008163574071360064)]
    for label, v, ref in out:
        print(f"{label:32s} {v:.12f}  test {ref:.12f}  difference {v - ref:+.1e}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        check()
    else:
        rows = run(2000, 1000)
        with open(os.path.join(here, "ref_sim_oc.csv"), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\r\n")
            w.writerow(["figure", "point", "look", "model", "quantity", "cutoff",
                        "value", "se", "n_sim"])
            for r in rows:
                w.writerow([r[0], r[1], r[2], r[3], r[4], r[5], repr(float(r[6])),
                            repr(float(r[7])), r[8]])
