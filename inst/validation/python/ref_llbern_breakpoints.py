"""Independent reference values for llbern_post when the posterior of sigma_d
has much mass near zero.

When sigma_d is small, the likelihood changes sharply in mu_d near the
logarithms of the follow-up times. This script integrates the four-dimensional
posterior without the box refinement and the Chebyshev interpolation of the
package:

- log(mu_p) and log(sigma_p): Gauss-Legendre (GL) nodes on the range between
  the prior quantiles q_tail and 1 - q_tail;
- sigma_d: its prior distribution-function scale v, with composite GL panels
  that are graded towards v = 0;
- mu_d itself: composite GL panels with breakpoints at log(t_i) + c * sigma_d
  (c = -8, -3, 0, 3, 8) for every distinct follow-up time, panels at most
  4 * sigma_d wide near these points and at most w_smooth wide elsewhere; the
  panel that contains the threshold mu_star (RR(t_star) = rr0) is split at
  mu_star, so that P(RR(t_star) > rr0 | data) has no discontinuous integrand.

Cases (output ref_llbern_breakpoints.csv):

- L1: the data set L1 of validate_convergence.R (sim_blinded with seed 201;
  data_llbern_L1.csv was written by R), computed with several settings to
  show that the reference does not depend on them;
- test_*: the data of test-llbern_post.R, to compare this script with the
  reference values of ref_llbern.py and ref_ll_q.py.

Run from this folder: python3 ref_llbern_breakpoints.py (about 15 minutes).
"""
import csv
import os
import numpy as np
from scipy import special, stats


def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * x + (b + a) / 2, w * (b - a) / 2


def load(path):
    """Read y and t and collapse identical (t, y) pairs into counts."""
    rows = list(csv.DictReader(open(path, newline="")))
    t = np.array([float(r["t"]) for r in rows])
    y = np.array([int(r["y"]) for r in rows])
    return collapse(t, y)


def collapse(t, y):
    uniq = {}
    for ti, yi in zip(np.round(t, 12), y):
        uniq[(ti, yi)] = uniq.get((ti, yi), 0) + 1
    tt = np.array([key[0] for key in uniq])
    yy = np.array([key[1] for key in uniq])
    nn = np.array(list(uniq.values()), dtype=float)
    return tt, yy, nn


def reference(tt, yy, nn, k, mu_hp, sg_hp, pt_hp, df_p, t_star, rr0=1.0,
              n_p=12, n_md=8, n_sd=8, q_tail=1e-12, md_tail=1e-14,
              win=(-8, -3, 0, 3, 8), w_smooth=1.0):
    """P(RR(t_star) > rr0 | data) and the posterior mean of F_d(t_star)."""
    wp, wd = 1 / (k + 1), k / (k + 1)
    g_mp = stats.gamma(pt_hp * mu_hp / df_p, scale=df_p / pt_hp)
    g_sp = stats.gamma(pt_hp * sg_hp / df_p, scale=df_p / pt_hp)
    g_md = stats.gamma(mu_hp, scale=1.0)
    g_sd = stats.gamma(sg_hp, scale=1.0)
    lt = np.log(tt)
    lts = np.log(t_star)
    # placebo parameters on the log scale: prior density times Jacobian
    ump, wmp = gl(np.log(g_mp.ppf(q_tail)), np.log(g_mp.ppf(1 - q_tail)), n_p)
    usp, wsp = gl(np.log(g_sp.ppf(q_tail)), np.log(g_sp.ppf(1 - q_tail)), n_p)
    MP, SP = [a.ravel() for a in np.meshgrid(np.exp(ump), np.exp(usp),
                                             indexing="ij")]
    WP = np.outer(wmp * g_mp.pdf(np.exp(ump)) * np.exp(ump),
                  wsp * g_sp.pdf(np.exp(usp)) * np.exp(usp)).ravel()
    zp = (lt[None, :] - MP[:, None]) / SP[:, None]
    A = np.where(yy[None, :] == 1, wp * special.expit(zp),
                 wp * special.expit(-zp))
    fp_star = special.expit((lts - MP) / SP)
    # sigma_d on its prior distribution-function scale, graded panels
    v_edges = [0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2,
               0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
               0.99, 0.999, 0.99999, 1 - 1e-8, 1 - 1e-12]
    vs, wv = [], []
    for a, b in zip(v_edges[:-1], v_edges[1:]):
        x, w = gl(a, b, n_sd)
        vs.append(x)
        wv.append(w)
    vs = np.concatenate(vs)
    wv = np.concatenate(wv)
    sds = np.maximum(g_sd.ppf(vs), 1e-300)
    md_lo, md_hi = g_md.ppf(md_tail), g_md.ppf(1 - md_tail)
    centers = np.unique(lt)
    n_tr = len(MP)
    acc = []
    for sd, w_s in zip(sds, wv):
        # breakpoints in mu_d
        cand = [md_lo, md_hi]
        cand += list(np.arange(np.floor(md_lo), md_hi, w_smooth))
        for c in win:
            cand += list(centers + c * sd)
        cand = np.unique(np.clip(np.array(cand), md_lo, md_hi))
        keep = [cand[0]]
        for x in cand[1:]:
            if x - keep[-1] >= max(sd / 2, 1e-14):
                keep.append(x)
        if keep[-1] < md_hi:
            keep[-1] = md_hi if md_hi - keep[-1] < max(sd / 2, 1e-14) else keep[-1]
            if keep[-1] < md_hi:
                keep.append(md_hi)
        bp = [keep[0]]
        for x in keep[1:]:
            a = bp[-1]
            # a panel inside a window [c - 8 sd, c + 8 sd] is at most 4 sd wide
            inside = np.min(np.abs((a + x) / 2 - centers)) <= max(abs(win[0]), win[-1]) * sd
            wmax = 4 * sd if inside else w_smooth
            m = max(int(np.ceil((x - a) / wmax)), 1)
            bp += list(a + (x - a) * np.arange(1, m + 1) / m)
        bp = np.array(bp)
        # GL nodes on all panels
        xg, wg = np.polynomial.legendre.leggauss(n_md)
        half = np.diff(bp) / 2
        midp = (bp[1:] + bp[:-1]) / 2
        M = (midp[:, None] + half[:, None] * xg[None, :]).ravel()
        Wm = (half[:, None] * wg[None, :]).ravel() * g_md.pdf(M)
        zd = (lt[None, :] - M[:, None]) / sd
        B = np.where(yy[None, :] == 1, wd * special.expit(zd),
                     wd * special.expit(-zd))
        # likelihood for each (placebo pair, mu_d node), product over subjects
        L = np.ones((n_tr, len(M)))
        for i in range(len(tt)):
            L *= (A[:, i][:, None] + B[:, i][None, :]) ** nn[i]
        G = L * Wm[None, :]
        panel = G.reshape(n_tr, len(bp) - 1, n_md).sum(2)
        cum = np.concatenate([np.zeros((n_tr, 1)), np.cumsum(panel, 1)], 1)
        full = cum[:, -1]
        fd_star = special.expit((lts - M) / sd)
        incd = (G * fd_star[None, :]).sum(1)
        # threshold: RR(t_star) > rr0 if and only if mu_d < mu_star
        arg = rr0 * fp_star
        ok = arg < 1
        mstar = np.where(ok, lts - sd * special.logit(np.minimum(arg, 1 - 1e-16)),
                         -np.inf)
        pr = np.zeros(n_tr)
        for p in range(n_tr):
            ms = mstar[p]
            if ms <= bp[0]:
                continue
            if ms >= bp[-1]:
                pr[p] = full[p]
                continue
            j = np.searchsorted(bp, ms) - 1
            xa, wa = gl(bp[j], ms, n_md)
            za = (lt[None, :] - xa[:, None]) / sd
            Ba = np.where(yy[None, :] == 1, wd * special.expit(za),
                          wd * special.expit(-za))
            La = np.prod((A[p][None, :] + Ba) ** nn[None, :], axis=1)
            pr[p] = cum[p, j] + np.sum(La * wa * g_md.pdf(xa))
        acc.append((w_s, (WP * full).sum(), (WP * pr).sum(), (WP * incd).sum()))
    acc = np.array(acc)
    tot = (acc[:, 0] * acc[:, 1]).sum()
    prob = (acc[:, 0] * acc[:, 2]).sum() / tot
    ipd = (acc[:, 0] * acc[:, 3]).sum() / tot
    return prob, ipd


# data of make_data() in tests/testthat/test-llbern_post.R
TEST_T = [66.5, 93.7, 81.6, 26.5, 34, 91.4, 4.5, 86.1, 83.7, 50.8, 34.3,
          31.8, 29.5, 48.5, 54.5, 59.3, 103.6, 83.3, 66.2, 102.9, 25.5, 20,
          65.3, 8.4, 7.6, 55.5, 50.6, 95.7, 66.9, 55.4, 53.7, 28.8, 5.2, 23.2,
          73.2, 24.1, 41, 4.4, 87, 19.4, 30.8, 92, 55, 88.7, 68, 78.2, 13.1,
          58.1, 54.8, 91.1, 40.1, 63.8, 9.9, 42.8, 36.3, 19, 85.6, 41.9, 101.9,
          63]
TEST_Y = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    # L1: n = 150, 104 weeks of follow-up, prior from mu = 7.3, sigma = 0.8,
    # PT_HP = 150 / 2 * 104 / 52, DF_p = 2, t_star = 104
    tt, yy, nn = load(os.path.join(here, "data_llbern_L1.csv"))
    l1 = dict(k=1.0, mu_hp=7.3, sg_hp=0.8, pt_hp=150.0, df_p=2.0, t_star=104.0)
    settings = [
        ("reference", dict(n_p=32)),
        ("n_p = 24", dict(n_p=24)),
        ("n_p = 40", dict(n_p=40)),
        ("n_p = 24, n_md = 12", dict(n_p=24, n_md=12)),
        ("n_p = 24, n_sd = 12", dict(n_p=24, n_sd=12)),
        ("n_p = 24, windows -10 to 10", dict(n_p=24, win=(-10, -5, -2, 0, 2, 5, 10))),
        ("n_p = 24, q_tail = 1e-14, md_tail = 1e-16",
         dict(n_p=24, q_tail=1e-14, md_tail=1e-16)),
        ("n_p = 24, w_smooth = 0.5", dict(n_p=24, w_smooth=0.5)),
    ]
    for label, s in settings:
        p, ipd = reference(tt, yy, nn, rr0=1.0, **l1, **s)
        rows.append(["L1", label, p, ipd, ""])
        print("L1", label, p, ipd, flush=True)
    # test data: reference values of ref_llbern.py and ref_ll_q.py in
    # test-llbern_post.R
    tt, yy, nn = collapse(np.array(TEST_T, dtype=float), np.array(TEST_Y))
    tests = [
        ("test_k1_rr1", dict(k=1.0, df_p=1.0, t_star=104.0, rr0=1.0), 0.6964229),
        ("test_k1_rr2", dict(k=1.0, df_p=1.0, t_star=104.0, rr0=2.0), 0.33036),
        ("test_k2_rr1.5", dict(k=2.0, df_p=2.0, t_star=52.0, rr0=1.5),
         0.6092996462522098),
    ]
    for case, s, other in tests:
        p, ipd = reference(tt, yy, nn, mu_hp=6.0, sg_hp=0.8, pt_hp=200.0,
                           n_p=24, **s)
        rows.append([case, "n_p = 24", p, ipd, other])
        print(case, p, ipd, other, flush=True)
    with open(os.path.join(here, "ref_llbern_breakpoints.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh, lineterminator="\r\n")
        w.writerow(["case", "setting", "prob", "ip_drug", "other_reference"])
        for r in rows:
            w.writerow([r[0], r[1], repr(float(r[2])), repr(float(r[3])),
                        "" if r[4] == "" else repr(float(r[4]))])
