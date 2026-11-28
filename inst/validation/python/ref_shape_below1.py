"""Reference values for binbeta_post() and exptime_post() with a control
prior whose shape parameter is below 1 (cases BS1 and ES1), for which the
functions integrate over the normal quantile of the control prior.

BS1 (binbeta_post): the binomial expansion of the likelihood
(wc theta_c + (1 - wc) theta_d)^y (wc (1 - theta_c) + (1 - wc) (1 - theta_d))^(n - y)
gives the posterior as a finite mixture of products of Beta distributions
with positive weights; P(theta_d > theta_c) of each product is a
one-dimensional integral computed with scipy.integrate.quad.

ES1 (exptime_post, data set E1 of ref_common.py): composite Gauss-Legendre
rules on (z_c, z_d), z = qnorm(F(lambda)) with F the prior distribution
function, on [-9, 9] x [-9, 9]; for each node of z_c the rule of z_d is
split exactly at the threshold (rate0 or rr0 * lambda_c). The computation
is repeated with more nodes, and the largest change is reported.
"""
import json
import numpy as np
from scipy import special, stats, integrate
import ref_common as rc


def binbeta_ref(y, n, k, ac, bc, ad, bd, p0):
    wc = 1 / (k + 1)
    i = np.arange(y + 1)[:, None]          # control part of the events
    j = np.arange(n - y + 1)[None, :]      # control part of the non-events
    lw = (special.gammaln(y + 1) - special.gammaln(i + 1) - special.gammaln(y - i + 1)
          + special.gammaln(n - y + 1) - special.gammaln(j + 1) - special.gammaln(n - y - j + 1)
          + (i + j) * np.log(wc) + (n - i - j) * np.log1p(-wc)
          + special.betaln(ac + i, bc + j) + special.betaln(ad + y - i, bd + n - y - j)
          - special.betaln(ac, bc) - special.betaln(ad, bd))
    w = np.exp(lw - lw.max())
    w = (w / w.sum()).ravel()
    A1 = (ac + i + 0 * j).ravel(); B1 = (bc + j + 0 * i).ravel()
    A2 = (ad + y - i + 0 * j).ravel(); B2 = (bd + n - y - j + 0 * i).ravel()
    keep = w > 1e-18
    w, A1, B1, A2, B2 = w[keep] / w[keep].sum(), A1[keep], B1[keep], A2[keep], B2[keep]
    p_rate = float(np.sum(w * stats.beta.sf(p0, A2, B2)))
    p_ratio = 0.0
    for wk, a1, b1, a2, b2 in zip(w, A1, B1, A2, B2):
        # P(theta_d > theta_c) = int f_c(x) P(theta_d > x) dx
        f = lambda x: stats.beta.pdf(x, a1, b1) * stats.beta.sf(x, a2, b2)
        val = integrate.quad(f, 0, 1, epsabs=1e-14, epsrel=1e-13, limit=400, points=[a1 / (a1 + b1)])[0]
        p_ratio += wk * val
    return dict(prob_rate=p_rate, prob_ratio=p_ratio, mean_theta_c=float(np.sum(w * A1 / (A1 + B1))),
                mean_theta_d=float(np.sum(w * A2 / (A2 + B2))), n_terms=int(len(w)))


def comp(a, b, n_pan, q):
    x, wq = np.polynomial.legendre.leggauss(q)
    e = np.linspace(a, b, n_pan + 1)
    X = ((e[1:] - e[:-1])[:, None] * x / 2 + (e[1:] + e[:-1])[:, None] / 2).ravel()
    W = ((e[1:] - e[:-1])[:, None] * wq / 2).ravel()
    return X, W


def exptime_ref(y, t, k, prior_c, prior_d, rate0, rr0, n_pan, q):
    wc = 1 / (k + 1)
    dc = stats.gamma(prior_c[0], scale=1 / prior_c[1])
    dd = stats.gamma(prior_d[0], scale=1 / prior_d[1])
    y = np.asarray(y, float); t = np.asarray(t, float)

    def lam_of(dist, z):
        return np.where(z <= 0, dist.ppf(special.ndtr(z)), dist.isf(special.ndtr(-z)))

    def z_of(dist, x):
        lf = dist.logcdf(x); ls = dist.logsf(x)
        return np.where(lf <= np.log(0.5), special.ndtri_exp(lf), -special.ndtri_exp(ls))

    def loglik(lc, LD):
        # sum over subjects of log(wc lc^y e^(-lc t) + (1 - wc) ld^y e^(-ld t))
        with np.errstate(divide="ignore"):
            la = np.log(wc) + np.where(y > 0, y * np.log(lc), 0.0) - lc * t
            lb = np.log1p(-wc) + np.where(y[:, None] > 0, y[:, None] * np.log(LD[None, :]), 0.0) - LD[None, :] * t[:, None]
        return np.logaddexp(la[:, None], lb).sum(0)

    zc, wzc = comp(-9, 9, n_pan, q)
    lc = lam_of(dc, zc)
    rows = []
    for i in range(len(zc)):
        out = {}
        for key, thr in (("rate", rate0), ("ratio", rr0 * lc[i])):
            zt = float(np.clip(z_of(dd, np.array([thr]))[0], -9, 9))
            x1, w1 = comp(-9, zt, n_pan // 2, q) if zt > -9 else (np.zeros(0), np.zeros(0))
            x2, w2 = comp(zt, 9, n_pan // 2, q) if zt < 9 else (np.zeros(0), np.zeros(0))
            X = np.concatenate([x1, x2]); W = np.concatenate([w1, w2])
            LD = lam_of(dd, X)
            L = loglik(lc[i], LD) + stats.norm.logpdf(X) + np.log(W)
            out[key] = (L, len(x1), LD)
        rows.append(out)
    lw_c = stats.norm.logpdf(zc) + np.log(wzc)
    m = max(r["ratio"][0].max() + lw_c[i] for i, r in enumerate(rows))
    tot = p_rate = p_ratio = m_c = m_d = 0.0
    for i, r in enumerate(rows):
        L, n1, LD = r["ratio"]
        e = np.exp(L + lw_c[i] - m)
        tot += e.sum(); p_ratio += e[n1:].sum(); m_c += lc[i] * e.sum(); m_d += (e * LD).sum()
        L, n1, _ = r["rate"]
        p_rate += np.exp(L + lw_c[i] - m)[n1:].sum()
    return dict(prob_rate=p_rate / tot, prob_ratio=p_ratio / tot, mean_lambda_c=m_c / tot, mean_lambda_d=m_d / tot)


if __name__ == "__main__":
    res = []
    b = binbeta_ref(12, 100, 1, 0.5, 5.0, 1.0, 1.0, 0.15)
    b.update(id="BS1", y=12, n=100, k=1, prior_c=[0.5, 5.0], prior_d=[1.0, 1.0], p0=0.15, rr0=1.0)
    print(json.dumps(b), flush=True)
    res.append(b)
    d = rc.make_data()["E1"]
    e1 = exptime_ref(d["y"], d["t"], 2, (0.7, 30.0), (2.0, 70.0), 0.04, 1.0, 18, 16)
    e2 = exptime_ref(d["y"], d["t"], 2, (0.7, 30.0), (2.0, 70.0), 0.04, 1.0, 28, 20)
    e2.update(id="ES1", data="E1", k=2, prior_c=[0.7, 30.0], prior_d=[2.0, 70.0], rate0=0.04, rr0=1.0,
              change=max(abs(e2[key] - e1[key]) for key in e1))
    print(json.dumps(e2), flush=True)
    res.append(e2)
    rc.dump({"note": "binbeta_post() and exptime_post() with a control prior shape below 1", "results": res},
            "ref_shape_below1.json")
