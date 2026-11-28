#include <Rcpp.h>
#include <cmath>
#include <vector>

// Factor x of a single subject: multiplied into q, or its logarithm added
// to s if x is below 1e-20.
static inline void single_term(const double x, double& s, double& q) {
  if (x < 1e-20) {
    s += std::log(x);
  } else {
    q *= x;
  }
}

// Moves a partial product below 'bound' to the sum of logarithms s.
static inline void move_small(double& q, const double bound, double& s) {
  if (q < bound) {
    s += std::log(q);
    q = 1.0;
  }
}

// Log-likelihood of the log-logistic Bernoulli (LL_BERN) mixture model.
// The cumulative incidence at time t is 1 - 1 / (1 + exp((log(t) - mu) / sigma)).
// t: follow-up times; y: event indicators (0 or 1); w: numbers of subjects
// with the same (y, t); wp: placebo weight 1 / (k + 1).
// mu_p, sigma_p: placebo parameter pairs (same length).
// mu_d, sigma_d: study drug parameter pairs (same length).
// Returns a matrix with one row per placebo pair and one column per study drug pair.
// [[Rcpp::export]]
Rcpp::NumericMatrix llbern_loglik_cpp(const Rcpp::NumericVector& t,
                                      const Rcpp::NumericVector& y,
                                      const Rcpp::NumericVector& w,
                                      const double wp,
                                      const Rcpp::NumericVector& mu_p,
                                      const Rcpp::NumericVector& sigma_p,
                                      const Rcpp::NumericVector& mu_d,
                                      const Rcpp::NumericVector& sigma_d) {
  const int n = t.size();
  const int np = mu_p.size();
  const int nd = mu_d.size();
  const double wd = 1.0 - wp;
  // single subjects (w = 1) first, then groups of several subjects
  std::vector<int> idx;
  for (int i = 0; i < n; ++i) {
    if (w[i] == 1.0) {
      idx.push_back(i);
    }
  }
  const int ns = idx.size();
  for (int i = 0; i < n; ++i) {
    if (w[i] != 1.0) {
      idx.push_back(i);
    }
  }
  std::vector<double> lt(n), wg(n);
  std::vector<int> ev(n);
  for (int k = 0; k < n; ++k) {
    lt[k] = std::log(t[idx[k]]);
    ev[k] = y[idx[k]] > 0.5;
    wg[k] = w[idx[k]];
  }
  // weighted probabilities of the observed outcome for each parameter pair
  std::vector<double> a(static_cast<size_t>(np) * n);
  std::vector<double> b(static_cast<size_t>(nd) * n);
  for (int p = 0; p < np; ++p) {
    for (int k = 0; k < n; ++k) {
      const double z = (lt[k] - mu_p[p]) / sigma_p[p];
      const double f = ev[k] ? 1.0 / (1.0 + std::exp(-z)) : 1.0 / (1.0 + std::exp(z));
      a[static_cast<size_t>(p) * n + k] = wp * f;
    }
  }
  for (int d = 0; d < nd; ++d) {
    for (int k = 0; k < n; ++k) {
      const double z = (lt[k] - mu_d[d]) / sigma_d[d];
      const double f = ev[k] ? 1.0 / (1.0 + std::exp(-z)) : 1.0 / (1.0 + std::exp(z));
      b[static_cast<size_t>(d) * n + k] = wd * f;
    }
  }
  const int ns4 = ns - ns % 4;
  Rcpp::NumericMatrix out(np, nd);
  for (int p = 0; p < np; ++p) {
    const double* ap = &a[static_cast<size_t>(p) * n];
    for (int d = 0; d < nd; ++d) {
      const double* bd = &b[static_cast<size_t>(d) * n];
      double s = 0.0;
      // the factors of single subjects are multiplied into four partial
      // products, so that each multiplication need not wait for the
      // previous one; every factor is at least 1e-20, and a product below
      // 1e-70 is moved to s after every four factors, so that the products
      // stay above 1e-150
      double q0 = 1.0;
      double q1 = 1.0;
      double q2 = 1.0;
      double q3 = 1.0;
      for (int i = 0; i < ns4; i += 4) {
        const double x0 = ap[i] + bd[i];
        const double x1 = ap[i + 1] + bd[i + 1];
        const double x2 = ap[i + 2] + bd[i + 2];
        const double x3 = ap[i + 3] + bd[i + 3];
        if (x0 < 1e-20 || x1 < 1e-20 || x2 < 1e-20 || x3 < 1e-20) {
          single_term(x0, s, q0);
          single_term(x1, s, q1);
          single_term(x2, s, q2);
          single_term(x3, s, q3);
        } else {
          q0 *= x0;
          q1 *= x1;
          q2 *= x2;
          q3 *= x3;
        }
        if ((i & 15) == 12) {
          move_small(q0, 1e-70, s);
          move_small(q1, 1e-70, s);
          move_small(q2, 1e-70, s);
          move_small(q3, 1e-70, s);
        }
      }
      for (int i = ns4; i < ns; ++i) {
        single_term(ap[i] + bd[i], s, q0);
      }
      // groups of several subjects
      for (int i = ns; i < n; ++i) {
        s += wg[i] * std::log(ap[i] + bd[i]);
      }
      // the products are now above 1e-190; after these moves they are above
      // 1e-140, so that the products of two stay above 1e-280
      move_small(q0, 1e-140, s);
      move_small(q1, 1e-140, s);
      move_small(q2, 1e-140, s);
      move_small(q3, 1e-140, s);
      const double q01 = q0 * q1;
      const double q23 = q2 * q3;
      if (q01 > 1e-150 && q23 > 1e-150) {
        s += std::log(q01 * q23);
      } else {
        s += std::log(q01) + std::log(q23);
      }
      out(p, d) = s;
    }
  }
  return out;
}
