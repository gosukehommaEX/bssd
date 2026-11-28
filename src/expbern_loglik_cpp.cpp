#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "fast_exp.h"

// Adds w * log(x) to s. The factors of single subjects (w = 1) are
// multiplied into prod, and the logarithm is taken only when needed.
static inline void add_log_term(const double x, const double w, double& s,
                                double& prod) {
  if (w != 1.0 || x < 1e-20) {
    s += w * std::log(x);
  } else {
    prod *= x;
    if (prod < 1e-280) {
      s += std::log(prod);
      prod = 1.0;
    }
  }
}

// Log-likelihood of the exponential Bernoulli (EXP_BERN) mixture model on a
// product grid of placebo hazard rates (lambda) and hazard ratios (r).
// t: follow-up times; y: event indicators (0 or 1); w: numbers of subjects
// with the same (y, t); wp: placebo weight 1 / (k + 1).
// Returns a matrix with one row per lambda value and one column per r value.
// [[Rcpp::export]]
Rcpp::NumericMatrix expbern_loglik_cpp(const Rcpp::NumericVector& t,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector& w,
                                       const double wp,
                                       const Rcpp::NumericVector& lambda,
                                       const Rcpp::NumericVector& r) {
  const int n = t.size();
  const int nl = lambda.size();
  const int nr = r.size();
  const double wd = 1.0 - wp;
  // subjects without (0) and with (1) the event, so that the loops over the
  // subjects have no branch on y
  std::vector<double> t0, w0, t1, w1;
  for (int i = 0; i < n; ++i) {
    if (y[i] > 0.5) {
      t1.push_back(t[i]);
      w1.push_back(w[i]);
    } else {
      t0.push_back(t[i]);
      w0.push_back(w[i]);
    }
  }
  const int n0 = t0.size();
  const int n1 = t1.size();
  Rcpp::NumericMatrix out(nl, nr);
  std::vector<double> a0(n0), a1(n1);
  for (int l = 0; l < nl; ++l) {
    // placebo parts of the probabilities of the observed outcomes
    for (int i = 0; i < n0; ++i) {
      a0[i] = wp * kernel_exp(-lambda[l] * t0[i]);
    }
    for (int i = 0; i < n1; ++i) {
      a1[i] = wp * -kernel_expm1(-lambda[l] * t1[i]);
    }
    for (int j = 0; j < nr; ++j) {
      const double rl = r[j] * lambda[l];
      double s = 0.0;
      double prod = 1.0;
      for (int i = 0; i < n0; ++i) {
        add_log_term(a0[i] + wd * kernel_exp(-rl * t0[i]), w0[i], s, prod);
      }
      for (int i = 0; i < n1; ++i) {
        add_log_term(a1[i] + wd * -kernel_expm1(-rl * t1[i]), w1[i], s,
                     prod);
      }
      out(l, j) = s + std::log(prod);
    }
  }
  return out;
}
