#include <Rcpp.h>
#include <cmath>
#include <vector>

// Log-likelihood of the exponential Bernoulli (EXP_BERN) mixture model on a
// product grid of placebo hazard rates (lambda) and hazard ratios (r).
// t: follow-up times; y: event indicators (0 or 1); wp: placebo weight 1 / (k + 1).
// Returns a matrix with one row per lambda value and one column per r value.
// [[Rcpp::export]]
Rcpp::NumericMatrix expbern_loglik_cpp(const Rcpp::NumericVector& t,
                                       const Rcpp::IntegerVector& y,
                                       const double wp,
                                       const Rcpp::NumericVector& lambda,
                                       const Rcpp::NumericVector& r) {
  const int n = t.size();
  const int nl = lambda.size();
  const int nr = r.size();
  const double wd = 1.0 - wp;
  Rcpp::NumericMatrix out(nl, nr);
  std::vector<double> a(n);
  for (int l = 0; l < nl; ++l) {
    for (int i = 0; i < n; ++i) {
      const double h = lambda[l] * t[i];
      a[i] = (y[i] == 1) ? -std::expm1(-h) : std::exp(-h);
    }
    for (int j = 0; j < nr; ++j) {
      const double rl = r[j] * lambda[l];
      double s = 0.0;
      double prod = 1.0;
      for (int i = 0; i < n; ++i) {
        const double h = rl * t[i];
        const double b = (y[i] == 1) ? -std::expm1(-h) : std::exp(-h);
        const double x = wp * a[i] + wd * b;
        // accumulate the product and take logarithms only when needed
        if (x < 1e-20) {
          s += std::log(x);
        } else {
          prod *= x;
          if (prod < 1e-280) {
            s += std::log(prod);
            prod = 1.0;
          }
        }
      }
      out(l, j) = s + std::log(prod);
    }
  }
  return out;
}
