#include <Rcpp.h>
#include <cmath>
#include <vector>

// Log-likelihood of the log-logistic Bernoulli (LL_BERN) mixture model.
// The cumulative incidence at time t is 1 - 1 / (1 + exp((log(t) - mu) / sigma)).
// t: follow-up times; y: event indicators (0 or 1); wp: placebo weight 1 / (k + 1).
// mu_p, sigma_p: placebo parameter pairs (same length).
// mu_d, sigma_d: study drug parameter pairs (same length).
// Returns a matrix with one row per placebo pair and one column per study drug pair.
// [[Rcpp::export]]
Rcpp::NumericMatrix llbern_loglik_cpp(const Rcpp::NumericVector& t,
                                      const Rcpp::IntegerVector& y,
                                      const double wp,
                                      const Rcpp::NumericVector& mu_p,
                                      const Rcpp::NumericVector& sigma_p,
                                      const Rcpp::NumericVector& mu_d,
                                      const Rcpp::NumericVector& sigma_d) {
  const int n = t.size();
  const int np = mu_p.size();
  const int nd = mu_d.size();
  const double wd = 1.0 - wp;
  std::vector<double> lt(n);
  for (int i = 0; i < n; ++i) {
    lt[i] = std::log(t[i]);
  }
  // weighted probabilities of the observed outcome for each parameter pair
  std::vector<double> a(static_cast<size_t>(np) * n);
  std::vector<double> b(static_cast<size_t>(nd) * n);
  for (int p = 0; p < np; ++p) {
    for (int i = 0; i < n; ++i) {
      const double z = (lt[i] - mu_p[p]) / sigma_p[p];
      const double f = (y[i] == 1) ? 1.0 / (1.0 + std::exp(-z)) : 1.0 / (1.0 + std::exp(z));
      a[static_cast<size_t>(p) * n + i] = wp * f;
    }
  }
  for (int d = 0; d < nd; ++d) {
    for (int i = 0; i < n; ++i) {
      const double z = (lt[i] - mu_d[d]) / sigma_d[d];
      const double f = (y[i] == 1) ? 1.0 / (1.0 + std::exp(-z)) : 1.0 / (1.0 + std::exp(z));
      b[static_cast<size_t>(d) * n + i] = wd * f;
    }
  }
  Rcpp::NumericMatrix out(np, nd);
  for (int p = 0; p < np; ++p) {
    const double* ap = &a[static_cast<size_t>(p) * n];
    for (int d = 0; d < nd; ++d) {
      const double* bd = &b[static_cast<size_t>(d) * n];
      double s = 0.0;
      double prod = 1.0;
      for (int i = 0; i < n; ++i) {
        const double x = ap[i] + bd[i];
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
      out(p, d) = s + std::log(prod);
    }
  }
  return out;
}
