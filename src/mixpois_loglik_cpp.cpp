#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

// Log of the Poisson kernel lambda^y * exp(-lambda * t) (the factor
// t^y / y! is the same for both arms and is omitted).
static inline double pois_kernel(const double y, const double lambda,
                                 const double t) {
  return (y > 0.0 ? y * std::log(lambda) : 0.0) - lambda * t;
}

// Log-likelihood of the blinded exposure-time Poisson model of Schnell and
// Ball (2016) on a product grid of control rates (lc) and study drug rates
// (ld). t: follow-up times; y: numbers of events; w: numbers of subjects with
// the same (y, t); wc: probability of the control arm 1 / (k + 1).
// Returns a matrix with one row per lc value and one column per ld value.
// [[Rcpp::export]]
Rcpp::NumericMatrix mixpois_loglik_cpp(const Rcpp::NumericVector& t,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector& w,
                                       const double wc,
                                       const Rcpp::NumericVector& lc,
                                       const Rcpp::NumericVector& ld) {
  const int n = t.size();
  const int nc = lc.size();
  const int nd = ld.size();
  const double ninf = -std::numeric_limits<double>::infinity();
  const double lwc = std::log(wc);
  const double lwd = std::log1p(-wc);
  Rcpp::NumericMatrix out(nc, nd);
  std::vector<double> la(n);
  std::vector<double> lb(static_cast<size_t>(n) * nd);
  for (int j = 0; j < nd; ++j) {
    for (int i = 0; i < n; ++i) {
      lb[static_cast<size_t>(j) * n + i] = lwd + pois_kernel(y[i], ld[j], t[i]);
    }
  }
  for (int l = 0; l < nc; ++l) {
    for (int i = 0; i < n; ++i) {
      la[i] = lwc + pois_kernel(y[i], lc[l], t[i]);
    }
    for (int j = 0; j < nd; ++j) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) {
        const double a = la[i];
        const double b = lb[static_cast<size_t>(j) * n + i];
        const double m = (a > b) ? a : b;
        if (m == ninf) {
          s = ninf;
          break;
        }
        s += w[i] * (m + std::log1p(std::exp(-std::fabs(a - b))));
      }
      out(l, j) = s;
    }
  }
  return out;
}
