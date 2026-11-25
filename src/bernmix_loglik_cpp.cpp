#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

// Log of the Bernoulli probability of y (0 or 1) when the event
// probability over time t is 1 - exp(-h * t).
static inline double bern_kernel(const double y, const double h,
                                 const double t) {
  return (y > 0.5) ? std::log(-std::expm1(-h * t)) : -h * t;
}

// Log-likelihood of a blinded Bernoulli mixture with constant hazards on a
// product grid of control hazards (hc) and study drug hazards (hd). With
// hazards h = -log(1 - pi), the event probability over time t is
// 1 - (1 - pi)^t, which is the model of Brock et al. (2023).
// t: times at risk; y: event indicators; w: numbers of subjects with the
// same (y, t); wc: probability of the control arm 1 / (k + 1).
// Returns a matrix with one row per hc value and one column per hd value.
// [[Rcpp::export]]
Rcpp::NumericMatrix bernmix_loglik_cpp(const Rcpp::NumericVector& t,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector& w,
                                       const double wc,
                                       const Rcpp::NumericVector& hc,
                                       const Rcpp::NumericVector& hd) {
  const int n = t.size();
  const int nc = hc.size();
  const int nd = hd.size();
  const double ninf = -std::numeric_limits<double>::infinity();
  const double lwc = std::log(wc);
  const double lwd = std::log1p(-wc);
  Rcpp::NumericMatrix out(nc, nd);
  std::vector<double> la(n);
  std::vector<double> lb(static_cast<size_t>(n) * nd);
  for (int j = 0; j < nd; ++j) {
    for (int i = 0; i < n; ++i) {
      lb[static_cast<size_t>(j) * n + i] = lwd + bern_kernel(y[i], hd[j], t[i]);
    }
  }
  for (int l = 0; l < nc; ++l) {
    for (int i = 0; i < n; ++i) {
      la[i] = lwc + bern_kernel(y[i], hc[l], t[i]);
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
