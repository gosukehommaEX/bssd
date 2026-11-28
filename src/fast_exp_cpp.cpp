#include <Rcpp.h>
#include "fast_exp.h"

// fast_exp(x) (m1 = false) or fast_expm1(x) (m1 = true) of fast_exp.h for a
// vector, on every system, for the tests of their accuracy.
// [[Rcpp::export]]
Rcpp::NumericVector fast_exp_cpp(const Rcpp::NumericVector& x,
                                 const bool m1) {
  const int n = x.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = m1 ? fast_expm1(x[i]) : fast_exp(x[i]);
  }
  return out;
}
