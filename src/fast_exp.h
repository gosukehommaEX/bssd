#ifndef BSSD_FAST_EXP_H
#define BSSD_FAST_EXP_H

#include <cmath>
#include <cstdint>
#include <cstring>

// exp(x) with an error of about one unit in the last place (at most 1.12 in
// a test of 4 million arguments against long double). x = k log(2) + r with
// |r| <= log(2) / 2: log(2) is split in two parts as in fdlibm, so that k
// times the first part is exact, and k is rounded to the nearest integer by
// adding and subtracting 1.5 * 2^52. exp(r) is its Taylor polynomial of
// degree 13 (truncation error below 5e-18) by Horner's rule, and 2^k is
// built from its exponent bits. Arguments outside (-708, 709), infinities
// and NaN are passed to std::exp.
static inline double fast_exp(const double x) {
  if (!(x > -708.0 && x < 709.0)) {
    return std::exp(x);
  }
  const double inv_ln2 = 1.44269504088896338700e+00;
  const double ln2_hi = 6.93147180369123816490e-01;
  const double ln2_lo = 1.90821492927058770002e-10;
  const double shift = 6755399441055744.0;
  const double kd = (x * inv_ln2 + shift) - shift;
  const double r = (x - kd * ln2_hi) - kd * ln2_lo;
  double p = 1.0 / 6227020800.0;
  p = p * r + 1.0 / 479001600.0;
  p = p * r + 1.0 / 39916800.0;
  p = p * r + 1.0 / 3628800.0;
  p = p * r + 1.0 / 362880.0;
  p = p * r + 1.0 / 40320.0;
  p = p * r + 1.0 / 5040.0;
  p = p * r + 1.0 / 720.0;
  p = p * r + 1.0 / 120.0;
  p = p * r + 1.0 / 24.0;
  p = p * r + 1.0 / 6.0;
  p = p * r + 0.5;
  p = p * r + 1.0;
  p = p * r + 1.0;
  const std::uint64_t bits =
    static_cast<std::uint64_t>(static_cast<std::int64_t>(kd) + 1023) << 52;
  double scale;
  std::memcpy(&scale, &bits, sizeof(double));
  return p * scale;
}

// expm1(x) with an error of about two units in the last place (at most 1.84
// in the same test): x times the Taylor polynomial of (exp(x) - 1) / x of
// degree 16 for |x| < log(2) (truncation error below 4e-19), and
// fast_exp(x) - 1 otherwise, where |exp(x) - 1| >= 0.5.
static inline double fast_expm1(const double x) {
  if (std::fabs(x) < 0.6931471805599453) {
    double p = 1.0 / 355687428096000.0;
    p = p * x + 1.0 / 20922789888000.0;
    p = p * x + 1.0 / 1307674368000.0;
    p = p * x + 1.0 / 87178291200.0;
    p = p * x + 1.0 / 6227020800.0;
    p = p * x + 1.0 / 479001600.0;
    p = p * x + 1.0 / 39916800.0;
    p = p * x + 1.0 / 3628800.0;
    p = p * x + 1.0 / 362880.0;
    p = p * x + 1.0 / 40320.0;
    p = p * x + 1.0 / 5040.0;
    p = p * x + 1.0 / 720.0;
    p = p * x + 1.0 / 120.0;
    p = p * x + 1.0 / 24.0;
    p = p * x + 1.0 / 6.0;
    p = p * x + 0.5;
    p = p * x + 1.0;
    return x * p;
  }
  return fast_exp(x) - 1.0;
}

// exp() and expm1() of the loops over the subjects in the likelihoods: the
// functions above on Windows, where std::exp() and std::expm1() took about
// 50 ns per call (2026-11), and the functions of the C library on other
// systems, where they are faster than the functions above.
static inline double kernel_exp(const double x) {
#ifdef _WIN32
  return fast_exp(x);
#else
  return std::exp(x);
#endif
}

static inline double kernel_expm1(const double x) {
#ifdef _WIN32
  return fast_expm1(x);
#else
  return std::expm1(x);
#endif
}

#endif
