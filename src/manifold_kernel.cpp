// The projection of positive states onto the conservation laws C x = T of a steady state (Pimpl):
// x exp(C' lambda), which minimises the relative entropy to x, with lambda by Newton and a halving
// line search on the squared residual. It ends where the residual is within the rounding of its
// sums, or where no step lowers it; a singular Newton matrix ends it as in solve() of R.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "pimpl_kernel.h"

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

namespace dmod {
namespace pimpl {

void onManifold(double* xa, int n, const double* C, int k, const double* total) {
  if (k == 0 || n == 0) return;
  const double eps = std::numeric_limits<double>::epsilon();
  const std::vector<double> x0(xa, xa + n);
  auto Cij = [&](int i, int j) { return C[i + (std::size_t) k * j]; };
  std::vector<double> lam(k, 0.0), e(n), r(k), rn(k), step(k), trial(k), H((std::size_t) k * k), work(4 * k);
  std::vector<int> ipiv(k), iwork(k);
  auto states = [&](const std::vector<double>& l, std::vector<double>& out) {
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int i = 0; i < k; ++i) s += Cij(i, j) * l[i];
      out[j] = x0[j] * std::exp(s);
    }
  };
  auto residual = [&](const std::vector<double>& l, std::vector<double>& out) {
    std::vector<double> x(n);
    states(l, x);
    for (int i = 0; i < k; ++i) {
      double s = 0.0;
      for (int j = 0; j < n; ++j) s += Cij(i, j) * x[j];
      out[i] = s - total[i];
    }
  };
  auto sq = [](const std::vector<double>& v) { double s = 0.0; for (double a : v) s += a * a; return s; };
  residual(lam, r);
  for (int it = 0; it < 100; ++it) {
    states(lam, e);
    bool within = true;
    for (int i = 0; i < k && within; ++i) {
      double s = std::fabs(total[i]);
      for (int j = 0; j < n; ++j) s += std::fabs(Cij(i, j)) * e[j];
      within = std::fabs(r[i]) <= (n + 1) * eps * s;
    }
    if (within) break;
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += Cij(a, j) * e[j] * Cij(b, j);
        H[a + (std::size_t) k * b] = s;
      }
    double anorm = 0.0;
    for (int b = 0; b < k; ++b) {
      double s = 0.0;
      for (int a = 0; a < k; ++a) s += std::fabs(H[a + (std::size_t) k * b]);
      anorm = std::max(anorm, s);
    }
    int info = 0, one = 1;
    F77_CALL(dgetrf)(&k, &k, H.data(), &k, ipiv.data(), &info);
    if (info != 0) break;
    double rcond = 0.0;
    F77_CALL(dgecon)("1", &k, H.data(), &k, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
    if (rcond < eps) break;
    step = r;
    F77_CALL(dgetrs)("N", &k, &one, H.data(), &k, ipiv.data(), step.data(), &k, &info FCONE);
    double a = 1.0;
    bool lowered = false;
    while (a >= 1e-8) {
      for (int i = 0; i < k; ++i) trial[i] = lam[i] - a * step[i];
      residual(trial, rn);
      bool finite = true;
      for (double v : rn) finite = finite && std::isfinite(v);
      if (finite && sq(rn) < sq(r)) { lowered = true; break; }
      a /= 2;
    }
    if (!lowered) break;
    lam = trial;
    r = rn;
  }
  states(lam, e);
  std::copy(e.begin(), e.end(), xa);
}

}  // namespace pimpl
}  // namespace dmod

// [[Rcpp::export]]
NumericVector pimpl_on_manifold(NumericVector xa, NumericMatrix C, NumericVector total) {
  NumericVector out = clone(xa);
  dmod::pimpl::onManifold(REAL(out), out.size(), REAL(C), C.nrow(), REAL(total));
  return out;
}
