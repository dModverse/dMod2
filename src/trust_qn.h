// Quasi-Newton update of a dense model Hessian, shared by the trust-region
// kernel and the partitioned updates of the multiple-shooting driver. Free of
// Rcpp, like trust_subproblem.h.

#ifndef DMOD_TRUST_QN_H
#define DMOD_TRUST_QN_H

#include <vector>
#include <cmath>
#include <cstddef>

namespace dmod { namespace trust_internal {

// Values of `kind`, the same as the kernel's HessianMethod.
enum QnKind { QN_BFGS = 1, QN_SR1 = 2 };

// Dense quasi-Newton update of `B` in place, all arguments in the minimised
// sign convention and `s` in the frame `B` lives in. BFGS is Powell-damped and
// stays positive definite; SR1 may be indefinite. Returns whether `B` changed.
inline bool qn_update_B(int kind, int K, std::vector<double>& B,
                        const std::vector<double>& s, const std::vector<double>& y) {
  std::vector<double> Bs(K, 0.0);
  for (int j = 0; j < K; ++j) {
    const double sj = s[j];
    if (sj != 0.0)
      for (int i = 0; i < K; ++i) Bs[i] += B[i + (std::size_t) j * K] * sj;
  }
  double sy = 0.0, sBs = 0.0;
  for (int i = 0; i < K; ++i) { sy += s[i] * y[i]; sBs += s[i] * Bs[i]; }

  if (kind == QN_BFGS) {
    if (!(sBs > 0.0)) return false;           // no curvature reference; keep B
    double theta = 1.0;
    if (sy < 0.2 * sBs) theta = (0.8 * sBs) / (sBs - sy);
    std::vector<double> rv(K);
    for (int i = 0; i < K; ++i) rv[i] = theta * y[i] + (1.0 - theta) * Bs[i];
    double sr = 0.0;
    for (int i = 0; i < K; ++i) sr += s[i] * rv[i];
    if (!(sr > 0.0)) return false;
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        B[i + (std::size_t) j * K] += rv[i] * rv[j] / sr - Bs[i] * Bs[j] / sBs;
    return true;
  }
  // SR1
  std::vector<double> w(K);
  for (int i = 0; i < K; ++i) w[i] = y[i] - Bs[i];
  double ws = 0.0, wn = 0.0, sn = 0.0;
  for (int i = 0; i < K; ++i) { ws += w[i] * s[i]; wn += w[i] * w[i]; sn += s[i] * s[i]; }
  if (std::fabs(ws) <= 1e-8 * std::sqrt(sn * wn) || !(wn > 0.0)) return false;
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      B[i + (std::size_t) j * K] += w[i] * w[j] / ws;
  return true;
}

}}  // namespace dmod::trust_internal

#endif  // DMOD_TRUST_QN_H
