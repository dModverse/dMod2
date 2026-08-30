// C++ kernels for the Laplace / clustered marginal-likelihood backend (R: nlmeLaplace.R, L1Clustering.R).
//
// Two compute cores, each a bit-comparable port of the validated R reference
// (R/nlmeLaplace.R normalLaplace, R/L1Clustering.R .aghqLogI/.clusterParamMarginal):
//
//   normalLaplaceCpp    - the 1-D normal-Laplace marginal (convolution of a data
//                         Gaussian with a Laplace prior), its log-derivatives and
//                         posterior moments, vectorised. Numerically stable via
//                         log erfcx (= log(2) + x^2 + pnorm(-x sqrt2, log)) and a
//                         two-branch log-sum-exp, exactly as the R helper.
//
//   clusterMarginalCpp  - the exact clustered marginal for one parameter over G
//                         groups: the mode-centred Gauss-Hermite node loop that
//                         dominates the R cost (nq^(G-1) nodes). Given the anchor
//                         basis B, the change-of-variables factor L and the node
//                         set (U, logw, z2, sgn), it evaluates both the data and
//                         the prior-normaliser integrals in one sweep and returns
//                         the marginal -2 log L plus the posterior moments the
//                         Fisher outer gradient and the lambda M-step consume.
//
// The R functions keep their pure-R implementations as the reference; these
// kernels are selected behind a `useCpp` toggle and asserted equal to the R path
// to ~1e-10 in tests/testthat/test-EMlaplace.R.

#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif

static const double TWO_OVER_SQRTPI = 2.0 / std::sqrt(M_PI);

// log erfcx(x) = log(exp(x^2) erfc(x)), stable for all real x.
static inline double log_erfcx(double x) {
  return std::log(2.0) + x * x + R::pnorm(-x * std::sqrt(2.0), 0.0, 1.0, 1, 1);
}
// log(exp(la) + exp(lb)), overflow-safe.
static inline double log_sum_exp2(double la, double lb) {
  double m = std::max(la, lb);
  return m + std::log1p(std::exp(-std::fabs(la - lb)));
}
// d/dx log erfcx(x) = 2x - (2/sqrt(pi)) / erfcx(x).
static inline double d_log_erfcx(double x, double le) {
  return 2.0 * x - TWO_OVER_SQRTPI * std::exp(-le);
}

// [[Rcpp::export]]
List normalLaplaceCpp(NumericVector a, NumericVector m, NumericVector lambda) {
  int n = std::max({a.size(), m.size(), lambda.size()});
  NumericVector logI(n), dlogI_dm(n), dlogI_da(n), dlogI_dl(n),
                Eeta(n), Eabs(n), Ecen2(n), Ppos(n);
  for (int i = 0; i < n; ++i) {
    double ai = a[i % a.size()], mi = m[i % m.size()], li = lambda[i % lambda.size()];
    double cc = std::sqrt(ai / 2.0);
    double u  = (li / ai - mi) * cc;   // eta > 0 branch
    double v  = (li / ai + mi) * cc;   // eta < 0 branch
    double lEu = log_erfcx(u), lEv = log_erfcx(v);
    double lse = log_sum_exp2(lEu, lEv);
    double base = std::log(li) - std::log(2.0)
                + 0.5 * (std::log(M_PI) - std::log(2.0 * ai)) - ai * mi * mi / 2.0;
    double wpos = std::exp(lEu - lse), wneg = std::exp(lEv - lse);
    double dEu = d_log_erfcx(u, lEu), dEv = d_log_erfcx(v, lEv);
    double dc_da = cc / (2.0 * ai);
    double du_dm = -cc,               dv_dm = cc;
    double du_dl = cc / ai,           dv_dl = cc / ai;
    double du_da = (-li / (ai * ai)) * cc + (li / ai - mi) * dc_da;
    double dv_da = (-li / (ai * ai)) * cc + (li / ai + mi) * dc_da;
    double dB_dm = -ai * mi;
    double dB_dl = 1.0 / li;
    double dB_da = -1.0 / (2.0 * ai) - mi * mi / 2.0;
    double gm = dB_dm + wpos * dEu * du_dm + wneg * dEv * dv_dm;
    double gl = dB_dl + wpos * dEu * du_dl + wneg * dEv * dv_dl;
    double ga = dB_da + wpos * dEu * du_da + wneg * dEv * dv_da;
    logI[i]     = base + lse;
    dlogI_dm[i] = gm; dlogI_da[i] = ga; dlogI_dl[i] = gl;
    Eeta[i]  = mi + gm / ai;
    Eabs[i]  = 1.0 / li - gl;
    Ecen2[i] = -2.0 * ga;
    Ppos[i]  = wpos;
  }
  return List::create(_["logI"] = logI, _["dlogI_dm"] = dlogI_dm,
    _["dlogI_da"] = dlogI_da, _["dlogI_dlambda"] = dlogI_dl,
    _["Eeta"] = Eeta, _["Eabs"] = Eabs, _["Ecen2"] = Ecen2, _["Ppos"] = Ppos);
}

// weighted complete-graph total variation sum_{i<j} sizes_i sizes_j |e_i - e_j|.
static inline double weighted_tv(const double* e, const double* sizes, int G) {
  double s = 0.0;
  for (int i = 0; i < G; ++i)
    for (int j = i + 1; j < G; ++j)
      s += sizes[i] * sizes[j] * std::fabs(e[i] - e[j]);
  return s;
}

// Exact clustered marginal for one parameter over G groups via the dense/adaptive
// Gauss-Hermite node loop. B (G x d), L (d x d, lower chol of the local cov), zc
// (d, data-mode centre), U (M x d node coords), logw/z2/sgn (M). Evaluates the
// data integral (centred at zc) and the prior normaliser (centred at 0) over the
// same nodes, and returns value = -2 logI_data + 2 logI_prior plus the posterior
// moments. Mirrors .aghqLogI + .clusterParamMarginal (all-positive weights path).
//
// The M per-node change-of-variables maps (each a tiny matrix-vector product) are
// batched into three BLAS dgemm calls over all nodes at once: LU = sqrt(2) L U^T
// (d x M), then the group values E = B (zc + LU) and EZ = B LU (G x M), which is
// where a level-3 BLAS is worth its call overhead. The remaining per-node work
// (the O(G^2) fusion TV, the log-sum-exp) has no useful matrix structure and stays
// a scalar loop.
// [[Rcpp::export]]
List clusterMarginalCpp(NumericVector H, NumericVector m, double lambda,
                        NumericVector sizes, NumericMatrix B, NumericMatrix L,
                        NumericVector zc, NumericMatrix U,
                        NumericVector logw, NumericVector z2, NumericVector sgn) {
  int G = B.nrow(), d = B.ncol(), M = U.nrow();
  const double s2 = std::sqrt(2.0), one = 1.0, zero = 0.0;

  // LU (d x M) = sqrt(2) * L (d x d) * U^T (d x M)
  std::vector<double> LU(static_cast<size_t>(d) * M);
  F77_CALL(dgemm)("N", "T", &d, &M, &d, &s2, L.begin(), &d, U.begin(), &M,
                  &zero, LU.data(), &d FCONE FCONE);
  // ZI (d x M) = zc (broadcast) + LU
  std::vector<double> ZI(static_cast<size_t>(d) * M);
  for (int q = 0; q < M; ++q)
    for (int i = 0; i < d; ++i) ZI[i + static_cast<size_t>(q) * d] = zc[i] + LU[i + static_cast<size_t>(q) * d];
  // EI = B * ZI  and  EZ = B * LU   (both G x M)
  std::vector<double> EI(static_cast<size_t>(G) * M), EZ(static_cast<size_t>(G) * M);
  F77_CALL(dgemm)("N", "N", &G, &M, &d, &one, B.begin(), &G, ZI.data(), &d,
                  &zero, EI.data(), &G FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &G, &M, &d, &one, B.begin(), &G, LU.data(), &d,
                  &zero, EZ.data(), &G FCONE FCONE);

  std::vector<double> ltI(M), ltZ(M);
  for (int q = 0; q < M; ++q) {
    const double* eI = &EI[static_cast<size_t>(q) * G];
    const double* eZ = &EZ[static_cast<size_t>(q) * G];
    double quad = 0.0;
    for (int g = 0; g < G; ++g) { double r = eI[g] - m[g]; quad += H[g] * r * r; }
    double phiI = 0.25 * quad + 0.5 * lambda * weighted_tv(eI, sizes.begin(), G);
    ltI[q] = logw[q] + z2[q] - phiI;
    double phiZ = 0.5 * lambda * weighted_tv(eZ, sizes.begin(), G);
    ltZ[q] = logw[q] + z2[q] - phiZ;
  }
  // signed log-sum-exp for both integrals
  double mxI = *std::max_element(ltI.begin(), ltI.end());
  double mxZ = *std::max_element(ltZ.begin(), ltZ.end());
  double SI = 0.0, SZ = 0.0;
  for (int q = 0; q < M; ++q) {
    SI += sgn[q] * std::exp(ltI[q] - mxI);
    SZ += sgn[q] * std::exp(ltZ[q] - mxZ);
  }
  double value;
  if (SI > 0.0 && SZ > 0.0)
    value = -2.0 * (mxI + std::log(SI)) + 2.0 * (mxZ + std::log(SZ));
  else
    value = R_NaN;                          // caller falls back to the R path
  // posterior moments from the data integral
  NumericVector Eeta(G), Ecen2(G);
  double Efus = 0.0;
  if (SI > 0.0) {
    for (int q = 0; q < M; ++q) {
      double post = sgn[q] * std::exp(ltI[q] - mxI) / SI;
      const double* eq = &EI[static_cast<size_t>(q) * G];
      for (int g = 0; g < G; ++g) {
        Eeta[g]  += post * eq[g];
        Ecen2[g] += post * (eq[g] - m[g]) * (eq[g] - m[g]);
      }
      Efus += post * weighted_tv(eq, sizes.begin(), G);
    }
  }
  return List::create(_["value"] = value, _["Eeta"] = Eeta,
                      _["Ecen2"] = Ecen2, _["Efus"] = Efus);
}
