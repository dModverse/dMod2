// The steady states of Pimpl(): G(x, p) = [f(x, p); C(p) x - T(p)] = 0 for the dependent states x,
// the conserved amounts C(p) x with coefficients that may hold volumes.
//
// A root is solved by cppDE's pseudo-transient continuation from fixed starts in a fixed order:
// the start the parameters give, the same with a tenfold tighter flow error (flow only), all
// states at 1 and the same tighter, then random starts drawn from a generator seeded by the
// seed, the parameters and the start. States whose producing fluxes vanish at the parameters
// and whose rate is negative rest at 0 and are not solved for. The first start that converges
// gives the root; with flow and stability, a stable one.
//
// Derivatives by the implicit function theorem on the stacked system [f_x; C] dx = -G_p,
// scaled by rows and by the states, through its singular value decomposition: the
// pseudo-inverse where it is rank-deficient. Second order from
//   [f_x; C] d2x[a, b] = -(dx_a' G_xx dx_b + G_xp[, a] dx_b + G_xp[, b] dx_a + G_pp[a, b]).
//
// A root depends on the parameters, the start and the controls only, and so do its derivatives:
// repeated values are answered from memory, the same for every thread count.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "pimpl_kernel.h"

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

namespace dmod {
namespace pimpl {

using EvalFn = void (*)(double*, double*, double*, int*, int*, int*);
using EvalAd = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*);
using EvalAd2 = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, double*, double*,
                         double*, int*, int*, int*, int*, int*);
using PtcFn = int (*)(EvalAd, int, int, int, double*, const double*, int, const int*, const int*, int,
                      const double*, const double*, const double*, int, int*, char*, int);

struct Ctrl {
  bool keepRoot = true;
  double reltol = 1e-10, abstol = 1e-14, flowTol = 1, dtInit = 1e-2;
  int maxit = 0;  // 0: ceiling(70 log(n + 1))
  bool positive = true, stability = true;
  double stabilityTol = 0;
  int nStarts = 20;
  double lo = -5, hi = 5, seed = 1;
  bool log10 = true;
};

const std::size_t kMemo = 64;

struct Entry {
  Root r;
  std::string how;
};

struct State {
  EvalFn val = nullptr, influx = nullptr;
  EvalAd ad = nullptr;
  EvalAd2 ad2 = nullptr;
  PtcFn ptc = nullptr;
  int nDep = 0, nPar = 0, nEq = 0, k = 0;
  bool flow = true;
  std::vector<int> inDep;     // the dependent state of each influx output
  std::vector<int> inPar;     // the parameters of the influx, indices into pv
  std::vector<std::string> depNames;
  Ctrl ctrl;
  std::unordered_map<std::string, Entry> memo;
  std::deque<std::string> lru;
  int calls = 0, memoHits = 0, solves = 0, failed = 0;
  double iter = 0;
  std::map<std::string, int> how;
};

namespace {

// One solve: its root, how it was found, its iterations, the reason where it failed.
struct Solved {
  Root r;
  std::string how, error;
  std::vector<std::string> warnings;
  int iter = 0;
  bool ok = false;
};

// G at the rows of `x` (n x nDep, column-major).
bool values(const State& s, std::vector<double> x, int n, const std::vector<double>& pv, std::vector<double>& G) {
  std::vector<double> p(pv);
  if (p.empty()) p.push_back(0.0);
  G.assign((std::size_t) n * s.nEq, 0.0);
  int nn = n, nv = s.nDep, no = s.nEq;
  s.val(x.data(), G.data(), p.data(), &nn, &nv, &no);
  for (double v : G) if (!std::isfinite(v)) return false;
  return true;
}

// The coefficients C (k x nDep, column-major) of the conservation rows, linear in the states with
// coefficients that may depend on the parameters: their Jacobian at x = 0.
bool conservation(const State& s, const std::vector<double>& pv, std::vector<double>& C) {
  const int n = s.nDep, np = s.nPar, m = s.nEq, k = s.k;
  C.assign((std::size_t) k * n, 0.0);
  if (k == 0 || n == 0) return true;
  std::vector<double> x(n, 0.0), p(pv), dX((std::size_t) n * n, 0.0), dP((std::size_t) std::max(np, 1) * n, 0.0);
  for (int j = 0; j < n; ++j) dX[j + (std::size_t) n * j] = 1.0;
  if (p.empty()) p.push_back(0.0);
  std::vector<double> y(m), dy((std::size_t) m * n);
  int one = 1, nv = n, npp = np, no = m, nd = n;
  s.ad(x.data(), p.data(), dX.data(), dP.data(), y.data(), dy.data(), &one, &nv, &npp, &no, &nd);
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < n; ++j) {
      const double v = dy[(n + i) + (std::size_t) m * j];
      if (!std::isfinite(v)) return false;
      C[i + (std::size_t) k * j] = v;
    }
  return true;
}

// States whose producing fluxes vanish and whose rate is negative at two fixed positive probes
// rest at 0, repeated with the states found set to 0 until no further one is found.
std::vector<char> zeroStates(const State& s, const std::vector<double>& pv) {
  std::vector<char> zero(s.nDep, 0);
  if (s.influx == nullptr || s.inDep.empty()) return zero;
  const int n = s.nDep, ni = s.inDep.size();
  std::vector<double> X0(2 * (std::size_t) n), pin(s.inPar.size());
  for (int j = 0; j < n; ++j) {
    X0[0 + 2 * (std::size_t) j] = std::exp(std::sin(j + 1.0));
    X0[1 + 2 * (std::size_t) j] = std::exp(std::cos(j + 1.0));
  }
  for (std::size_t q = 0; q < s.inPar.size(); ++q) pin[q] = pv[s.inPar[q]];
  if (pin.empty()) pin.push_back(0.0);
  for (;;) {
    std::vector<double> X(X0), F(2 * (std::size_t) ni), G;
    for (int j = 0; j < n; ++j)
      if (zero[j]) X[0 + 2 * (std::size_t) j] = X[1 + 2 * (std::size_t) j] = 0.0;
    int two = 2, nv = n, no = ni;
    std::vector<double> Xi(X);
    s.influx(Xi.data(), F.data(), pin.data(), &two, &nv, &no);
    if (!values(s, X, 2, pv, G)) return zero;
    bool added = false;
    for (int i = 0; i < ni; ++i) {
      const int j = s.inDep[i];
      if (zero[j]) continue;
      if (F[0 + 2 * (std::size_t) i] == 0 && F[1 + 2 * (std::size_t) i] == 0 &&
          G[0 + 2 * (std::size_t) j] < 0 && G[1 + 2 * (std::size_t) j] < 0) {
        zero[j] = 1;
        added = true;
      }
    }
    if (!added) return zero;
  }
}

// A deterministic stream of uniforms on [0, 1) from the bytes of `key`.
struct Stream {
  std::mt19937_64 g;
  explicit Stream(const std::string& key) {
    std::uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : key) { h ^= c; h *= 1099511628211ULL; }
    g.seed(h);
  }
  double next() { return (g() >> 11) * (1.0 / 9007199254740992.0); }
};

std::string bytes(const std::vector<double>& a, const std::vector<double>& b) {
  std::string s(sizeof(double) * (a.size() + b.size()), '\0');
  if (!a.empty()) std::memcpy(&s[0], a.data(), sizeof(double) * a.size());
  if (!b.empty()) std::memcpy(&s[sizeof(double) * a.size()], b.data(), sizeof(double) * b.size());
  return s;
}

// A X = B [m x n, m x q] in the least-squares sense with X = diag(cs) Y, rows scaled to unit
// absolute sum; the pseudo-inverse where A is rank-deficient, with a warning naming the null
// directions. False where the decomposition fails.
bool solveScaled(const std::vector<double>& A0, int m, int n, const std::vector<double>& B0, int q,
                 const std::vector<double>& cs, const std::vector<std::string>& names, std::vector<double>& X,
                 std::vector<std::string>& warnings) {
  X.assign((std::size_t) n * q, 0.0);
  if (n == 0 || q == 0) return true;
  if (m < n) return false;
  std::vector<double> A(A0), B(B0), rs(m, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      A[i + (std::size_t) m * j] *= cs[j];
      rs[i] += std::fabs(A[i + (std::size_t) m * j]);
    }
  for (int i = 0; i < m; ++i) {
    rs[i] = std::max(rs[i], std::numeric_limits<double>::min());
    for (int j = 0; j < n; ++j) A[i + (std::size_t) m * j] /= rs[i];
    for (int c = 0; c < q; ++c) B[i + (std::size_t) m * c] /= rs[i];
  }
  const int mn = std::min(m, n);
  std::vector<double> sv(mn), U((std::size_t) m * mn), VT((std::size_t) mn * n), Aw(A);
  int lwork = -1, info = 0;
  double wq = 0;
  F77_CALL(dgesvd)("S", "S", &m, &n, Aw.data(), &m, sv.data(), U.data(), &m, VT.data(), &mn, &wq, &lwork, &info
                   FCONE FCONE);
  lwork = std::max(1, (int) wq);
  std::vector<double> w(lwork);
  F77_CALL(dgesvd)("S", "S", &m, &n, Aw.data(), &m, sv.data(), U.data(), &m, VT.data(), &mn, w.data(), &lwork,
                   &info FCONE FCONE);
  if (info != 0) return false;
  const double tol = std::max(m, n) * sv[0] * std::numeric_limits<double>::epsilon();
  int rank = 0;
  for (double d : sv) if (d > tol) ++rank;
  if (rank < n) {
    std::string msg;
    char buf[200];
    std::snprintf(buf, sizeof buf, "df/dx is rank-deficient (rank %d of %d; smallest sv %.2e of largest %.2e); ",
                  rank, n, sv[mn - 1], sv[0]);
    msg += buf;
    msg += "using SVD pseudoinverse (minimum-norm sensitivity on the constraint manifold).\n";
    msg += "Null-space direction(s) in dependent-state coordinates:";
    for (int j = rank; j < n; ++j) {
      double vmax = 0;
      for (int t = 0; t < n; ++t) vmax = std::max(vmax, std::fabs(VT[j + (std::size_t) mn * t]));
      std::snprintf(buf, sizeof buf, "\n  null #%d:", j - rank + 1);
      msg += buf;
      for (int t = 0; t < n; ++t) {
        const double v = VT[j + (std::size_t) mn * t];
        if (std::fabs(v) <= 0.05 * vmax && std::fabs(v) != vmax) continue;
        std::snprintf(buf, sizeof buf, " %+.3f*%s", v, names[t].c_str());
        msg += buf;
      }
    }
    msg += "\nLikely cause: an unmodelled conserved quantity, redundant equations, or a continuum of steady "
           "states. Pass an `eqnlist` so `Pimpl` can auto-detect CQs for cleaner sensitivities.";
    warnings.push_back(msg);
  }
  // Y = V diag(1/s) U' B over the singular values kept
  std::vector<double> UB((std::size_t) mn * q, 0.0);
  for (int c = 0; c < q; ++c)
    for (int a = 0; a < rank; ++a) {
      double s = 0;
      for (int i = 0; i < m; ++i) s += U[i + (std::size_t) m * a] * B[i + (std::size_t) m * c];
      UB[a + (std::size_t) mn * c] = s / sv[a];
    }
  for (int c = 0; c < q; ++c)
    for (int j = 0; j < n; ++j) {
      double s = 0;
      for (int a = 0; a < rank; ++a) s += VT[a + (std::size_t) mn * j] * UB[a + (std::size_t) mn * c];
      X[j + (std::size_t) n * c] = cs[j] * s;
      if (!std::isfinite(X[j + (std::size_t) n * c])) return false;
    }
  return true;
}

// The derivatives of the root r.x up to `order`. False where the model or the decomposition fails.
bool derivatives(const State& s, const std::vector<double>& pv, int order, Root& r,
                 std::vector<std::string>& warnings) {
  if (order < 1 || (r.hasD && (order < 2 || r.hasD2))) return true;
  const int n = s.nDep, np = s.nPar, m = s.nEq, nd = n + np;
  std::vector<double> x(r.x), p(pv), dX((std::size_t) std::max(n, 1) * nd, 0.0), dP((std::size_t) std::max(np, 1) * nd, 0.0);
  for (int j = 0; j < n; ++j) dX[j + (std::size_t) n * j] = 1.0;
  for (int a = 0; a < np; ++a) dP[a + (std::size_t) np * (n + a)] = 1.0;
  if (x.empty()) x.push_back(0.0);
  if (p.empty()) p.push_back(0.0);
  std::vector<double> y(m), dy((std::size_t) m * nd), d2y;
  int one = 1, nv = n, npp = np, no = m, ndd = nd;
  const bool second = order >= 2;
  if (second) {
    if (s.ad2 == nullptr) return false;
    d2y.assign((std::size_t) m * nd * nd, 0.0);
    int hx = 0, hp = 0;
    double none = 0;
    s.ad2(x.data(), p.data(), dX.data(), dP.data(), &none, &none, &hx, &hp, y.data(), dy.data(), d2y.data(), &one,
          &nv, &npp, &no, &ndd);
  } else {
    s.ad(x.data(), p.data(), dX.data(), dP.data(), y.data(), dy.data(), &one, &nv, &npp, &no, &ndd);
  }
  for (double v : dy) if (!std::isfinite(v)) return false;
  std::vector<double> A((std::size_t) m * n), B((std::size_t) m * np), cs(n);
  for (int j = 0; j < n; ++j) {
    cs[j] = r.x[j] > 0 ? r.x[j] : 1.0;
    for (int i = 0; i < m; ++i) A[i + (std::size_t) m * j] = dy[i + (std::size_t) m * j];
  }
  for (int a = 0; a < np; ++a)
    for (int i = 0; i < m; ++i) B[i + (std::size_t) m * a] = -dy[i + (std::size_t) m * (n + a)];
  std::vector<double> D;
  if (!solveScaled(A, m, n, B, np, cs, s.depNames, D, warnings)) return false;
  r.D = D;
  r.hasD = true;
  if (!second) return true;
  for (double v : d2y) if (!std::isfinite(v)) return false;
  // RHS[i, a, b] = dx_a' G_xx dx_b + G_xp[, a] dx_b + G_xp[, b] dx_a + G_pp[a, b]
  auto H = [&](int i, int u, int v) { return d2y[i + (std::size_t) m * (u + (std::size_t) nd * v)]; };
  std::vector<double> R((std::size_t) m * np * np, 0.0), HD((std::size_t) n * np);
  for (int i = 0; i < m; ++i) {
    // HD[j, b] = sum_l G_xx[j, l] dx[l, b] + G_xp[j, b]
    for (int b = 0; b < np; ++b)
      for (int j = 0; j < n; ++j) {
        double v = H(i, j, n + b);
        for (int l = 0; l < n; ++l) v += H(i, j, l) * D[l + (std::size_t) n * b];
        HD[j + (std::size_t) n * b] = v;
      }
    for (int b = 0; b < np; ++b)
      for (int a = 0; a < np; ++a) {
        double v = H(i, n + a, n + b);
        for (int j = 0; j < n; ++j) v += D[j + (std::size_t) n * a] * HD[j + (std::size_t) n * b] +
                                         H(i, n + a, j) * D[j + (std::size_t) n * b];
        R[i + (std::size_t) m * (a + (std::size_t) np * b)] = -v;
      }
  }
  std::vector<double> D2;
  std::vector<std::string> again;
  if (!solveScaled(A, m, n, R, np * np, cs, s.depNames, D2, again)) return false;
  r.D2 = D2;
  r.hasD2 = true;
  return true;
}

// The eigenvalues `re` + i `im` of J = df/dx at the root x on ker C, and for each one whose real
// part exceeds `onset` (1-based in `active`) the gradient of that real part in the parameters,
// Re w^H (J_p + J_x D) u / (w^H u), u and w its right and left eigenvectors in J and D = dx/dp.
bool spectrum(const State& s, const std::vector<double>& pv, const std::vector<double>& x, double onset,
              std::vector<double>& re, std::vector<double>& im, std::vector<int>& active,
              std::vector<double>& grad, std::string& error) {
  using cd = std::complex<double>;
  const int n = s.nDep, np = s.nPar, m = s.nEq, k = s.k, nd = n + np, q = n - k;
  if (s.ad2 == nullptr) {
    error = "the stability of the steady state needs its second derivatives, Pimpl(deriv2 = TRUE)";
    return false;
  }
  Root r;
  r.x = x;
  std::vector<std::string> warn;
  if (np > 0 && !derivatives(s, pv, 1, r, warn)) {
    error = "the derivatives of the steady state failed";
    return false;
  }
  std::vector<double> xx(x), p(pv), dX((std::size_t) std::max(n, 1) * nd, 0.0),
      dP((std::size_t) std::max(np, 1) * nd, 0.0);
  for (int j = 0; j < n; ++j) dX[j + (std::size_t) n * j] = 1.0;
  for (int a = 0; a < np; ++a) dP[a + (std::size_t) np * (n + a)] = 1.0;
  if (xx.empty()) xx.push_back(0.0);
  if (p.empty()) p.push_back(0.0);
  std::vector<double> y(m), dy((std::size_t) m * nd), d2y((std::size_t) m * nd * nd, 0.0);
  int one = 1, nv = n, npp = np, no = m, ndd = nd, hx = 0, hp = 0;
  double none = 0;
  s.ad2(xx.data(), p.data(), dX.data(), dP.data(), &none, &none, &hx, &hp, y.data(), dy.data(), d2y.data(), &one,
        &nv, &npp, &no, &ndd);
  for (double v : dy)
    if (!std::isfinite(v)) {
      error = "the Jacobian is not finite at the steady state";
      return false;
    }
  std::vector<double> J((std::size_t) n * n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) J[i + (std::size_t) n * j] = dy[i + (std::size_t) m * j];
  // an orthonormal basis N of ker C: the right singular vectors of C beyond its rank
  std::vector<double> N((std::size_t) n * std::max(q, 1), 0.0);
  if (k == 0) {
    for (int j = 0; j < n; ++j) N[j + (std::size_t) n * j] = 1.0;
  } else {
    std::vector<double> C((std::size_t) k * n), sv(std::min(k, n)), U(1), VT((std::size_t) n * n);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < k; ++i) C[i + (std::size_t) k * j] = dy[(n + i) + (std::size_t) m * j];
    int kk = k, nn = n, ldu = 1, lwork = -1, info = 0;
    double wq = 0;
    F77_CALL(dgesvd)("N", "A", &kk, &nn, C.data(), &kk, sv.data(), U.data(), &ldu, VT.data(), &nn, &wq, &lwork,
                     &info FCONE FCONE);
    lwork = std::max(1, (int) wq);
    std::vector<double> w(lwork);
    F77_CALL(dgesvd)("N", "A", &kk, &nn, C.data(), &kk, sv.data(), U.data(), &ldu, VT.data(), &nn, w.data(),
                     &lwork, &info FCONE FCONE);
    if (info != 0) {
      error = "the decomposition of the conservation rows failed";
      return false;
    }
    for (int t = 0; t < q; ++t)
      for (int j = 0; j < n; ++j) N[j + (std::size_t) n * t] = VT[(k + t) + (std::size_t) n * j];
  }
  // eigen decomposition with left and right vectors; false where it fails
  auto eigen = [](std::vector<double> M, int d, bool vectors, std::vector<double>& wr, std::vector<double>& wi,
                  std::vector<double>& vl, std::vector<double>& vr) {
    wr.assign(d, 0.0);
    wi.assign(d, 0.0);
    vl.assign(vectors ? (std::size_t) d * d : 1, 0.0);
    vr.assign(vectors ? (std::size_t) d * d : 1, 0.0);
    if (d == 0) return true;
    const char* job = vectors ? "V" : "N";
    int dd = d, ldv = vectors ? d : 1, lwork = -1, info = 0;
    double wq = 0;
    F77_CALL(dgeev)(job, job, &dd, M.data(), &dd, wr.data(), wi.data(), vl.data(), &ldv, vr.data(), &ldv, &wq,
                    &lwork, &info FCONE FCONE);
    lwork = std::max(1, (int) wq);
    std::vector<double> w(lwork);
    F77_CALL(dgeev)(job, job, &dd, M.data(), &dd, wr.data(), wi.data(), vl.data(), &ldv, vr.data(), &ldv, w.data(),
                    &lwork, &info FCONE FCONE);
    return info == 0;
  };
  // the spectrum on ker C: of A = N' J N, since C J = 0 puts the range of J into ker C
  std::vector<double> A((std::size_t) q * q, 0.0);
  for (int b = 0; b < q; ++b)
    for (int a = 0; a < q; ++a) {
      double v = 0.0;
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          v += N[i + (std::size_t) n * a] * J[i + (std::size_t) n * j] * N[j + (std::size_t) n * b];
      A[a + (std::size_t) q * b] = v;
    }
  std::vector<double> vl, vr, wrJ, wiJ, vlJ, vrJ;
  if (!eigen(A, q, false, re, im, vl, vr) || !eigen(J, n, true, wrJ, wiJ, vlJ, vrJ)) {
    error = "the eigenvalues of the Jacobian at the steady state failed";
    return false;
  }
  auto H = [&](int i, int j, int l) { return d2y[i + (std::size_t) m * (j + (std::size_t) nd * l)]; };
  active.clear();
  grad.clear();
  for (int t = 0; t < q; ++t) {
    if (!(re[t] > onset)) continue;
    // the same eigenvalue of J, whose vectors the derivative takes
    int c = 0;
    double best = std::numeric_limits<double>::infinity();
    for (int e = 0; e < n; ++e) {
      const double dist = std::hypot(wrJ[e] - re[t], wiJ[e] - im[t]);
      if (dist < best) { best = dist; c = e; }
    }
    if (best > 1e-6 * (1.0 + std::hypot(re[t], im[t]))) {
      error = "an eigenvalue on ker C has no match in the Jacobian";
      return false;
    }
    std::vector<cd> u(n), w(n);
    for (int i = 0; i < n; ++i) {
      if (wiJ[c] == 0.0) {
        u[i] = vrJ[i + (std::size_t) n * c];
        w[i] = vlJ[i + (std::size_t) n * c];
      } else if (wiJ[c] > 0.0) {
        u[i] = cd(vrJ[i + (std::size_t) n * c], vrJ[i + (std::size_t) n * (c + 1)]);
        w[i] = cd(vlJ[i + (std::size_t) n * c], vlJ[i + (std::size_t) n * (c + 1)]);
      } else {
        u[i] = cd(vrJ[i + (std::size_t) n * (c - 1)], -vrJ[i + (std::size_t) n * c]);
        w[i] = cd(vlJ[i + (std::size_t) n * (c - 1)], -vlJ[i + (std::size_t) n * c]);
      }
    }
    cd den = 0.0;
    for (int i = 0; i < n; ++i) den += std::conj(w[i]) * u[i];
    // Q_l = w^H (dJ/dz_l) u for every variable and parameter z_l
    std::vector<cd> Q(nd, 0.0);
    for (int l = 0; l < nd; ++l)
      for (int j = 0; j < n; ++j) {
        if (u[j] == cd(0.0)) continue;
        cd v = 0.0;
        for (int i = 0; i < n; ++i) v += std::conj(w[i]) * H(i, j, l);
        Q[l] += v * u[j];
      }
    active.push_back(t + 1);
    for (int a = 0; a < np; ++a) {
      cd g = Q[n + a];
      for (int l = 0; l < n; ++l) g += Q[l] * r.D[l + (std::size_t) n * a];
      grad.push_back(std::real(g / den));
    }
  }
  return true;
}

// The root from the starts in their order; the iterations, how, or the reasons of every start.
Solved solveOne(const State& s, const Ctrl& c, const std::vector<double>& pv, const std::vector<double>& guess,
                int order) {
  Solved out;
  const int n = s.nDep, k = s.k;
  std::vector<char> zero = zeroStates(s, pv);
  std::vector<int> act;
  for (int j = 0; j < n; ++j) if (!zero[j]) act.push_back(j);
  out.r.x.assign(n, 0.0);
  if (act.empty()) {
    out.ok = true;
    out.how = "all zero";
    if (!derivatives(s, pv, order, out.r, out.warnings)) out.r.hasD = out.r.hasD2 = false;
    return out;
  }
  // the totals the dependent states hold: -G at x = 0 on the conservation rows
  std::vector<double> G0, T(k);
  if (!values(s, std::vector<double>(std::max(n, 1), 0.0), 1, pv, G0)) {
    out.error = "the system is not finite at the parameters";
    return out;
  }
  std::vector<double> C;
  if (!conservation(s, pv, C)) {
    out.error = "the conserved quantities are not finite at the parameters";
    return out;
  }
  const int na = act.size();
  std::vector<double> Ca((std::size_t) k * na);
  for (int i = 0; i < k; ++i) {
    T[i] = -G0[n + i];
    bool held = false;
    for (int t = 0; t < na; ++t) {
      Ca[i + (std::size_t) k * t] = C[i + (std::size_t) k * act[t]];
      held = held || Ca[i + (std::size_t) k * t] != 0;
    }
    if (T[i] > 0 && !held) {
      out.error = "Pimpl: a conserved total is positive but every species holding it rests at 0.";
      return out;
    }
  }
  auto prep = [&](std::vector<double> x) {
    for (int j = 0; j < n; ++j) if (zero[j]) x[j] = 0;
    std::vector<double> xa(na);
    for (int t = 0; t < na; ++t) {
      xa[t] = x[act[t]];
      if (c.positive && !(xa[t] > 0)) xa[t] = 1;
    }
    if (k > 0) onManifold(xa.data(), na, Ca.data(), k, T.data());
    for (int t = 0; t < na; ++t) x[act[t]] = xa[t];
    return x;
  };
  struct Start { std::string name; std::vector<double> x; bool fine; };
  std::vector<Start> starts;
  starts.push_back({"guess", guess, false});
  if (s.flow) starts.push_back({"guess_fine", guess, true});
  bool ones = false;
  for (int t = 0; t < na; ++t) ones = ones || guess[act[t]] != 1;
  if (ones) {
    starts.push_back({"ones", std::vector<double>(n, 1.0), false});
    if (s.flow) starts.push_back({"ones_fine", std::vector<double>(n, 1.0), true});
  }
  if (c.nStarts > 0) {
    std::vector<double> key(pv);
    key.push_back(c.seed);
    Stream rng(bytes(key, guess));
    for (int i = 0; i < c.nStarts; ++i) {
      std::vector<double> x(n);
      for (int j = 0; j < n; ++j) {
        const double u = c.lo + (c.hi - c.lo) * rng.next();
        x[j] = c.log10 ? std::pow(10.0, u) : u;
      }
      starts.push_back({"start" + std::to_string(i + 1), x, false});
    }
  }
  const int maxit = c.maxit > 0 ? c.maxit : (int) std::ceil(70 * std::log(n + 1.0));
  std::string log;
  for (const Start& st : starts) {
    std::vector<double> x = prep(st.x);
    const double ctrl[10] = {c.reltol, c.abstol, c.flowTol / (st.fine ? 10 : 1), c.dtInit, 5, (double) maxit,
                             s.flow ? 1.0 : 0.0, c.positive ? 1.0 : 0.0, c.stability ? 1.0 : 0.0, c.stabilityTol};
    int iterations = 0;
    char msg[256];
    const int ok = s.ptc(s.ad, n, s.nPar, s.nEq, x.data(), pv.empty() ? nullptr : pv.data(), na, act.data(),
                         act.data(), k, Ca.data(), T.data(), ctrl, 10, &iterations, msg, sizeof msg);
    out.iter += iterations;
    log += "\n  " + st.name + ": " + msg;
    if (!ok) continue;
    for (int j = 0; j < n; ++j) out.r.x[j] = zero[j] ? 0.0 : x[j];
    out.ok = true;
    out.how = st.name + (std::strstr(msg, "unstable") != nullptr ? "+escape" : "");
    if (!derivatives(s, pv, order, out.r, out.warnings)) {
      out.r.hasD = out.r.hasD2 = false;
      out.warnings.push_back("Pimpl: IFT-based sensitivities unavailable at the current root. Returning value only.");
    }
    return out;
  }
  out.error = std::string("Pimpl: no ") + (s.flow ? "stable " : "") + "root found. Attempts:" + log;
  return out;
}

void remember(State& s, const std::string& key, const Root& r, const std::string& how) {
  auto it = s.memo.find(key);
  if (it == s.memo.end()) {
    s.memo.emplace(key, Entry{r, how});
    s.lru.push_back(key);
    if (s.lru.size() > kMemo) {
      s.memo.erase(s.lru.front());
      s.lru.pop_front();
    }
  } else {
    it->second.r = r;
  }
}

double num(const List& l, const char* name, double def) {
  if (!l.containsElementNamed(name)) return def;
  SEXP v = l[name];
  if (Rf_isNull(v) || Rf_length(v) == 0) return def;
  return Rf_asReal(v);
}

}  // namespace

bool roots(const std::vector<Request>& req, std::vector<Root>& out, int threads, std::string& error,
           std::vector<std::string>& warnings) {
  const int N = req.size();
  out.assign(N, Root());
  // per request: a memory entry, the task of another request, or a task of its own
  std::vector<int> task(N, -1);
  std::vector<std::string> keys(N);
  std::vector<char> hit(N, 0);
  std::vector<Entry> hitEntry(N);
  struct Task { int r; bool fromMemory; };
  std::vector<Task> tasks;
  std::map<std::pair<State*, std::string>, int> pending;
  for (int i = 0; i < N; ++i) {
    State& s = *req[i].st;
    ++s.calls;
    keys[i] = bytes(req[i].pv, req[i].guess);
    auto at = pending.find({req[i].st, keys[i]});
    if (at != pending.end()) { task[i] = at->second; continue; }
    if (s.ctrl.keepRoot) {
      auto m = s.memo.find(keys[i]);
      if (m != s.memo.end()) {
        hit[i] = 1;
        hitEntry[i] = m->second;
        const Root& r = m->second.r;
        if (req[i].order < 1 || (r.hasD && (req[i].order < 2 || r.hasD2))) continue;
      }
    }
    task[i] = tasks.size();
    pending[{req[i].st, keys[i]}] = tasks.size();
    tasks.push_back({i, hit[i] != 0});
  }
  // the highest order any request of a task asks for
  std::vector<int> order(tasks.size(), 0);
  for (int i = 0; i < N; ++i) if (task[i] >= 0) order[task[i]] = std::max(order[task[i]], req[i].order);
  std::vector<Solved> res(tasks.size());
  const int T = tasks.size();
  int nt = std::max(1, std::min(threads, T));
#ifdef _OPENMP
  if (omp_in_parallel()) nt = 1;
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int t = 0; t < T; ++t) {
    const Request& q = req[tasks[t].r];
    if (tasks[t].fromMemory) {
      Solved& o = res[t];
      o.r = hitEntry[tasks[t].r].r;
      o.how = hitEntry[tasks[t].r].how;
      o.ok = true;
      if (!derivatives(*q.st, q.pv, order[t], o.r, o.warnings)) {
        o.r.hasD = o.r.hasD2 = false;
        o.warnings.push_back("Pimpl: IFT-based sensitivities unavailable at the current root. Returning value only.");
      }
      continue;
    }
    res[t] = solveOne(*q.st, q.st->ctrl, q.pv, q.guess, order[t]);
  }
  // in the order of the requests: memory, statistics, the first failure
  std::vector<char> counted(T, 0);
  for (int i = 0; i < N; ++i) {
    State& s = *req[i].st;
    if (task[i] < 0) {
      ++s.memoHits;
      out[i] = hitEntry[i].r;
      continue;
    }
    const int t = task[i];
    const Solved& o = res[t];
    const bool owner = tasks[t].r == i;
    if (owner && !counted[t]) {
      counted[t] = 1;
      for (const std::string& w : o.warnings)
        if (std::find(warnings.begin(), warnings.end(), w) == warnings.end()) warnings.push_back(w);
      if (tasks[t].fromMemory) ++s.memoHits;
      else {
        ++s.solves;
        s.iter += o.iter;
        if (o.ok) ++s.how[o.how];
        else ++s.failed;
      }
      if (!o.ok) { error = o.error; return false; }
      if (s.ctrl.keepRoot) remember(s, keys[i], o.r, o.how);
    } else {
      if (s.ctrl.keepRoot) ++s.memoHits;
      else ++s.solves;
    }
    out[i] = o.r;
  }
  return true;
}

State* stateOf(SEXP xp) {
  if (TYPEOF(xp) != EXTPTRSXP) return nullptr;
  return static_cast<State*>(R_ExternalPtrAddr(xp));
}

}  // namespace pimpl
}  // namespace dmod

using dmod::pimpl::State;

namespace {

void finalizeState(SEXP xp) {
  delete static_cast<State*>(R_ExternalPtrAddr(xp));
  R_ClearExternalPtr(xp);
}

State* need(SEXP xp) {
  State* s = dmod::pimpl::stateOf(xp);
  if (s == nullptr) stop("Pimpl: its native state is gone");
  return s;
}

DL_FUNC fnOf(SEXP x) {
  if (TYPEOF(x) != EXTPTRSXP) return nullptr;
  return R_ExternalPtrAddrFn(x);
}

}  // namespace

// The native state of a Pimpl(): the entries of its system (`val`, `ad`, `ad2`, `influx`),
// cppDE's solver (`ptc`), the sizes, the number of conservation rows `nTotals`, per influx output
// its dependent state `inDep` and the parameters it reads `inPar` (0-based), and `flow`.
// [[Rcpp::export]]
SEXP pimpl_state(List spec) {
  State* s = new State();
  s->val = reinterpret_cast<dmod::pimpl::EvalFn>(fnOf(spec["val"]));
  s->ad = reinterpret_cast<dmod::pimpl::EvalAd>(fnOf(spec["ad"]));
  s->ad2 = reinterpret_cast<dmod::pimpl::EvalAd2>(fnOf(spec["ad2"]));
  s->influx = reinterpret_cast<dmod::pimpl::EvalFn>(fnOf(spec["influx"]));
  s->ptc = reinterpret_cast<dmod::pimpl::PtcFn>(fnOf(spec["ptc"]));
  if (s->val == nullptr || s->ad == nullptr || s->ptc == nullptr) {
    delete s;
    stop("Pimpl: the compiled system or cppDE's solver is not available");
  }
  s->nDep = as<int>(spec["nDep"]);
  s->nPar = as<int>(spec["nPar"]);
  s->k = as<int>(spec["nTotals"]);
  s->nEq = s->nDep + s->k;
  s->inDep = as<std::vector<int>>(spec["inDep"]);
  s->inPar = as<std::vector<int>>(spec["inPar"]);
  s->depNames = as<std::vector<std::string>>(spec["depNames"]);
  s->flow = as<bool>(spec["flow"]);
  SEXP xp = PROTECT(R_MakeExternalPtr(s, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(xp, finalizeState, TRUE);
  UNPROTECT(1);
  return xp;
}

// Sets the controls of the state and forgets its roots.
// [[Rcpp::export]]
void pimpl_controls(SEXP xp, bool keepRoot, List ctrl) {
  State* s = need(xp);
  dmod::pimpl::Ctrl c;
  c.keepRoot = keepRoot;
  c.reltol = dmod::pimpl::num(ctrl, "reltol", c.reltol);
  c.abstol = dmod::pimpl::num(ctrl, "abstol", c.abstol);
  c.flowTol = dmod::pimpl::num(ctrl, "flowTol", c.flowTol);
  c.dtInit = dmod::pimpl::num(ctrl, "dtInit", c.dtInit);
  c.maxit = (int) dmod::pimpl::num(ctrl, "maxit", 0);
  c.positive = dmod::pimpl::num(ctrl, "positive", 1) != 0;
  c.stability = dmod::pimpl::num(ctrl, "stability", 1) != 0;
  c.stabilityTol = dmod::pimpl::num(ctrl, "stabilityTol", 0);
  c.nStarts = (int) dmod::pimpl::num(ctrl, "nStarts", c.nStarts);
  c.seed = dmod::pimpl::num(ctrl, "seed", c.seed);
  if (ctrl.containsElementNamed("startRange")) {
    NumericVector r = ctrl["startRange"];
    c.lo = r[0];
    c.hi = r[1];
  }
  if (ctrl.containsElementNamed("startScale"))
    c.log10 = as<std::string>(ctrl["startScale"]) == "log10";
  s->ctrl = c;
  s->memo.clear();
  s->lru.clear();
}

// The statistics of the state: calls, roots from memory, solves, iterations, how the roots were
// found and failed solves; with `reset` set to zero afterwards.
// [[Rcpp::export]]
List pimpl_stats(SEXP xp, bool reset = false) {
  State* s = need(xp);
  IntegerVector how(s->how.size());
  CharacterVector nm(s->how.size());
  int i = 0;
  for (const auto& h : s->how) { nm[i] = h.first; how[i] = h.second; ++i; }
  how.attr("names") = nm;
  List out = List::create(Named("calls") = s->calls, Named("memo") = s->memoHits, Named("solves") = s->solves,
                          Named("iter") = s->iter, Named("how") = how, Named("failed") = s->failed);
  if (reset) {
    s->calls = s->memoHits = s->solves = s->failed = 0;
    s->iter = 0;
    s->how.clear();
  }
  return out;
}

// The root at `pv` from `guess` with derivatives of order `order`: list(x, D, D2, warnings), D and
// D2 NULL where not asked for or not available, the warnings for R to report.
// [[Rcpp::export]]
List pimpl_root(SEXP xp, NumericVector pv, NumericVector guess, int order) {
  State* s = need(xp);
  dmod::pimpl::Request q;
  q.st = s;
  q.pv.assign(pv.begin(), pv.end());
  q.guess.assign(guess.begin(), guess.end());
  q.order = order;
  std::vector<dmod::pimpl::Root> out;
  std::string error;
  std::vector<std::string> warnings;
  const bool ok = dmod::pimpl::roots({q}, out, 1, error, warnings);
  if (!ok) stop(error);
  const dmod::pimpl::Root& r = out[0];
  const int n = s->nDep, np = s->nPar;
  // RObject keeps D protected while D2 and the result list are allocated
  RObject D = R_NilValue, D2 = R_NilValue;
  if (order >= 1 && r.hasD) {
    NumericMatrix d(n, np);
    std::copy(r.D.begin(), r.D.end(), d.begin());
    D = d;
  }
  if (order >= 2 && r.hasD2) {
    NumericVector d2(r.D2.begin(), r.D2.end());
    d2.attr("dim") = IntegerVector::create(n, np, np);
    D2 = d2;
  }
  return List::create(Named("x") = NumericVector(r.x.begin(), r.x.end()), Named("D") = D, Named("D2") = D2,
                      Named("warnings") = wrap(warnings));
}

// The spectrum of the Jacobian at the root `x` of `pv` on ker C: list(re, im, active, grad), `grad`
// the gradients in the parameters of the real parts of the eigenvalues `active` above `onset`.
// [[Rcpp::export]]
List pimpl_spectrum(SEXP xp, NumericVector pv, NumericVector x, double onset) {
  State* s = need(xp);
  std::vector<double> re, im, grad;
  std::vector<int> active;
  std::string error;
  if (!dmod::pimpl::spectrum(*s, as<std::vector<double>>(pv), as<std::vector<double>>(x), onset, re, im, active,
                             grad, error))
    stop("Pimpl: " + error);
  const int na = active.size(), np = s->nPar;
  NumericMatrix G(na, np);
  for (int t = 0; t < na; ++t)
    for (int a = 0; a < np; ++a) G(t, a) = grad[(std::size_t) t * np + a];
  return List::create(Named("re") = wrap(re), Named("im") = wrap(im), Named("active") = wrap(active),
                      Named("grad") = G);
}

// Whether the native state behind `xp` still exists; not after the session it was built in.
// [[Rcpp::export]]
bool pimpl_live(SEXP xp) {
  return dmod::pimpl::stateOf(xp) != nullptr;
}
