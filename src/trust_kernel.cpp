// trust-region optimizer.
//
// Two boundary schemes, selected by `boundary`:
//
//   "reflective"  Coleman-Li interior trust-region-reflective (see
//                 trust_subproblem.h). Iterates stay strictly inside the box and
//                 the bound acts as a change of metric rather than a projection,
//                 so the predicted and the achieved reduction always describe
//                 the same step.
//   "clip"        Moré-Sorensen on a reduced working space that drops
//                 coordinates pinned at a bound, followed by componentwise
//                 clipping. Kept so fits made before the reflective scheme
//                 remain reproducible; frozen, not developed further.

#include <Rcpp.h>
#include "trust_subproblem.h"
#include "trust_driver.h"
#include "trust_qn.h"
#include "trust_model.h"
#include <vector>
#include <deque>
#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>
#include <climits>

using namespace Rcpp;
using dmod::trust_internal::affine_scaling;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::model_value;
using dmod::trust_internal::qn_update_B;
using dmod::trust_internal::stepback;
using dmod::trust_internal::trust_sub;
using dmod::trust_driver::Blather;
using dmod::trust_driver::Reporter;
using dmod::trust_driver::CV_EXACT;
using dmod::trust_driver::CV_GN;
using dmod::trust_driver::CV_GRADIENT;
using dmod::trust_driver::eval_objfun;
using dmod::trust_driver::stop_parinit_failed;
using dmod::trust_driver::fill_bound;
using dmod::trust_driver::fill_parscale;
using dmod::trust_driver::kInf;
using dmod::trust_driver::kStallLimit;
using dmod::trust_driver::kReseedLimit;
using dmod::trust_driver::push_interior;
using dmod::trust_driver::read_hessian;
using dmod::trust_driver::subproblem_label;

namespace {

// Interchangeable source of the model Hessian. The quasi-Newton approximation
// lives in the working Hessian `H_full`, so the scaling, eigendecomposition
// and subproblem downstream are untouched.
enum HessianMethod { HM_NONE = -1, HM_GN = 0, HM_BFGS = 1, HM_SR1 = 2,
                     HM_EXACT = 3, HM_SBFGS = 4, HM_SSR1 = 5, HM_GNSBFGS = 6,
                     HM_GNSSR1 = 7 };

// Structured methods keep the Gauss-Newton Hessian of the objective and update
// only a correction A; the hybrid ones choose between both per accepted step.
inline bool hm_structured(int m) { return m >= HM_SBFGS; }
inline bool hm_hybrid(int m)     { return m == HM_GNSBFGS || m == HM_GNSSR1; }
inline bool hm_sr1_type(int m)   { return m == HM_SR1 || m == HM_SSR1 || m == HM_GNSSR1; }

// Seed of the quasi-Newton approximation. Identity is all a gradient-only
// derivative scheme can supply.
enum HessianInit { HI_GN = 0, HI_IDENTITY = 1, HI_EXACT = 2, HI_TRACKED = 3 };
enum HessianReseed { HR_NEVER = 0, HR_STALL = 1 };

// Quasi-Newton settings of one phase.
struct QnParams {
  int init;
  int reseed;
  int memory;
  double cautious;
  bool rejected;
  bool scaled;        // A scaled with the residual norm (Huschens)
  double hybrid_tol;  // threshold of the switch rule of the hybrid methods
  bool hybrid_shift;  // rho I added on the Gauss-Newton branch of the hybrid methods
};

// One phase of a run: the Hessian source, a separate iteration budget, the sweep
// the objective is evaluated with, and the quasi-Newton settings.
struct Phase {
  int method;
  int iterlim;
  int radius_stall;  // handover after this many iterations at an unchanged radius, 0 off
  bool reverse;
  bool keep_radius;  // the phase starts at the radius of the handover instead of rinit
  bool stop_ends;    // a value, model, step or stagnation stop ends the run instead of handing over
  QnParams qn;
};

// Same update, applied to the working Hessian `H` of `val` rather than of phi.
void qn_update(int kind, int K, std::vector<double>& H,
               const std::vector<double>& s, const std::vector<double>& y,
               double sgn) {
  const std::size_t KK = (std::size_t) K * K;
  std::vector<double> B(KK);
  for (std::size_t t = 0; t < KK; ++t) B[t] = sgn * H[t];
  qn_update_B(kind, K, B, s, y);
  for (std::size_t t = 0; t < KK; ++t) H[t] = sgn * B[t];
}

// Limited memory: rebuild `H` from gamma*I and the kept pairs, oldest first, so
// it depends on neither the seed nor curvature outside the window.
void qn_assemble(int kind, int K, std::vector<double>& H, double gamma,
                 const std::deque<std::vector<double>>& S,
                 const std::deque<std::vector<double>>& Y, double sgn) {
  const std::size_t KK = (std::size_t) K * K;
  std::vector<double> B(KK, 0.0);
  for (int i = 0; i < K; ++i) B[i + (std::size_t) i * K] = gamma;
  for (std::size_t k = 0; k < S.size(); ++k) qn_update_B(kind, K, B, S[k], Y[k]);
  for (std::size_t t = 0; t < KK; ++t) H[t] = sgn * B[t];
}

// A BFGS update keeps the inertia of the matrix it starts from, so an indefinite seed stays
// indefinite for the whole phase. Such a seed takes the absolute values of the eigenvalues, floored
// at the noise level of the decomposition; a positive semidefinite seed is left as it is.
void bfgs_seed_definite(int K, std::vector<double>& H, double sgn) {
  const std::size_t KK = (std::size_t) K * K;
  std::vector<double> B(KK), vals(K), vecs(KK);
  for (std::size_t t = 0; t < KK; ++t) B[t] = sgn * H[t];
  eigen_sym_local(B.data(), K, vals.data(), vecs.data());
  if (!(vals[0] < 0.0)) return;
  double vmax = 0.0;
  for (int j = 0; j < K; ++j) vmax = std::max(vmax, std::fabs(vals[j]));
  const double cut = 8.0 * K * std::numeric_limits<double>::epsilon() * vmax;
  std::fill(B.begin(), B.end(), 0.0);
  for (int k = 0; k < K; ++k) {
    const double lam = std::max(std::fabs(vals[k]), cut);
    const double* q = vecs.data() + (std::size_t) k * K;
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) B[i + (std::size_t) j * K] += lam * q[i] * q[j];
  }
  for (std::size_t t = 0; t < KK; ++t) H[t] = sgn * B[t];
}

// Li-Fukushima cautious update: a pair with negligible curvature along the step
// contains no information. SR1 is exempt, negative curvature is legal there.
bool qn_pair_informative(int K, const std::vector<double>& s,
                         const std::vector<double>& y, double cautious) {
  if (!(cautious > 0.0)) return true;
  double sy = 0.0, ss = 0.0, yy = 0.0;
  for (int i = 0; i < K; ++i) { sy += s[i] * y[i]; ss += s[i] * s[i]; yy += y[i] * y[i]; }
  return sy > cautious * std::sqrt(ss * yy);
}

// Pairs kept by limited-memory BFGS on the sparse model where qnMemory is 0.
const int kSparseMemory = 10;

// y = H x for the Hessian of the sparse model in the sign of the value: S from the upper triangle
// in compressed columns plus U diag(d) U'.
void sparse_hess_apply(int K, const dmod::trust_model::SparseModel::Hess& h, const double* x, double* y) {
  std::fill(y, y + K, 0.0);
  for (int j = 0; j < K; ++j)
    for (int k = h.p[j]; k < h.p[j + 1]; ++k) {
      const int i = h.i[k];
      if (i > j) continue;
      y[i] += h.x[k] * x[j];
      if (i != j) y[j] += h.x[k] * x[i];
    }
  const int m = (int) h.d.size();
  for (int l = 0; l < m; ++l) {
    const double* Ul = h.U.data() + (std::size_t) l * K;
    double t = 0.0;
    for (int i = 0; i < K; ++i) t += Ul[i] * x[i];
    t *= h.d[l];
    for (int i = 0; i < K; ++i) y[i] += t * Ul[i];
  }
}

// A BFGS update keeps the inertia of the matrix it starts from. A sparse seed that is not positive
// definite, by the inertia of the factorization, takes tau I added, tau the smallest of a geometric
// sequence that makes it definite, refined by bisection (Nocedal and Wright, Algorithm 3.3).
void sparse_seed_definite(int K, dmod::trust_model::SparseModel::Hess& h, double sgn) {
  dmod::trust_model::SparseModel sm(K);
  std::vector<double> w(K, 1.0), zero(K, 0.0);
  sm.build(h, w.data(), sgn, zero.data());
  if (sm.factor(0.0)) return;
  double lo = 0.0, hi = 1e-8 * sm.scale;
  while (!sm.factor(hi)) {
    lo = hi;
    hi *= 10.0;
    if (!(hi < 1e3 * sm.scale)) return;
  }
  for (int b = 0; b < 8 && lo > 0.0; ++b) {
    const double mid = std::sqrt(lo * hi);
    if (sm.factor(mid)) hi = mid; else lo = mid;
  }
  // tau I joins the sparse part, with a diagonal entry added to a column that has none
  std::vector<int> p(K + 1, 0), i;
  std::vector<double> x;
  for (int j = 0; j < K; ++j) {
    bool has = false;
    for (int k = h.p[j]; k < h.p[j + 1]; ++k) {
      i.push_back(h.i[k]);
      x.push_back(h.x[k] + (h.i[k] == j ? sgn * hi : 0.0));
      has = has || h.i[k] == j;
    }
    if (!has) { i.push_back(j); x.push_back(sgn * hi); }
    p[j + 1] = (int) i.size();
  }
  h.p = std::move(p); h.i = std::move(i); h.x = std::move(x);
  if (!Rf_isNull(h.hessian)) {
    Environment mat = Environment::namespace_env("Matrix");
    Function diagonal = mat["Diagonal"], plus("+");
    h.hessian = plus(h.hessian, diagonal(K, sgn * hi));
  }
}

// Limited-memory BFGS on the sparse model: H = H0 plus the compact correction of Byrd, Nocedal and
// Schnabel (1994) from the pairs (S, Y) in the minimized sign, -W M^-1 W' with W = [B0 S, Y] and
// M = [S'B0 S, L; L', -D], which joins the low-rank part of H0 as W P diag(-1/theta) P' W' from
// M = P diag(theta) P'. False where M is singular to rounding.
bool lbfgs_compose(int K, const dmod::trust_model::SparseModel::Hess& H0,
                   const std::deque<std::vector<double>>& S, const std::deque<std::vector<double>>& Y,
                   double sgn, const CharacterVector& parnames, dmod::trust_model::SparseModel::Hess& H) {
  const int m = (int) S.size(), n = 2 * m;
  std::vector<double> W((std::size_t) K * n);
  for (int l = 0; l < m; ++l) {
    double* BSl = W.data() + (std::size_t) l * K;
    sparse_hess_apply(K, H0, S[l].data(), BSl);
    for (int i = 0; i < K; ++i) BSl[i] *= sgn;
    std::copy(Y[l].begin(), Y[l].end(), W.data() + (std::size_t) (m + l) * K);
  }
  std::vector<double> M((std::size_t) n * n, 0.0);
  auto dot = [&](const double* a, const double* b) {
    double t = 0.0;
    for (int i = 0; i < K; ++i) t += a[i] * b[i];
    return t;
  };
  for (int a = 0; a < m; ++a)
    for (int b = 0; b < m; ++b) {
      M[a + (std::size_t) b * n] = dot(S[a].data(), W.data() + (std::size_t) b * K);
      const double sy = dot(S[a].data(), Y[b].data());
      if (a > b) { M[a + (std::size_t) (m + b) * n] = sy; M[(m + b) + (std::size_t) a * n] = sy; }
      if (a == b) M[(m + a) + (std::size_t) (m + a) * n] = -sy;
    }
  std::vector<double> theta(n), P((std::size_t) n * n);
  eigen_sym_local(M.data(), n, theta.data(), P.data());
  double tmax = 0.0;
  for (int a = 0; a < n; ++a) tmax = std::max(tmax, std::fabs(theta[a]));
  for (int a = 0; a < n; ++a) if (!(std::fabs(theta[a]) > 1e-12 * tmax)) return false;
  const int m0 = (int) H0.d.size();
  H.p = H0.p; H.i = H0.i; H.x = H0.x; H.hessian = H0.hessian;
  H.U = H0.U;
  H.U.resize((std::size_t) K * (m0 + n), 0.0);
  H.d = H0.d;
  H.d.resize(m0 + n);
  for (int e = 0; e < n; ++e) {
    double* Ue = H.U.data() + (std::size_t) (m0 + e) * K;
    for (int a = 0; a < n; ++a) {
      const double pa = P[a + (std::size_t) e * n];
      const double* Wa = W.data() + (std::size_t) a * K;
      for (int i = 0; i < K; ++i) Ue[i] += Wa[i] * pa;
    }
    H.d[m0 + e] = -sgn / theta[e];
  }
  NumericMatrix U(K, m0 + n);
  std::copy(H.U.begin(), H.U.end(), U.begin());
  U.attr("dimnames") = List::create(parnames, R_NilValue);
  H.lowrank = List::create(Named("U") = U, Named("d") = NumericVector(H.d.begin(), H.d.end()));
  return true;
}

// `exact` and `gn` are both pass-throughs of the objective's own Hessian and
// differ only in which one was asked for, so they are told apart by the
// running phase's method rather than by the quasi-Newton state.
// Slot of the source currently producing the model Hessian, ordered as
// kHessianSourceNames, which follows the HessianMethod values.
const char* const kHessianSourceNames[] = {"gn", "bfgs", "sr1", "exact", "sbfgs",
                                           "ssr1", "gnsbfgs", "gnssr1"};
const int kHessianSourceCount = 8;

inline int hessian_source_slot(bool qn_active, int qn_kind, int method) {
  if (!qn_active)        return (method == HM_EXACT) ? HM_EXACT : HM_GN;
  return qn_kind;
}

inline const char* hessian_source_label(bool qn_active, int qn_kind,
                                        int method) {
  return kHessianSourceNames[hessian_source_slot(qn_active, qn_kind, method)];
}

// BFGS update A + z z' / z's - u u' / s'u, u = A s, of a symmetric `A` in place, so that A s = z
// afterwards also where s'As <= 0, and from A = 0. A pair without positive curvature z's, or with
// s'As = 0 but As != 0, leaves `A` unchanged. Returns whether `A` changed.
bool bfgs_update_secant(int K, std::vector<double>& A, const std::vector<double>& s,
                        const std::vector<double>& z) {
  std::vector<double> u(K, 0.0);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) u[i] += A[i + (std::size_t) j * K] * s[j];
  double c = 0.0, b = 0.0, un = 0.0;
  for (int i = 0; i < K; ++i) { c += u[i] * s[i]; b += z[i] * s[i]; un += u[i] * u[i]; }
  if (!(b > 0.0)) return false;
  if (c == 0.0 && un > 0.0) return false;
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) {
      double d = z[i] * z[j] / b;
      if (c != 0.0) d -= u[i] * u[j] / c;
      A[i + (std::size_t) j * K] += d;
    }
  return true;
}

// The data rows an objective returns with `structure` (see objlist()), in the sign of the
// objective: the Jacobian J of the predictions and sigmas, column-major with one row per data row,
// the derivatives u of the value in them, and the Gauss-Newton Hessian G_LS of the data term.
// n < 0 where the evaluation returned none.
struct Rows {
  std::vector<double> J, u, GLS;
  int n = -1;
  bool ok() const { return n >= 0; }
};

bool read_rows(const List& out, int K, Rows& r) {
  r.n = -1;
  if (!out.containsElementNamed("structure")) return false;
  SEXP st = out["structure"];
  if (TYPEOF(st) != VECSXP) return false;
  List L(st);
  if (!L.containsElementNamed("jacobian") || !L.containsElementNamed("seeds") ||
      !L.containsElementNamed("hessian")) return false;
  SEXP J = L["jacobian"], u = L["seeds"], H = L["hessian"];
  if (TYPEOF(J) != REALSXP || TYPEOF(u) != REALSXP || TYPEOF(H) != REALSXP) return false;
  const int n = Rf_length(u);
  if ((std::size_t) Rf_length(J) != (std::size_t) n * K || Rf_length(H) != K * K) return false;
  r.J.assign(REAL(J), REAL(J) + (std::size_t) n * K);
  r.u.assign(REAL(u), REAL(u) + n);
  r.GLS.assign(REAL(H), REAL(H) + (std::size_t) K * K);
  r.n = n;
  return true;
}

// The residual norm of the scaled methods: |u| with data rows, sqrt(f) without, which agree up to
// a constant factor for a sum of squares.
double residual_norm(const Rows& r, double f) {
  if (!r.ok()) return std::sqrt(std::max(f, 0.0));
  double uu = 0.0;
  for (double v : r.u) uu += v * v;
  return std::sqrt(uu);
}

// State of a structured phase, in the minimized sign: the Gauss-Newton Hessian G and the
// residual norm rho at the current iterate, the correction A and, for a hybrid method, the
// branch of the last accepted step. With residual scaling, A approximates S / rho.
struct Structured {
  std::vector<double> G, A;
  double rho = 0.0;
  bool gn_branch = false;

  void reset(int K, const std::vector<double>& H, double sgn, double rho0) {
    const std::size_t KK = (std::size_t) K * K;
    G.resize(KK); A.assign(KK, 0.0);
    for (std::size_t t = 0; t < KK; ++t) G[t] = sgn * H[t];
    rho = rho0;
    gn_branch = false;
  }
  // The model Hessian in the sign of the objective: G + c A, c = rho with scaling and 1
  // otherwise, or on the Gauss-Newton branch of a hybrid method G, or G + rho I with `shift`.
  void compose(int K, std::vector<double>& H, double sgn, int method, bool scaled,
               bool shift) const {
    const std::size_t KK = (std::size_t) K * K;
    H.resize(KK);
    const bool gn = hm_hybrid(method) && gn_branch;
    const double c = gn ? 0.0 : (scaled ? rho : 1.0);
    for (std::size_t t = 0; t < KK; ++t) H[t] = sgn * (G[t] + c * A[t]);
    if (gn && shift)
      for (int i = 0; i < K; ++i) H[i + (std::size_t) i * K] += sgn * rho;
  }
};

// Update of the structured correction from the pair (s, y), y = g+ - g, with G+ and rho+ at the
// trial point, all in the minimized sign, and the data rows `rc` at the iterate and `rt` at the
// trial point. The structured secant is y# = y - J'(u+ - u) - (G+ - G_LS+) s, which is (J+ - J)'u+
// for the data term and y - G+ s for the other terms; without data rows y - G+ s throughout. The
// scaled methods use z = y# rho+ / rho. Returns 1 for an update, 0 for a skipped pair, and for a
// hybrid method 2 where the switch rule takes the Gauss-Newton branch.
int structured_update(int method, const QnParams& q, int K, Structured& st,
                      const std::vector<double>& s, const std::vector<double>& y,
                      const std::vector<double>& Gp, double rho_p,
                      const Rows& rc, const Rows& rt, double sgn) {
  const std::size_t KK = (std::size_t) K * K;
  std::vector<double> ysh(y);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) ysh[i] -= Gp[i + (std::size_t) j * K] * s[j];
  if (rc.ok() && rt.ok() && rc.n == rt.n) {
    const int n = rc.n;
    std::vector<double> du(n);
    for (int r = 0; r < n; ++r) du[r] = rt.u[r] - rc.u[r];
    for (int k = 0; k < K; ++k) {
      const double* Jk = rc.J.data() + (std::size_t) k * n;
      double v = 0.0;
      for (int r = 0; r < n; ++r) v += Jk[r] * du[r];
      ysh[k] -= sgn * v;
    }
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) ysh[i] += sgn * rt.GLS[i + (std::size_t) j * K] * s[j];
  }
  const bool hybrid = hm_hybrid(method);
  const bool scaled = hybrid || q.scaled;
  std::vector<double> z(ysh);
  if (scaled) {
    if (!(st.rho > 0.0)) return 0;
    for (int i = 0; i < K; ++i) z[i] *= rho_p / st.rho;
  }

  if (hybrid) {
    double zs = 0.0, ss = 0.0;
    for (int i = 0; i < K; ++i) { zs += z[i] * s[i]; ss += s[i] * s[i]; }
    const double ratio = (method == HM_GNSSR1) ? std::fabs(zs) / ss : zs / ss;
    if (!(ratio >= q.hybrid_tol)) return 2;
    if (method == HM_GNSBFGS) return bfgs_update_secant(K, st.A, s, z) ? 1 : 0;
    return qn_update_B(HM_SR1, K, st.A, s, z) ? 1 : 0;
  }

  // Huschens: the update acts on rho+ A, the estimate of S at the trial point.
  if (scaled && !(rho_p > 0.0)) return 0;
  const double c = scaled ? rho_p : 1.0;
  std::vector<double> At(KK);
  for (std::size_t t = 0; t < KK; ++t) At[t] = c * st.A[t];
  bool changed;
  if (method == HM_SSR1) {
    changed = qn_update_B(HM_SR1, K, At, s, z);
  } else {
    // BFGS on G+ + A with the secant G+ s + z, undamped: G+ changes from iterate to iterate,
    // so G+ + A need not be positive definite, and the update restores the secant along s.
    std::vector<double> Bs(KK), ys(z);
    for (std::size_t t = 0; t < KK; ++t) Bs[t] = Gp[t] + At[t];
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) ys[i] += Gp[i + (std::size_t) j * K] * s[j];
    if (!qn_pair_informative(K, s, ys, q.cautious)) return 0;
    changed = bfgs_update_secant(K, Bs, s, ys);
    if (changed) for (std::size_t t = 0; t < KK; ++t) At[t] = Bs[t] - Gp[t];
  }
  if (!changed) return 0;
  for (std::size_t t = 0; t < KK; ++t) st.A[t] = At[t] / c;
  return 1;
}

// ---- Coleman-Li interior trust-region-reflective --------------------------
// One kernel for both curvature models of trust_model.h: the dense one and the sparse one with
// a low-rank part take the same steps; quasi-Newton phases update a dense Hessian.
template <class Model>
List trust_reflective(Function objfun, NumericVector parinit,
                      double rinit, double rmax,
                      Nullable<NumericVector> parscale,
                      int iterlim,
                      double ftol, double mtol,
                      double gtol, double xtol,
                      double rmin, double thetamax,
                      const std::vector<Phase>& phases, int switches,
                      double nonmonotone,
                      bool minimize, bool blather_on,
                      Nullable<NumericVector> parupper,
                      Nullable<NumericVector> parlower,
                      bool printIter, Nullable<CharacterVector> traceFile) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  for (int i = 0; i < K; ++i)
    if (!(pl[i] < pu[i])) stop("trust: parlower must lie strictly below parupper");
  fill_parscale(parscale, ps, K, "trust");

  // Work in the parscale frame z = parscale * theta, where the radius lives.
  std::vector<double> lbz(K), ubz(K), z(K);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    lbz[i] = ps[i] * pl[i];
    ubz[i] = ps[i] * pu[i];
    double zi = ps[i] * parinit[i];
    if (zi > ubz[i]) { any_above = true; zi = ubz[i]; }
    if (zi < lbz[i]) { any_below = true; zi = lbz[i]; }
    z[i] = zi;
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");
  push_interior(K, z, lbz, ubz);

  Reporter report;
  report.init(printIter, iterlim, traceFile, K, parnames);

  NumericVector x_named(K);
  x_named.names() = parnames;
  for (int i = 0; i < K; ++i) x_named[i] = z[i] / ps[i];

  // Only a run that starts quasi-Newton can do without the objective's Hessian;
  // gn needs it every iteration and a later phase seeds the handover from it.
  const int hessianMethod = phases[0].method;
  const QnParams& qpri = phases[0].qn;
  const bool qn_start = (hessianMethod == HM_BFGS || hessianMethod == HM_SR1);
  const bool seed_gn  = !(qn_start && qpri.init == HI_IDENTITY);
  // An exact seed is the whole argument for the reverse mode: one expensive
  // evaluation at the start buys a true curvature, and the descent then runs on
  // gradients alone. `hessianMethod = "exact"` asks for it at every iterate.
  const bool seed_exact =
      (qpri.init == HI_EXACT) || (hessianMethod == HM_EXACT);

  List out_init;
  if (!eval_objfun(objfun, x_named, out_init,
                   seed_exact ? CV_EXACT : (seed_gn ? CV_GN : CV_GRADIENT), 0))
    stop_parinit_failed();
  double val = as<double>(out_init["value"]);
  if (!std::isfinite(val)) stop("parinit not feasible: value is not finite");
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);

  std::vector<double> grad_full(grad0.begin(), grad0.end());
  typename Model::Hess H_full;
  Model model(K);
  bool seeded = false;
  if (seed_gn) {
    seeded = Model::read(out_init, K, H_full);
    // gn reads a Hessian at every iterate, so an objective that declines one
    // cannot drive it at all; a quasi-Newton start only loses the seed.
    if (!seeded && !qn_start)
      stop("objfun returned no Hessian at parinit, which this hessianMethod "
           "needs at every iterate");
    if (!seeded)
      Rf_warning("objfun returned no Hessian at parinit; seeding the "
                 "quasi-Newton approximation with the identity instead");
  }
  if (!seeded) Model::identity(K, minimize ? 1.0 : -1.0, H_full);
  if constexpr (Model::dense)
    if (seeded && hessianMethod == HM_BFGS) bfgs_seed_definite(K, H_full.H, minimize ? 1.0 : -1.0);
  Structured st;
  int n_hybrid_gn = 0;
  // The data rows at the current iterate and at the trial point; `rows_seen` once the objective
  // returned any.
  Rows rows_cur, rows_try;
  bool rows_seen = read_rows(out_init, K, rows_cur);
  if constexpr (Model::dense)
    if (hm_structured(hessianMethod))
      st.reset(K, H_full.H, minimize ? 1.0 : -1.0, residual_norm(rows_cur, minimize ? val : -val));

  // The sparse model of a BFGS phase: the Hessian of the objective at the start of the phase, made
  // positive definite, to which the kept pairs add the compact correction.
  typename Model::Hess H_lb0;
  if constexpr (!Model::dense) if (hessianMethod == HM_BFGS) {
    sparse_seed_definite(K, H_full, minimize ? 1.0 : -1.0);
    H_lb0 = H_full;
  }

  int neval = 1;
  report(neval, val, x_named, /*head=*/true);

  Blather trace;
  std::vector<double> absv(K), jv(K), sqrtv(K), g_z(K), ghat(K), wscale(K), cdiag(K);
  std::vector<double> shat(K), shat_step(K), s_step(K), z_try(K);
  std::vector<unsigned char> at_bound(K, 0);

  double r = rinit;
  double f_used = minimize ? val : -val;
  double opt_measure = kInf;
  bool   accept = true, converged = false, bail = false;
  int    n_iter = 0, n_fail = 0, n_stall = 0;
  std::string stop_reason = "iterlim";
  // The convergence stop held back by newton_small() since the last accepted step that changed the
  // value by ftol or more. A run that then ends on the radius ends on that stop.
  const char* vetoed = nullptr;
  auto end_on_radius = [&]() {
    if (vetoed != nullptr) { converged = true; stop_reason = vetoed; }
    else stop_reason = "radius";
  };

  // Zhang-Hager reference value for acceptance. nonmonotone = 0 holds zh_C at
  // f_used, which is the monotone rule.
  double zh_C = f_used, zh_Q = 1.0;
  // Acceptance may take an increase, so the last iterate need not be the best
  // one visited. Snapshots are taken where all of these agree on one iterate.
  double best_f = kInf, best_val = val;
  std::vector<double> best_z = z, best_grad = grad_full;
  typename Model::Hess best_H = H_full;
  std::vector<unsigned char> best_at_bound = at_bound;

  // Hessian source: gn passes the objective's J^T J through, bfgs and sr1
  // maintain their own update seeded from it (the seed already filled H_full).
  // At a soft stop the next phase of the list takes over, cyclically, while
  // switches remain.
  const int n_phases = (int) phases.size();
  int  cur = 0, method = hessianMethod;
  bool qn_dirty = false;
  // Settings of the phase running now.
  const QnParams* qp = &qpri;
  // What an iterate of the running phase that wants curvature asks for.
  int  want_level = (method == HM_EXACT) ? CV_EXACT : CV_GN;
  // An evaluation at the current iterate before the next step: the Hessian of a
  // handover or reseed at `pending_level`, and with `pending_refresh` value and
  // gradient as the new phase computes them.
  int  pending_level = -1;
  bool pending_refresh = false;
  // Iterations of the running phase, against the phase budget, and of every
  // phase run so far in run order. The iteration of a handover counts on the
  // phase it ends.
  int  phase_iter = 0;
  std::vector<int> phase_iters(1, 0);
  // Consecutive iterations of the running phase that left the radius unchanged.
  int  radius_same = 0;
  // "exact" is a Newton run: the objective's own Hessian at every iterate, the
  // way gn passes J^T J through, and the subproblem solver takes an indefinite
  // matrix natively.
  bool qn_active = (method != HM_GN && method != HM_EXACT);
  int  qn_kind  = hm_structured(method) ? method : ((method == HM_SR1) ? HM_SR1 : HM_BFGS);
  int  qn_neval = 0, qn_skipped = 0, n_switch = 0;
  // At a standing iterate a reseed refetches the same matrix. The allowance
  // kReseedLimit holds per phase, n_reseed counts over the run.
  int  n_reseed = 0, phase_reseed = 0;
  bool moved_since_reseed = true;
  // The initial evaluation belongs to whichever source the run starts on.
  std::vector<int> eval_src(kHessianSourceCount, 0);
  eval_src[hessian_source_slot(qn_active, qn_kind, method)] = 1;
  const double qn_sgn = minimize ? 1.0 : -1.0;
  std::vector<double> s_qn(K), y_qn(K);
  // The BFGS matrix of a phase seeded "tracked": the start Hessian made positive definite, updated
  // with every accepted pair while no tracked phase runs, and taken over by such a phase.
  int track_k = -1;
  for (int k = 0; k < (int) phases.size() && track_k < 0; ++k)
    if (phases[k].qn.init == HI_TRACKED) track_k = k;
  std::vector<double> H_track;
  if constexpr (Model::dense) if (track_k >= 0) {
    H_track = H_full.H;
    bfgs_seed_definite(K, H_track, qn_sgn);
  }
  // qnMemory = 0 accumulates onto the seed, otherwise rebuild from the window.
  std::deque<std::vector<double>> qn_S, qn_Y;
  double qn_gamma = 1.0;
  // A stalled quasi-Newton phase can fetch a fresh curvature instead of stopping. A reseed, not a
  // handover: the source stays and only the matrix it updates from is replaced, so the stored
  // pairs go with it. Under "stall" it fetches the kind the phase was seeded with.
  auto can_reseed = [&]() {
    return qp->reseed == HR_STALL && qn_active && moved_since_reseed &&
           phase_reseed < kReseedLimit;
  };
  auto reseed = [&]() {
    pending_level = (qp->init == HI_EXACT) ? CV_EXACT : CV_GN;
    n_reseed++; phase_reseed++; moved_since_reseed = false;
    n_stall = 0; r = rinit; radius_same = 0;
    qn_S.clear(); qn_Y.clear(); qn_gamma = 1.0;
  };

  // Scaling, optimality measure, bound flags and the best-iterate snapshot at
  // the current point. Runs at the head of an iteration and once more after the
  // loop, which is where the soft stops and the iteration limit leave it.
  auto refresh_current = [&]() {
    f_used = minimize ? val : -val;
    const double sgn = minimize ? 1.0 : -1.0;
    for (int i = 0; i < K; ++i) g_z[i] = sgn * grad_full[i] / ps[i];

    affine_scaling(K, z.data(), g_z.data(), lbz.data(), ubz.data(),
                   absv.data(), jv.data());
    for (int i = 0; i < K; ++i) sqrtv[i] = std::sqrt(absv[i]);

    // First-order optimality for the box problem: |v| * |g| vanishes both at
    // an interior stationary point and at a bound the gradient pushes into.
    opt_measure = 0.0;
    for (int i = 0; i < K; ++i)
      opt_measure = std::max(opt_measure, std::fabs(absv[i] * g_z[i]));
    // A coordinate counts as bound-active when only the scaling makes it
    // stationary, the gradient entry of the coordinate still well away from zero. Set before
    // the convergence break, which is exactly when it matters.
    const double btol = std::max(gtol, 1e-10);
    for (int i = 0; i < K; ++i)
      at_bound[i] = (jv[i] > 0.0 && std::fabs(absv[i] * g_z[i]) <= btol &&
                     std::fabs(g_z[i]) > btol) ? 1 : 0;

    if (f_used <= best_f) {
      best_f = f_used; best_val = val; best_z = z; best_grad = grad_full;
      best_H = H_full; best_at_bound = at_bound;
    }
  };

  // The scaled model (7) of the methods article at the current iterate.
  auto build_model = [&]() {
    const double sgn = minimize ? 1.0 : -1.0;
    for (int i = 0; i < K; ++i) {
      ghat[i] = sqrtv[i] * g_z[i];
      wscale[i] = sqrtv[i] / ps[i];
      cdiag[i] = std::fabs(g_z[i]) * jv[i];
    }
    model.build(H_full, wscale.data(), sgn, cdiag.data());
  };
  // A convergence stop also needs the full Newton step of the built model to predict at most
  // mtol; mtol = 0 drops the condition.
  auto newton_small = [&]() {
    return mtol <= 0.0 || model.newton_decrease(ghat.data()) <= mtol;
  };

  for (int iter = 1; iter <= iterlim; ++iter) {
    R_CheckUserInterrupt();

    // The Hessian a handover or reseed needs at the current iterate, and after
    // a change of sweep value and gradient as the new phase computes them. A
    // failed evaluation or a declined Hessian keeps the approximation.
    if (pending_level >= 0) {
      const int level = pending_level;
      const bool refresh = pending_refresh;
      pending_level = -1; pending_refresh = false;
      NumericVector x_now(K);
      x_now.names() = parnames;
      for (int i = 0; i < K; ++i) x_now[i] = z[i] / ps[i];
      // A Gauss-Newton Hessian is computed forward, so a reverse phase takes
      // the value and gradient from a second evaluation.
      const bool split = refresh && phases[cur].reverse && level == CV_GN;
      auto take = [&](const List& out) {
        if (!out.containsElementNamed("value") ||
            !out.containsElementNamed("gradient")) return;
        const double v = as<double>(out["value"]);
        if (!std::isfinite(v)) return;
        NumericVector g = as<NumericVector>(out["gradient"]);
        if (g.size() != K) return;
        val = v;
        grad_full.assign(g.begin(), g.end());
        zh_C = minimize ? val : -val; zh_Q = 1.0;
      };
      List out_now;
      if (eval_objfun(objfun, x_now, out_now, level, cur)) {
        if (level >= CV_GN && Model::read(out_now, K, H_full)) {
          if constexpr (Model::dense)
            if (qn_active && qn_kind == HM_BFGS) bfgs_seed_definite(K, H_full.H, qn_sgn);
          if constexpr (!Model::dense)
            if (qn_active && qn_kind == HM_BFGS) {
              sparse_seed_definite(K, H_full, qn_sgn);
              H_lb0 = H_full; qn_S.clear(); qn_Y.clear();
            }
        }
        if (level >= CV_GN) {
          Rows fresh;
          if (read_rows(out_now, K, fresh)) { rows_cur = std::move(fresh); rows_seen = true; }
        }
        if (refresh && !split) take(out_now);
        // A structured phase starts from A = 0 at the fresh Gauss-Newton Hessian.
        if constexpr (Model::dense)
          if (qn_active && hm_structured(qn_kind))
            st.reset(K, H_full.H, qn_sgn, residual_norm(rows_cur, qn_sgn * val));
        neval++;
        eval_src[hessian_source_slot(qn_active, qn_kind, method)]++;
      }
      if (split) {
        List out_g;
        if (eval_objfun(objfun, x_now, out_g, CV_GRADIENT, cur)) {
          take(out_g);
          neval++;
          eval_src[hessian_source_slot(qn_active, qn_kind, method)]++;
        }
      }
      qn_dirty = true;
    }

    // A rejected-step SR1 update changes H_full without moving the iterate, so
    // qn_dirty rebuilds the model too.
    if (accept || qn_dirty) {
      qn_dirty = false;
      refresh_current();
      build_model();
      if (opt_measure <= gtol) {
        if (newton_small()) { converged = true; stop_reason = "gradient"; break; }
        vetoed = "gradient";
      }
    }

    if (blather_on) {
      for (int i = 0; i < K; ++i) trace.argpath.push_back(z[i] / ps[i]);
      trace.r.push_back(r);
      trace.valpath.push_back(val);
    }

    const bool hyb_gn_now = st.gn_branch;
    const double r_head = r;
    bool is_newton = false, is_hard = false, is_easy = false;
    model.solve(ghat.data(), r, shat.data(), &is_newton, &is_hard, &is_easy);

    // Let the iterate approach a face as the optimality measure falls, but
    // never reach it: |v| = 0 would freeze that coordinate for good.
    const double theta_frac =
        std::min(std::max(thetamax, 1.0 - opt_measure), 1.0 - 1e-12);

    const char* sb_label = "full";
    double m_value = 0.0;
    stepback(K, z.data(), lbz.data(), ubz.data(), sqrtv.data(),
             ghat.data(), model, shat.data(), r, theta_frac,
             shat_step.data(), s_step.data(), &m_value, &sb_label);

    double stepnorm = 0.0;
    for (int i = 0; i < K; ++i) stepnorm += shat_step[i] * shat_step[i];
    stepnorm = std::sqrt(stepnorm);

    for (int i = 0; i < K; ++i) z_try[i] = z[i] + s_step[i];
    // theta_frac keeps the step interior in exact arithmetic, but once |v| is
    // small the remaining gap underflows and z_try rounds onto the bound,
    // setting |v| = 0 and freezing the coordinate.
    push_interior(K, z_try, lbz, ubz);

    // A step below the resolution of the iterate would evaluate the iterate itself, and every
    // later step is smaller still: the radius has collapsed whatever rmin is. A quasi-Newton
    // phase that may still reseed stalls instead.
    bool moved = false;
    for (int i = 0; i < K && !moved; ++i) moved = z_try[i] != z[i];
    if (!moved && can_reseed()) { reseed(); continue; }
    if (!moved) { end_on_radius(); break; }

    NumericVector x_try(K);
    x_try.names() = parnames;
    for (int i = 0; i < K; ++i) x_try[i] = z_try[i] / ps[i];

    // The quasi-Newton phase maintains a separate Hessian, so it asks the
    // objective not to build J^T J; only the gn pass-through reads it back.
    const bool structured = qn_active && hm_structured(qn_kind);
    const bool want_h = !qn_active || structured;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try,
                               want_h ? want_level : CV_GRADIENT, cur);
    double val_try = kInf;
    NumericVector grad_try;
    typename Model::Hess Htry;
    if (eval_ok) {
      val_try  = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      // A gn iteration without a Hessian has no model to build, so the point
      // counts as a failed evaluation rather than as a silent fallback.
      if (want_h && !Model::read(out_try, K, Htry)) eval_ok = false;
      if (!std::isfinite(val_try)) eval_ok = false;
    }
    rows_try.n = -1;
    if (eval_ok && want_h && read_rows(out_try, K, rows_try)) rows_seen = true;
    neval++;
    if (qn_active) qn_neval++;
    eval_src[hessian_source_slot(qn_active, qn_kind, method)]++;
    report(neval, val_try, x_try, /*head=*/false);

    const double pred_pos  = -m_value;
    const double ftry_used = minimize ? val_try : -val_try;
    const double dval      = std::fabs(ftry_used - f_used);
    const double rho = (eval_ok && pred_pos > 0.0)
                         ? (zh_C - ftry_used) / pred_pos : -kInf;

    if (!eval_ok) {
      n_fail++;
      accept = false;
      r *= 0.25;
      if (n_fail >= 3) { bail = true; stop_reason = "objfun"; }
    } else {
      n_fail = 0;
      if (rho < 0.25) {
        accept = false;
        r = std::min(0.25 * r, 0.25 * stepnorm);
      } else {
        accept = true;
        // Only grow when the solution actually pressed against the radius.
        if (rho > 0.75 && stepnorm >= 0.9 * r) r = std::min(2.0 * r, rmax);
      }
    }

    // Count rejected steps that also left the objective flat.
    if (accept || !eval_ok || dval >= ftol) n_stall = 0;
    else                                                              n_stall++;
    if (accept && eval_ok && dval >= ftol) vetoed = nullptr;

    // SR1 also takes the pair from a rejected trial point, where the curvature
    // along s is still valid (Nocedal and Wright, sec. 6.2). BFGS needs an
    // accepted step for the curvature condition.
    const bool qn_pair = qn_active && eval_ok &&
                         (accept || (qp->rejected && hm_sr1_type(qn_kind)));
    const double rho_try = residual_norm(rows_try, ftry_used);

    // (s, y) in minimized sign, before z/grad_full are overwritten.
    if (qn_pair)
      for (int i = 0; i < K; ++i) {
        s_qn[i] = (z_try[i] - z[i]) / ps[i];
        y_qn[i] = qn_sgn * (grad_try[i] - grad_full[i]);
      }

    if constexpr (Model::dense)
      if (track_k >= 0 && accept && eval_ok && !(qn_active && qp->init == HI_TRACKED)) {
        std::vector<double> s_tr(K), y_tr(K);
        for (int i = 0; i < K; ++i) {
          s_tr[i] = (z_try[i] - z[i]) / ps[i];
          y_tr[i] = qn_sgn * (grad_try[i] - grad_full[i]);
        }
        if (qn_pair_informative(K, s_tr, y_tr, phases[track_k].qn.cautious))
          qn_update(HM_BFGS, K, H_track, s_tr, y_tr, qn_sgn);
        else qn_skipped++;
      }

    if (accept && eval_ok) {
      moved_since_reseed = true;
      z = z_try;
      val = val_try;
      grad_full.assign(grad_try.begin(), grad_try.end());
      const double q_new = nonmonotone * zh_Q + 1.0;
      zh_C = (nonmonotone * zh_Q * zh_C + ftry_used) / q_new;
      zh_Q = q_new;
      if (!qn_active) H_full = std::move(Htry);
    }

    // The structured correction is updated with G+ of the trial point; an accepted step also
    // moves G and rho, and the model is composed afresh.
    if constexpr (Model::dense) if (structured && eval_ok && (qn_pair || accept)) {
      std::vector<double> Gp(Htry.H.size());
      for (std::size_t t = 0; t < Gp.size(); ++t) Gp[t] = qn_sgn * Htry.H[t];
      if (qn_pair) {
        const int res = structured_update(qn_kind, *qp, K, st, s_qn, y_qn, Gp, rho_try,
                                          rows_cur, rows_try, qn_sgn);
        if (res == 0) qn_skipped++;
        if (accept && hm_hybrid(qn_kind)) {
          st.gn_branch = (res == 2);
          if (st.gn_branch) n_hybrid_gn++;
        }
      }
      if (accept) { st.G = std::move(Gp); st.rho = rho_try; }
      st.compose(K, H_full.H, qn_sgn, qn_kind, qp->scaled, qp->hybrid_shift);
      if (!accept) qn_dirty = true;
    }

    // The data rows move with the iterate; a trial point without them leaves none.
    if (accept && eval_ok) { rows_cur = std::move(rows_try); rows_try = Rows(); }

    // The pair is damped as in the dense update, against the model Hessian it updates, and kept so.
    if constexpr (!Model::dense) if (qn_pair && !structured) {
      std::vector<double> Bs(K);
      sparse_hess_apply(K, H_full, s_qn.data(), Bs.data());
      double sy = 0.0, sBs = 0.0;
      for (int i = 0; i < K; ++i) { Bs[i] *= qn_sgn; sy += s_qn[i] * y_qn[i]; sBs += s_qn[i] * Bs[i]; }
      const double th = (sy < 0.2 * sBs) ? 0.8 * sBs / (sBs - sy) : 1.0;
      std::vector<double> yd(K);
      double syd = 0.0;
      for (int i = 0; i < K; ++i) { yd[i] = th * y_qn[i] + (1.0 - th) * Bs[i]; syd += s_qn[i] * yd[i]; }
      if (!qn_pair_informative(K, s_qn, y_qn, qp->cautious) || !(sBs > 0.0) || !(syd > 0.0)) {
        qn_skipped++;
      } else {
        const int mem = qp->memory > 0 ? qp->memory : kSparseMemory;
        qn_S.push_back(s_qn);
        qn_Y.push_back(yd);
        if ((int) qn_S.size() > mem) { qn_S.pop_front(); qn_Y.pop_front(); }
        // a singular middle matrix drops the oldest pairs
        while (!lbfgs_compose(K, H_lb0, qn_S, qn_Y, qn_sgn, parnames, H_full)) {
          qn_S.pop_front(); qn_Y.pop_front();
          if (qn_S.empty()) { H_full = H_lb0; break; }
        }
        if (!accept) qn_dirty = true;
      }
    }

    if constexpr (Model::dense) if (qn_pair && !structured) {
      const bool keep = (qn_kind == HM_SR1) ||
                        qn_pair_informative(K, s_qn, y_qn, qp->cautious);
      if (!keep) {
        qn_skipped++;
      } else {
        double sy = 0.0, yy = 0.0;
        for (int i = 0; i < K; ++i) { sy += s_qn[i] * y_qn[i]; yy += y_qn[i] * y_qn[i]; }
        if (sy > 0.0 && yy > 0.0) qn_gamma = yy / sy;
        if (qp->memory > 0) {
          qn_S.push_back(s_qn);
          qn_Y.push_back(y_qn);
          if ((int) qn_S.size() > qp->memory) { qn_S.pop_front(); qn_Y.pop_front(); }
          qn_assemble(qn_kind, K, H_full.H, qn_gamma, qn_S, qn_Y, qn_sgn);
        } else {
          qn_update(qn_kind, K, H_full.H, s_qn, y_qn, qn_sgn);
        }
        if (!accept) qn_dirty = true;
      }
    }

    if (blather_on) {
      for (int i = 0; i < K; ++i) trace.argtry.push_back(z_try[i] / ps[i]);
      trace.valtry.push_back(val_try);
      trace.accept.push_back(accept ? 1 : 0);
      trace.preddiff.push_back(m_value);
      trace.stepnorm.push_back(stepnorm);
      trace.rho.push_back(rho);
      trace.steptype.push_back(subproblem_label(is_newton, is_hard, is_easy));
      trace.stepback.push_back(sb_label);
      // a hybrid phase on the Gauss-Newton branch reports "gn"
      trace.hsource.push_back(qn_active && hm_hybrid(qn_kind) && hyb_gn_now ? "gn" :
          hessian_source_label(qn_active, qn_kind, method));
    }
    n_iter = iter;

    if (bail) break;

    // Soft stops read the last accepted step rather than first-order
    // optimality, so a misleading model can trip them far from a stationary
    // point. The gradient test above is not one and stays terminal.
    const char* soft_stop = nullptr;
    bool model_stop = false;
    if (accept && eval_ok) {
      if (dval < ftol)                        { soft_stop = "fvalue";   model_stop = true; }
      else if (std::fabs(m_value) < mtol)     { soft_stop = "preddiff"; model_stop = true; }
      else if (xtol > 0.0 && stepnorm < xtol) { soft_stop = "step";     model_stop = true; }
    }
    // A rejected Newton step inside the radius that predicts less than mtol: no smaller step
    // predicts more, the model is at the minimum.
    if (!accept && eval_ok && is_newton && std::fabs(m_value) < mtol) {
      soft_stop = "preddiff"; model_stop = true;
    }
    // "fvalue" and "preddiff" end the run only where the full Newton step predicts at most mtol;
    // a quasi-Newton phase drops them further down.
    const bool ends_run = model_stop && (n_switch >= switches || phases[cur].stop_ends) && !qn_active;
    if (ends_run && std::strcmp(soft_stop, "step") != 0) {
      refresh_current();
      build_model();
      if (!newton_small()) { vetoed = soft_stop; soft_stop = nullptr; model_stop = false; }
    }
    if (soft_stop == nullptr && r >= rmin && n_stall >= kStallLimit)
      soft_stop = "stagnation";

    if (soft_stop != nullptr && std::strcmp(soft_stop, "stagnation") == 0 && can_reseed()) {
      reseed();
      soft_stop = nullptr;
    }
    // A phase that used up the phase budget ends like a soft stop, but without the
    // claim of convergence.
    bool budget_out = false;
    phase_iters.back()++;
    if (++phase_iter >= phases[cur].iterlim && soft_stop == nullptr) {
      soft_stop = "iterlim";
      budget_out = true;
    }
    // A radius unchanged for radius_stall iterations hands over; with no handover left it is no stop.
    radius_same = (r == r_head) ? radius_same + 1 : 0;
    if (soft_stop == nullptr && phases[cur].radius_stall > 0 &&
        radius_same >= phases[cur].radius_stall && n_switch < switches)
      soft_stop = "radiusStall";

    // A soft stop hands over to the next phase and restarts the trust region, unless the phase
    // ends the run on a stop. A quasi-Newton phase seeded "gn" starts from the last accepted
    // J^T J, which is why H_full is not reset for it.
    const bool by_stop = soft_stop != nullptr && std::strcmp(soft_stop, "iterlim") != 0 &&
                         std::strcmp(soft_stop, "radiusStall") != 0;
    if (soft_stop != nullptr && n_switch < switches && !(by_stop && phases[cur].stop_ends)) {
      if constexpr (Model::dense) if (qn_active && qp->init == HI_TRACKED) H_track = H_full.H;
      const bool was_qn = qn_active;
      const int prev_method = method;
      const int next = (cur + 1) % n_phases;
      const bool refreshOnSwitch = phases[next].reverse != phases[cur].reverse;
      cur = next;
      qp = &phases[cur].qn;
      method = phases[cur].method;
      want_level = (method == HM_EXACT) ? CV_EXACT : CV_GN;
      qn_active = (method != HM_GN && method != HM_EXACT);
      if (qn_active && hm_structured(method)) {
        // From a gn phase, H_full is the Gauss-Newton Hessian at the iterate; otherwise it is
        // fetched, with value and gradient where the sweep changes.
        qn_kind = method;
        if constexpr (Model::dense) {
          st.reset(K, H_full.H, qn_sgn, residual_norm(rows_cur, qn_sgn * val));
          if (prev_method != HM_GN || (rows_seen && !rows_cur.ok())) pending_level = CV_GN;
          else if (refreshOnSwitch)  pending_level = CV_GRADIENT;
        }
        qn_dirty = true;
      } else if (qn_active) {
        qn_kind = (method == HM_SR1) ? HM_SR1 : HM_BFGS;
        if (qp->init == HI_TRACKED) {
          if constexpr (Model::dense) H_full.H = H_track;
          qn_dirty = true;
        } else if (qp->init == HI_IDENTITY) {
          Model::identity(K, qn_sgn, H_full);
          qn_dirty = true;
        } else if constexpr (Model::dense) {
          if (qn_kind == HM_BFGS) { bfgs_seed_definite(K, H_full.H, qn_sgn); qn_dirty = true; }
        } else {
          sparse_seed_definite(K, H_full, qn_sgn);
          H_lb0 = H_full;
          qn_dirty = true;
        }
        if (qp->init == HI_EXACT)                 pending_level = CV_EXACT;
        else if (qp->init == HI_GN && was_qn)     pending_level = CV_GN;
        else if (refreshOnSwitch)                 pending_level = CV_GRADIENT;
      } else {
        pending_level = want_level;
      }
      pending_refresh = refreshOnSwitch;
      phase_iter = 0; phase_reseed = 0; radius_same = 0;
      phase_iters.push_back(0);
      n_switch++; n_stall = 0;
      if (!phases[cur].keep_radius) r = rinit;
      qn_S.clear(); qn_Y.clear(); qn_gamma = 1.0;
      soft_stop = nullptr;
      budget_out = false;
    }
    // The three tests above read the quadratic model, which in a quasi-Newton
    // phase describes the approximation rather than the iterate. Such a run
    // ends on the gradient; stagnation reads accepted steps and still applies.
    if (soft_stop != nullptr && qn_active && model_stop) soft_stop = nullptr;
    if (soft_stop != nullptr) {
      converged = !budget_out; stop_reason = soft_stop;
      // Stagnation is convergence only where the model predicts no decrease above mtol, or a
      // convergence stop was held back since the last step that changed the value by ftol.
      if (converged && std::strcmp(soft_stop, "stagnation") == 0 && vetoed == nullptr) {
        refresh_current();
        build_model();
        converged = newton_small();
      }
      break;
    }
    if (r < rmin) { end_on_radius(); break; }
  }

  // fvalue, preddiff, step and the iteration limit end the run right after an
  // accepted step, which the snapshot at the head of the next iteration never
  // reaches. The other exits follow a rejection and are already covered.
  refresh_current();

  if (stop_reason == "objfun")
    Rf_warning("trust: objfun evaluation failed 3 times in a row");
  else if (stop_reason == "radius")
    Rf_warning("Trust radius fell below rmin or below the resolution of the parameters. "
               "Fit is not converged.");
  else if (!converged && n_iter >= iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  // The best iterate visited, which is the last one unless nonmonotone > 0.
  NumericVector arg_out(K), grad_out(best_grad.begin(), best_grad.end());
  for (int i = 0; i < K; ++i) arg_out[i] = best_z[i] / ps[i];
  arg_out.names()  = parnames;
  grad_out.names() = parnames;
  SEXP Hess_out = Model::hessian(best_H, K, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (best_at_bound[i] != 0);
  at_bound_out.names() = parnames;

  IntegerVector eval_by_source(kHessianSourceCount);
  CharacterVector src_names(kHessianSourceCount);
  for (int i = 0; i < kHessianSourceCount; ++i) {
    eval_by_source[i] = eval_src[i];
    src_names[i]      = kHessianSourceNames[i];
  }
  eval_by_source.names() = src_names;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = best_val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = n_iter,
      Named("neval")      = neval,
      Named("qnEval")     = qn_neval,
      Named("qnSkipped")  = qn_skipped,
      Named("nSwitch")    = n_switch,
      Named("nReseed")    = n_reseed,
      Named("hybridGN")   = n_hybrid_gn,
      Named("phaseIterations") = IntegerVector(phase_iters.begin(), phase_iters.end()),
      Named("evalBySource") = eval_by_source,
      Named("converged")  = converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = stop_reason);

  if (!Rf_isNull(Model::lowrank(best_H))) result["lowrank"] = Model::lowrank(best_H);
  if (blather_on) trace.attach(result, n_iter, K, parnames, minimize);
  return result;
}

// ---- "clip": active-set reduction plus componentwise clipping --------------
List trust_clip(Function objfun, NumericVector parinit,
                double rinit, double rmax,
                Nullable<NumericVector> parscale,
                int iterlim, double ftol, double mtol,
                bool minimize, bool blather_on,
                Nullable<NumericVector> parupper,
                Nullable<NumericVector> parlower,
                bool printIter, Nullable<CharacterVector> traceFile) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  const bool rescale = parscale.isNotNull();
  fill_parscale(parscale, ps, K, "trust");

  std::vector<double> theta(K);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    if (parinit[i] > pu[i])      { any_above = true; theta[i] = pu[i]; }
    else if (parinit[i] < pl[i]) { any_below = true; theta[i] = pl[i]; }
    else                          theta[i] = parinit[i];
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");

  std::vector<unsigned char> at_upper(K, 0), at_lower(K, 0);
  for (int i = 0; i < K; ++i) {
    at_upper[i] = (theta[i] >= pu[i]) ? 1 : 0;
    at_lower[i] = (theta[i] <= pl[i]) ? 1 : 0;
  }

  Reporter report;
  report.init(printIter, iterlim, traceFile, K, parnames);

  NumericVector x_named(theta.begin(), theta.end());
  x_named.names() = parnames;
  // Through eval_objfun like every other call: an objective that throws at
  // parinit should say "not feasible" rather than surface as a bare R error
  // from inside a conversion.
  List out_init;
  if (!eval_objfun(objfun, x_named, out_init, CV_GN))
    stop_parinit_failed();
  double val = as<double>(out_init["value"]);
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  if (!std::isfinite(val)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_full(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K, 0.0);
  if (!read_hessian(out_init, K, H_full))
    stop("objfun returned no Hessian at parinit, which boundary = \"clip\" "
         "needs at every iterate");

  int neval = 1;
  report(neval, val, x_named, /*head=*/true);

  Blather trace;
  std::vector<int>    active;
  std::vector<double> g_red, H_red, eigvals_red, eigvecs_red;
  double f_used = val;

  bool accept = true, converged = false, is_terminate = false;
  double r = rinit;
  int iter = 0, n_fail = 0;
  std::vector<double> theta_try(K), p_full(K), p_red;

  for (iter = 1; iter <= iterlim; ++iter) {

    if (blather_on) {
      trace.argpath.insert(trace.argpath.end(), theta.begin(), theta.end());
      trace.r.push_back(r);
      trace.valpath.push_back(val);
    }

    if (accept) {
      active.clear();
      for (int i = 0; i < K; ++i) {
        double gi = grad_full[i];
        bool drop;
        if (minimize) drop = (at_upper[i] && gi < 0.0) || (at_lower[i] && gi > 0.0);
        else          drop = (at_upper[i] && gi > 0.0) || (at_lower[i] && gi < 0.0);
        if (!drop) active.push_back(i);
      }
      int Kred = static_cast<int>(active.size());
      g_red.assign(Kred, 0.0);
      H_red.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        int i = active[ii];
        g_red[ii] = grad_full[i];
        for (int jj = 0; jj < Kred; ++jj) {
          int j = active[jj];
          H_red[ii + (std::size_t) jj * Kred] = H_full[i + (std::size_t) j * K];
        }
      }
      if (rescale) {
        for (int ii = 0; ii < Kred; ++ii) g_red[ii] /= ps[active[ii]];
        for (int jj = 0; jj < Kred; ++jj)
          for (int ii = 0; ii < Kred; ++ii)
            H_red[ii + (std::size_t) jj * Kred] /= ps[active[ii]] * ps[active[jj]];
      }
      f_used = val;
      if (!minimize) {
        for (auto& x : g_red) x = -x;
        for (auto& x : H_red) x = -x;
        f_used = -val;
      }
      if (Kred > 0) {
        eigvals_red.assign(Kred, 0.0);
        eigvecs_red.assign((std::size_t) Kred * Kred, 0.0);
        eigen_sym_local(H_red.data(), Kred, eigvals_red.data(), eigvecs_red.data());
      }
    }

    int Kred = static_cast<int>(active.size());

    p_red.assign(Kred, 0.0);
    bool is_newton = false, is_hard = false, is_easy = false;
    double m_value = 0.0, pred_pos = 0.0;
    if (Kred > 0) {
      trust_sub(Kred, g_red.data(), eigvals_red.data(), eigvecs_red.data(), r,
                p_red.data(), &pred_pos, &is_newton, &is_hard, &is_easy);
      m_value  = model_value(Kred, g_red.data(), H_red.data(), p_red.data());
      pred_pos = -m_value;
    }

    std::fill(p_full.begin(), p_full.end(), 0.0);
    for (int ii = 0; ii < Kred; ++ii) {
      int i = active[ii];
      double s = p_red[ii];
      if (rescale) s /= ps[i];
      p_full[i] = s;
    }
    double stepnorm = 0.0;
    for (int i = 0; i < K; ++i) stepnorm += p_full[i] * p_full[i];
    stepnorm = std::sqrt(stepnorm);

    for (int i = 0; i < K; ++i) theta_try[i] = theta[i] + p_full[i];
    for (int i = 0; i < K; ++i) {
      at_upper[i] = !(theta_try[i] < pu[i]) ? 1 : 0;
      at_lower[i] = !(theta_try[i] > pl[i]) ? 1 : 0;
      if (at_upper[i]) theta_try[i] = pu[i];
      if (at_lower[i]) theta_try[i] = pl[i];
    }

    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try, CV_GN);
    double val_try = kInf;
    NumericVector grad_try;
    std::vector<double> Htry;
    if (eval_ok) {
      val_try  = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      if (!read_hessian(out_try, K, Htry)) eval_ok = false;
      if (!std::isfinite(val_try)) eval_ok = false;
    }
    neval++;
    report(neval, val_try, x_try, /*head=*/false);

    double ftry_used = minimize ? val_try : -val_try;
    double rho;
    if (eval_ok && pred_pos > 0.0) rho = (f_used - ftry_used) / pred_pos;
    else                           rho = -kInf;

    is_terminate = eval_ok && (std::fabs(ftry_used - f_used) < ftol ||
                               std::fabs(m_value) < mtol);

    if (!eval_ok) {
      n_fail++;
      accept = false;
      r *= 0.25;
      if (n_fail >= 3) {
        Rf_warning("trust: objfun evaluation failed 3 times in a row");
        break;
      }
    } else {
      n_fail = 0;
      if (is_terminate) {
        accept = (ftry_used < f_used);
      } else if (rho < 0.25) {
        accept = false;
        r *= 0.25;
      } else {
        accept = true;
        if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
      }
    }

    if (accept && eval_ok) {
      for (int i = 0; i < K; ++i) theta[i] = theta_try[i];
      val = val_try;
      grad_full.assign(grad_try.begin(), grad_try.end());
      H_full = Htry;
    }

    if (blather_on) {
      trace.argtry.insert(trace.argtry.end(), theta_try.begin(), theta_try.end());
      trace.valtry.push_back(val_try);
      trace.accept.push_back(accept ? 1 : 0);
      trace.preddiff.push_back(m_value);
      trace.stepnorm.push_back(stepnorm);
      trace.rho.push_back(rho);
      trace.steptype.push_back(subproblem_label(is_newton, is_hard, is_easy));
      trace.stepback.push_back("full");
    }

    if (is_terminate) { converged = true; break; }
  }

  int final_iter = (iter <= iterlim) ? iter : iterlim;
  if (!converged && final_iter == iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  NumericVector arg_out(theta.begin(), theta.end());
  NumericVector grad_out(grad_full.begin(), grad_full.end());
  arg_out.names()  = parnames;
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) Hess_out(i, j) = H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (at_upper[i] || at_lower[i]);
  at_bound_out.names() = parnames;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = final_iter,
      Named("converged")  = converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = std::string(converged ? "fvalue" : "iterlim"));

  if (blather_on) trace.attach(result, final_iter, K, parnames, minimize);
  return result;
}

}  // namespace


// [[Rcpp::export]]
List trust_impl(Function objfun,
                NumericVector parinit,
                double rinit,
                double rmax,
                Nullable<NumericVector> parscale,
                int iterlim,
                double ftol,
                double mtol,
                double gtol,
                double xtol,
                double rmin,
                double thetamax,
                std::string boundary,
                CharacterVector phaseMethod,
                IntegerVector phaseIterlim,
                IntegerVector phaseRadiusStall,
                LogicalVector phaseReverse,
                LogicalVector phaseKeepRadius,
                LogicalVector phaseStopEnds,
                CharacterVector phaseInit,
                CharacterVector phaseReseed,
                IntegerVector phaseMemory,
                NumericVector phaseCautious,
                LogicalVector phaseRejected,
                LogicalVector phaseScaled,
                NumericVector phaseHybridTol,
                LogicalVector phaseHybridShift,
                int    switches = 0,
                double nonmonotone = 0.0,
                bool   sparse    = false,
                bool   minimize  = true,
                bool   blather   = false,
                Nullable<NumericVector>  parupper  = R_NilValue,
                Nullable<NumericVector>  parlower  = R_NilValue,
                bool   printIter = false,
                Nullable<CharacterVector> traceFile = R_NilValue) {

  const int K = parinit.size();
  if (K == 0) stop("trust: parinit must be non-empty");
  if (!parinit.hasAttribute("names"))
    stop("trust: parinit must be a named numeric vector");
  for (int i = 0; i < K; ++i)
    if (!std::isfinite(parinit[i])) stop("trust: parinit not all finite");

  const int P = phaseMethod.size();
  if (P < 1) stop("trust: at least one phase is needed");
  if (phaseIterlim.size() != P || phaseRadiusStall.size() != P || phaseReverse.size() != P ||
      phaseKeepRadius.size() != P || phaseStopEnds.size() != P ||
      phaseInit.size() != P || phaseReseed.size() != P ||
      phaseMemory.size() != P || phaseCautious.size() != P ||
      phaseRejected.size() != P || phaseScaled.size() != P ||
      phaseHybridTol.size() != P || phaseHybridShift.size() != P)
    stop("trust: the phase settings differ in length");
  if (switches < 0) stop("trust: switches must be >= 0");

  std::vector<Phase> phases(P);
  for (int k = 0; k < P; ++k) {
    Phase& ph = phases[k];
    const std::string m = as<std::string>(phaseMethod[k]);
    if      (m == "gn")    ph.method = HM_GN;
    else if (m == "bfgs")  ph.method = HM_BFGS;
    else if (m == "sr1")   ph.method = HM_SR1;
    else if (m == "exact") ph.method = HM_EXACT;
    else if (m == "sbfgs")   ph.method = HM_SBFGS;
    else if (m == "ssr1")    ph.method = HM_SSR1;
    else if (m == "gnsbfgs") ph.method = HM_GNSBFGS;
    else if (m == "gnssr1")  ph.method = HM_GNSSR1;
    else stop("trust: hessianMethod must be one of \"gn\", \"bfgs\", \"sr1\", "
              "\"exact\", \"sbfgs\", \"ssr1\", \"gnsbfgs\", \"gnssr1\"");
    // NA: no phase budget, only the total iterlim.
    ph.iterlim = (phaseIterlim[k] == NA_INTEGER) ? INT_MAX : phaseIterlim[k];
    if (ph.iterlim < 1) stop("trust: the iterlim of a phase must be >= 1");
    // NA: no handover on the radius.
    ph.radius_stall = (phaseRadiusStall[k] == NA_INTEGER) ? 0 : phaseRadiusStall[k];
    if (ph.radius_stall < 0) stop("trust: the radiusStall of a phase must be >= 1");
    ph.reverse = phaseReverse[k] == TRUE;
    ph.keep_radius = phaseKeepRadius[k] == TRUE;
    ph.stop_ends = phaseStopEnds[k] == TRUE;

    QnParams& q = ph.qn;
    const std::string init = as<std::string>(phaseInit[k]);
    const std::string reseed = as<std::string>(phaseReseed[k]);
    if      (init == "gn")       q.init = HI_GN;
    else if (init == "identity") q.init = HI_IDENTITY;
    else if (init == "exact")    q.init = HI_EXACT;
    else if (init == "tracked")  q.init = HI_TRACKED;
    else stop("trust: hessianInit must be one of \"gn\", \"identity\", "
              "\"exact\", \"tracked\"");
    if      (reseed == "never") q.reseed = HR_NEVER;
    else if (reseed == "stall") q.reseed = HR_STALL;
    else stop("trust: hessianReseed must be one of \"never\", \"stall\"");
    q.memory = phaseMemory[k];
    q.cautious = phaseCautious[k];
    q.rejected = phaseRejected[k] == TRUE;
    q.scaled = phaseScaled[k] == TRUE;
    q.hybrid_tol = phaseHybridTol[k];
    q.hybrid_shift = phaseHybridShift[k] == TRUE;
    if (q.memory < 0) stop("trust: qnMemory must be >= 0");
    if (!(q.hybrid_tol > 0.0)) stop("trust: qnHybridTol must be > 0");
    if (hm_structured(ph.method) && (q.memory != 0 || q.init != HI_GN))
      stop("trust: a structured hessianMethod starts from the Gauss-Newton Hessian and keeps "
           "no pairs; it takes qnControl$hessianInit = \"gn\" and qnMemory = 0");
    if (q.cautious < 0.0) stop("trust: qnCautious must be >= 0");
    if (q.init == HI_TRACKED && (k == 0 || ph.method != HM_BFGS || q.memory != 0 || sparse))
      stop("trust: hessianInit = \"tracked\" takes a \"bfgs\" phase after the first, with "
           "qnMemory = 0, on the dense path");
  }
  if (!(nonmonotone >= 0.0 && nonmonotone < 1.0))
    stop("trust: nonmonotone must be in [0, 1)");

  if (sparse)
    for (const Phase& ph : phases)
      if (ph.method != HM_GN && ph.method != HM_EXACT && !(ph.method == HM_BFGS && ph.qn.init == HI_GN))
        stop("trust(sparse = TRUE): the Hessian comes from the objective, \"gn\" or \"exact\", or "
             "from limited-memory BFGS on it, \"bfgs\" with hessianInit = \"gn\"");
  if (boundary == "clip") {
    if (sparse) stop("trust(sparse = TRUE) takes boundary = \"reflective\"");
    if (P != 1 || phases[0].method != HM_GN || switches > 0)
      stop("trust: boundary = \"clip\" takes a single \"gn\" phase; other "
           "Hessian sources and phases require boundary = \"reflective\"");
    return trust_clip(objfun, parinit, rinit, rmax, parscale, iterlim,
                      ftol, mtol, minimize, blather,
                      parupper, parlower, printIter, traceFile);
  }
  if (boundary != "reflective")
    stop("trust: boundary must be one of \"reflective\", \"clip\"");

  if (sparse)
    return trust_reflective<dmod::trust_model::SparseModel>(objfun, parinit, rinit, rmax, parscale, iterlim,
                          ftol, mtol, gtol, xtol, rmin, thetamax,
                          phases, switches, nonmonotone, minimize, blather,
                          parupper, parlower, printIter, traceFile);
  return trust_reflective<dmod::trust_model::DenseModel>(objfun, parinit, rinit, rmax, parscale, iterlim,
                          ftol, mtol, gtol, xtol, rmin, thetamax,
                          phases, switches, nonmonotone, minimize, blather,
                          parupper, parlower, printIter, traceFile);
}


// The cost of the factor of a sparse symmetric pattern, the upper triangle in compressed columns:
// flops of one numeric factorization and the entries of L, from the symbolic analysis alone.
// [[Rcpp::export]]
NumericVector sparse_factor_cost(IntegerVector p, IntegerVector i, int n) {
  dmod::trust_model::SparseModel::Hess h;
  h.p.assign(p.begin(), p.end());
  h.i.assign(i.begin(), i.end());
  dmod::trust_model::SparseModel model(n);
  model.pattern(h);
  return NumericVector::create(Named("flops") = model.c.fl, Named("lnz") = model.c.lnz);
}


// One subproblem step of the curvature models on H + U diag(d) U', `kernel` as the R side passes it
// (see trust_model.h) and `H` the same matrix dense: the step, the step kind, the predicted
// reduction and the reduction the full Newton step predicts.
// [[Rcpp::export]]
List trust_model_step(List kernel, NumericMatrix H, NumericVector g, double r, bool sparse) {
  const int K = g.size();
  const std::vector<double> ones(K, 1.0), zeros(K, 0.0);
  std::vector<double> s(K);
  bool newton = false, hard = false, easy = false;
  double pred = 0.0, newton_dec = 0.0;
  List out = List::create(Named("kernel") = kernel, Named("hessian") = H);
  if (sparse) {
    dmod::trust_model::SparseModel model(K);
    dmod::trust_model::SparseModel::Hess h;
    if (!dmod::trust_model::SparseModel::read(out, K, h)) stop("trust_model_step: no kernel");
    model.build(h, ones.data(), 1.0, zeros.data());
    model.solve(g.begin(), r, s.data(), &newton, &hard, &easy, &pred);
    newton_dec = model.newton_decrease(g.begin());
  } else {
    dmod::trust_model::DenseModel model(K);
    dmod::trust_model::DenseModel::Hess h;
    if (!dmod::trust_model::DenseModel::read(out, K, h)) stop("trust_model_step: no dense Hessian");
    model.build(h, ones.data(), 1.0, zeros.data());
    model.solve(g.begin(), r, s.data(), &newton, &hard, &easy, &pred);
    newton_dec = model.newton_decrease(g.begin());
  }
  return List::create(Named("p") = NumericVector(s.begin(), s.end()), Named("newton") = newton,
                      Named("hard") = hard, Named("easy") = easy, Named("pred") = pred,
                      Named("newtonDecrease") = newton_dec);
}
