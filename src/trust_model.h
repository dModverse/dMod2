// The curvature models of the trust-region kernels: Bhat = sgn W H W + diag(c) dense (eigendecomposition
// and trust_sub) or sparse plus a low-rank part U diag(d) U' (CHOLMOD and Woodbury), with the product
// Bhat x and the subproblem solved to rounding on both.

#ifndef DMOD_TRUST_MODEL_H
#define DMOD_TRUST_MODEL_H

#include <Rcpp.h>
#define R_MATRIX_INLINE inline
#include <Matrix/Matrix.h>
#include <Matrix/stubs.c>
#include "trust_subproblem.h"
#include "trust_driver.h"
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace dmod { namespace trust_model {

using namespace Rcpp;

// ---- dense --------------------------------------------------------------------------------

struct DenseModel {
  // The Hessian of an iterate.
  struct Hess {
    std::vector<double> H;
  };
  static constexpr bool dense = true;

  int K;
  std::vector<double> B, vals, vecs;

  explicit DenseModel(int K_) : K(K_), B((std::size_t) K_ * K_), vals(K_), vecs((std::size_t) K_ * K_) {}

  static bool read(const List& out, int K, Hess& h) {
    return dmod::trust_driver::read_hessian(out, K, h.H);
  }
  static void identity(int K, double sgn, Hess& h) {
    h.H.assign((std::size_t) K * K, 0.0);
    for (int i = 0; i < K; ++i) h.H[i + (std::size_t) i * K] = sgn;
  }

  void build(const Hess& h, const double* w, double sgn, const double* c) {
    // LAPACK does not return on NaN input
    for (double v : h.H) if (!std::isfinite(v)) throw std::runtime_error("the Hessian is not finite.");
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        B[i + (std::size_t) j * K] = w[i] * w[j] * sgn * h.H[i + (std::size_t) j * K];
    for (int i = 0; i < K; ++i) B[i + (std::size_t) i * K] += c[i];
    dmod::trust_internal::eigen_sym_local(B.data(), K, vals.data(), vecs.data());
  }

  void apply(const double* x, double* y) const { dmod::trust_internal::DenseOp{K, B.data()}.apply(x, y); }

  // also the predicted reduction -m(s) in `pred`
  void solve(const double* g, double r, double* s, bool* newton, bool* hard, bool* easy,
             double* pred = nullptr) {
    double pr = 0.0;
    dmod::trust_internal::trust_sub(K, g, vals.data(), vecs.data(), r, s, &pr, newton, hard, easy);
    if (pred) *pred = pr;
  }

  // the decrease the full Newton step on the positive part of Bhat predicts
  double newton_decrease(const double* g) const {
    return dmod::trust_internal::newton_decrease_eigen(K, g, vals.data(), vecs.data());
  }

  static SEXP hessian(const Hess& h, int K, const CharacterVector& names) {
    NumericMatrix M(K, K);
    std::copy(h.H.begin(), h.H.end(), M.begin());
    M.attr("dimnames") = List::create(names, names);
    return M;
  }
  static SEXP lowrank(const Hess&) { return R_NilValue; }
};

// ---- sparse -------------------------------------------------------------------------------

struct SparseModel {
  // The Hessian of an iterate: the upper triangle of S in compressed columns, the low-rank part
  // U (K x m, column-major) diag(d), and the R objects it came from, which the fit returns.
  struct Hess {
    std::vector<int> p, i;
    std::vector<double> x, U, d;
    RObject hessian, lowrank;
  };
  static constexpr bool dense = false;

  int K, m = 0;
  cholmod_common c;
  CHM_SP A = nullptr;          // Bhat_S, upper triangle with every diagonal entry
  CHM_FR L = nullptr;          // LL' of Bhat_S + lambda I
  CHM_FR Q = nullptr;          // LDL' of an indefinite Bhat_S + lambda I, analyzed where needed
  CHM_FR F = nullptr;          // the factor of the last factorization
  std::vector<int> pat_p, pat_i, diag;   // pattern of the last analysis, diagonal positions
  std::vector<int> map;                  // entry k of S to the position in A
  std::vector<double> V, E;              // the low-rank part of Bhat: V diag(E) V'
  // at the last factorization: lambda, Z = (Bhat_S + lambda I)^-1 V, Kinv = (diag(1/E) + V'Z)^-1
  double lam_f = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> Z, Kinv;
  double scale = 1.0;          // a bound of |Bhat|

  explicit SparseModel(int K_) : K(K_) {
    M_cholmod_start(&c);
    c.final_ll = 1;
    c.quick_return_if_not_posdef = 1;
    c.error_handler = nullptr;
    c.print = 0;
  }
  ~SparseModel() {
    if (Q) M_cholmod_free_factor(&Q, &c);
    if (L) M_cholmod_free_factor(&L, &c);
    if (A) M_cholmod_free_sparse(&A, &c);
    M_cholmod_finish(&c);
  }
  SparseModel(const SparseModel&) = delete;
  SparseModel& operator=(const SparseModel&) = delete;

  // `out$kernel` from the R side: list(p, i, x, U, d), the upper triangle of S over the
  // parameters and the low-rank part.
  static bool read(const List& out, int K, Hess& h) {
    if (!out.containsElementNamed("kernel")) return false;
    SEXP kr = out["kernel"];
    if (TYPEOF(kr) != VECSXP) return false;
    List k(kr);
    IntegerVector p = k["p"], i = k["i"];
    NumericVector x = k["x"], U = k["U"], d = k["d"];
    if (p.size() != K + 1 || i.size() != x.size()) return false;
    if (U.size() != (R_xlen_t) K * d.size()) return false;
    h.p.assign(p.begin(), p.end());
    h.i.assign(i.begin(), i.end());
    h.x.assign(x.begin(), x.end());
    h.U.assign(U.begin(), U.end());
    h.d.assign(d.begin(), d.end());
    h.hessian = out["hessian"];
    h.lowrank = out.containsElementNamed("lowrank") ? RObject(out["lowrank"]) : RObject(R_NilValue);
    return true;
  }
  static void identity(int, double, Hess&) {
    Rcpp::stop("trust(sparse = TRUE) takes Hessians of the objective; a quasi-Newton phase is dense");
  }

  // A and the symbolic analysis of A for the pattern of `h`; the map of the entries of S into A.
  void pattern(const Hess& h) {
    if (A && h.p == pat_p && h.i == pat_i) return;
    pat_p = h.p; pat_i = h.i;
    std::vector<int> Ap(K + 1, 0), Ai;
    Ai.reserve(h.i.size() + K);
    map.assign(h.i.size(), -1);
    diag.assign(K, -1);
    for (int j = 0; j < K; ++j) {
      bool has = false;
      for (int k = h.p[j]; k < h.p[j + 1]; ++k) {
        const int row = h.i[k];
        if (row > j) continue;
        if (row == j) has = true;
        map[k] = (int) Ai.size();
        if (row == j) diag[j] = (int) Ai.size();
        Ai.push_back(row);
      }
      if (!has) { diag[j] = (int) Ai.size(); Ai.push_back(j); }
      Ap[j + 1] = (int) Ai.size();
    }
    if (Q) M_cholmod_free_factor(&Q, &c);
    if (L) M_cholmod_free_factor(&L, &c);
    if (A) M_cholmod_free_sparse(&A, &c);
    A = M_cholmod_allocate_sparse(K, K, Ai.size(), 1, 1, 1, CHOLMOD_REAL, &c);
    std::copy(Ap.begin(), Ap.end(), (int*) A->p);
    std::copy(Ai.begin(), Ai.end(), (int*) A->i);
    std::fill((double*) A->x, (double*) A->x + Ai.size(), 0.0);
    L = M_cholmod_analyze(A, &c);
    if (!L) Rcpp::stop("trust(sparse = TRUE): the analysis of the sparse Hessian failed");
  }

  void build(const Hess& h, const double* w, double sgn, const double* cd) {
    pattern(h);
    double* Ax = (double*) A->x;
    const int* Ap = (const int*) A->p;
    const int* Ai = (const int*) A->i;
    std::fill(Ax, Ax + Ap[K], 0.0);
    for (int j = 0; j < K; ++j)
      for (int k = h.p[j]; k < h.p[j + 1]; ++k)
        if (map[k] >= 0) Ax[map[k]] += sgn * w[h.i[k]] * w[j] * h.x[k];
    for (int j = 0; j < K; ++j) Ax[diag[j]] += cd[j];
    m = (int) h.d.size();
    V.assign((std::size_t) K * m, 0.0);
    E.assign(m, 0.0);
    for (int l = 0; l < m; ++l) {
      E[l] = sgn * h.d[l];
      for (int i = 0; i < K; ++i) V[i + (std::size_t) l * K] = w[i] * h.U[i + (std::size_t) l * K];
    }
    // |Bhat| <= max column sum of |Bhat_S| + sum_l |E_l| |V_l|^2
    std::vector<double> col(K, 0.0);
    for (int j = 0; j < K; ++j)
      for (int k = Ap[j]; k < Ap[j + 1]; ++k) {
        col[j] += std::fabs(Ax[k]);
        if (Ai[k] != j) col[Ai[k]] += std::fabs(Ax[k]);
      }
    scale = *std::max_element(col.begin(), col.end());
    for (int l = 0; l < m; ++l) {
      double vv = 0.0;
      for (int i = 0; i < K; ++i) vv += V[i + (std::size_t) l * K] * V[i + (std::size_t) l * K];
      scale += std::fabs(E[l]) * vv;
    }
    if (!(scale > 0.0)) scale = 1.0;
    lam_f = std::numeric_limits<double>::quiet_NaN();
  }

  void apply(const double* x, double* y) const {
    const int* Ap = (const int*) A->p;
    const int* Ai = (const int*) A->i;
    const double* Ax = (const double*) A->x;
    for (int i = 0; i < K; ++i) y[i] = 0.0;
    for (int j = 0; j < K; ++j)
      for (int k = Ap[j]; k < Ap[j + 1]; ++k) {
        const int i = Ai[k];
        y[i] += Ax[k] * x[j];
        if (i != j) y[j] += Ax[k] * x[i];
      }
    for (int l = 0; l < m; ++l) {
      const double* Vl = V.data() + (std::size_t) l * K;
      double t = 0.0;
      for (int i = 0; i < K; ++i) t += Vl[i] * x[i];
      t *= E[l];
      for (int i = 0; i < K; ++i) y[i] += t * Vl[i];
    }
  }

  // (Bhat_S + lambda I)^-1 b, in place over the columns of the K x n matrix b.
  void solveS(double* b, int n) {
    CHM_DN B = M_cholmod_allocate_dense(K, n, K, CHOLMOD_REAL, &c);
    std::copy(b, b + (std::size_t) K * n, (double*) B->x);
    CHM_DN X = M_cholmod_solve(CHOLMOD_A, F, B, &c);
    std::copy((double*) X->x, (double*) X->x + (std::size_t) K * n, b);
    M_cholmod_free_dense(&X, &c);
    M_cholmod_free_dense(&B, &c);
  }

  // The LDL' factor of Bhat_S + lambda I and the number of negative pivots; -1 for a zero pivot.
  int factorLDL(double lambda) {
    if (!Q) {
      const int super = c.supernodal, ll = c.final_ll;
      c.supernodal = CHOLMOD_SIMPLICIAL;
      c.final_ll = 0;
      Q = M_cholmod_analyze(A, &c);
      c.supernodal = super;
      c.final_ll = ll;
      if (!Q) return -1;
    }
    double beta[2] = {lambda, 0.0};
    const int quick = c.quick_return_if_not_posdef, ll = c.final_ll;
    c.quick_return_if_not_posdef = 0;
    c.final_ll = 0;
    M_cholmod_factorize_p(A, beta, nullptr, 0, Q, &c);
    c.quick_return_if_not_posdef = quick;
    c.final_ll = ll;
    if (c.status != CHOLMOD_OK || Q->minor < (size_t) K || Q->is_ll) { c.status = 0; return -1; }
    const int* Qp = (const int*) Q->p;
    const double* Qx = (const double*) Q->x;
    double dmax = 0.0;
    for (int j = 0; j < K; ++j) dmax = std::max(dmax, std::fabs(Qx[Qp[j]]));
    int neg = 0;
    for (int j = 0; j < K; ++j) {
      const double dj = Qx[Qp[j]];
      if (!(std::fabs(dj) > 1e-14 * dmax)) return -1;
      neg += dj < 0.0;
    }
    return neg;
  }

  // Factors of Bhat + lambda I, false where it is not positive definite by the inertia neg(S) +
  // pos(Kmat) - pos(E), Kmat = diag(1/E) + V'S^-1 V, S = Bhat_S + lambda I by Cholesky or, where S is
  // indefinite and E has a positive entry, by LDL'.
  bool factor(double lambda) {
    if (lambda == lam_f) return true;
    lam_f = std::numeric_limits<double>::quiet_NaN();
    double beta[2] = {lambda, 0.0};
    M_cholmod_factorize_p(A, beta, nullptr, 0, L, &c);
    int negS = 0;
    F = L;
    if (c.status == CHOLMOD_NOT_POSDEF || L->minor < (size_t) K) {
      c.status = 0;
      if (std::none_of(E.begin(), E.begin() + m, [](double e) { return e > 0.0; })) return false;
      negS = factorLDL(lambda);
      if (negS < 0) return false;
      F = Q;
    }
    if (m == 0 && negS > 0) return false;
    if (m > 0) {
      Z = V;
      solveS(Z.data(), m);
      std::vector<double> Km((std::size_t) m * m, 0.0);
      for (int a = 0; a < m; ++a)
        for (int b = 0; b < m; ++b) {
          double t = 0.0;
          for (int i = 0; i < K; ++i) t += V[i + (std::size_t) a * K] * Z[i + (std::size_t) b * K];
          Km[a + (std::size_t) b * m] = t + (a == b ? 1.0 / E[a] : 0.0);
        }
      std::vector<double> ev(m), evec((std::size_t) m * m);
      dmod::trust_internal::eigen_sym_local(Km.data(), m, ev.data(), evec.data());
      double amax = 0.0;
      int pos = 0, posE = 0;
      for (int a = 0; a < m; ++a) {
        amax = std::max(amax, std::fabs(ev[a]));
        pos += ev[a] > 0.0;
        posE += E[a] > 0.0;
      }
      for (int a = 0; a < m; ++a) if (!(std::fabs(ev[a]) > 1e-14 * amax)) return false;
      if (negS + pos != posE) return false;
      Kinv.assign((std::size_t) m * m, 0.0);
      for (int a = 0; a < m; ++a)
        for (int b = 0; b < m; ++b) {
          double t = 0.0;
          for (int e = 0; e < m; ++e) t += evec[a + (std::size_t) e * m] * evec[b + (std::size_t) e * m] / ev[e];
          Kinv[a + (std::size_t) b * m] = t;
        }
    }
    lam_f = lambda;
    return true;
  }

  // (Bhat + lambda I)^-1 b at the last factorization, in place.
  void solveB(double* b) {
    solveS(b, 1);
    if (m == 0) return;
    std::vector<double> t(m, 0.0), u(m, 0.0);
    for (int a = 0; a < m; ++a)
      for (int i = 0; i < K; ++i) t[a] += V[i + (std::size_t) a * K] * b[i];
    for (int a = 0; a < m; ++a)
      for (int e = 0; e < m; ++e) u[a] += Kinv[a + (std::size_t) e * m] * t[e];
    for (int a = 0; a < m; ++a)
      for (int i = 0; i < K; ++i) b[i] -= Z[i + (std::size_t) a * K] * u[a];
  }

  // the Euclidean norm, scaled against under- and overflow
  static double norm(const double* x, int n) {
    double m = 0.0;
    for (int i = 0; i < n; ++i) m = std::max(m, std::fabs(x[i]));
    if (m == 0.0 || !std::isfinite(m)) return m;
    double s = 0.0;
    for (int i = 0; i < n; ++i) { const double t = x[i] / m; s += t * t; }
    return m * std::sqrt(s);
  }
  static double norm(const std::vector<double>& x) { return norm(x.data(), (int) x.size()); }

  // Moré-Sorensen: the Newton step inside the region, else |s(lambda)| = r by safeguarded Newton on
  // 1/|s|; a bracket that closes with |s| < r is the hard case, extended to the boundary along the
  // eigenvector of lambda_min from inverse iteration unless lambda_min is zero to rounding.
  void solve(const double* g, double r, double* s_out, bool* newton, bool* hard, bool* easy,
             double* pred = nullptr) {
    solveStep(g, r, s_out, newton, hard, easy);
    if (pred) *pred = -dmod::trust_internal::model_value(K, g, *this, s_out);
  }

  // The decrease the full Newton step on the positive part of Bhat predicts: from the factorization
  // where Bhat is positive definite, from the dense eigendecomposition of Bhat otherwise.
  double newton_decrease(const double* g) {
    if (!factor(0.0)) {
      std::vector<double> B((std::size_t) K * K), e(K, 0.0), vals(K), vecs((std::size_t) K * K);
      for (int j = 0; j < K; ++j) {
        e[j] = 1.0;
        apply(e.data(), B.data() + (std::size_t) j * K);
        e[j] = 0.0;
      }
      dmod::trust_internal::eigen_sym_local(B.data(), K, vals.data(), vecs.data());
      return dmod::trust_internal::newton_decrease_eigen(K, g, vals.data(), vecs.data());
    }
    std::vector<double> s(g, g + K);
    solveB(s.data());
    double d = 0.0;
    for (int i = 0; i < K; ++i) d += g[i] * s[i];
    return 0.5 * d;
  }

  void solveStep(const double* g, double r, double* s_out, bool* newton, bool* hard, bool* easy) {
    *newton = *hard = *easy = false;
    std::vector<double> s(K), w(K);
    auto step = [&](double lambda) {
      if (!factor(lambda)) return false;
      for (int i = 0; i < K; ++i) s[i] = -g[i];
      solveB(s.data());
      return true;
    };
    if (step(0.0) && norm(s) <= r) {
      std::copy(s.begin(), s.end(), s_out);
      *newton = true;
      return;
    }
    const double gn = norm(g, K);
    if (gn == 0.0) { std::fill(s_out, s_out + K, 0.0); *newton = true; return; }
    const double eps = 8.0 * K * std::numeric_limits<double>::epsilon() * scale;
    double lo = 0.0, hi = gn / r + scale;
    double lam = std::max(eps, gn / r);
    bool done = false;
    double lam_in = -1.0;   // the smallest lambda known with |s| < r
    for (int it = 0; it < 300 && !done; ++it) {
      if (!(lam > lo && lam < hi)) lam = (lo > 0.0) ? std::sqrt(lo * hi) : 0.5 * hi;
      if (!(lam > lo) || hi - lo <= 1e-15 * hi) break;
      if (!step(lam)) { lo = lam; lam = 0.0; continue; }
      const double ns = norm(s);
      if (std::fabs(ns - r) <= 1e-12 * r) { done = true; break; }
      if (ns < r) { hi = lam; lam_in = lam; } else lo = lam;
      // Newton on 1/|s| = 1/r: d|s|/dlambda = -s'(B + lambda I)^-1 s / |s|
      w = s;
      solveB(w.data());
      double q = 0.0;
      for (int i = 0; i < K; ++i) q += s[i] * w[i];
      lam = lam + (ns - r) / r * ns * ns / q;
    }
    // a bracket that closes before the root resolves is the hard case to rounding: the step at the
    // upper end of the bracket lies inside the region
    if (!done) step(lam_in > 0.0 ? lam_in : hi);
    std::copy(s.begin(), s.end(), s_out);
    const double ns = norm(s);
    if (ns >= r * (1.0 - 1e-10)) { *easy = true; return; }
    // the hard case: s inside the region at lambda next to -lambda_min
    *hard = true;
    if (lam_f <= eps) return;
    std::vector<double> z(K, 1.0);
    for (int i = 0; i < K; ++i) z[i] = 1.0 + 1e-3 * i;
    for (int it = 0; it < 30; ++it) {
      solveB(z.data());
      const double nz = norm(z);
      for (int i = 0; i < K; ++i) z[i] /= nz;
    }
    // s + tau z on the boundary, the root of the smaller model value
    double sz = 0.0;
    for (int i = 0; i < K; ++i) sz += s[i] * z[i];
    const double disc = std::sqrt(std::max(0.0, sz * sz + r * r - ns * ns));
    const double t1 = -sz + disc, t2 = -sz - disc;
    std::vector<double> a(K), b(K);
    for (int i = 0; i < K; ++i) { a[i] = s[i] + t1 * z[i]; b[i] = s[i] + t2 * z[i]; }
    const double ma = dmod::trust_internal::model_value(K, g, *this, a.data());
    const double mb = dmod::trust_internal::model_value(K, g, *this, b.data());
    const std::vector<double>& best = (ma <= mb) ? a : b;
    std::copy(best.begin(), best.end(), s_out);
  }

  static SEXP hessian(const Hess& h, int, const CharacterVector&) { return h.hessian; }
  static SEXP lowrank(const Hess& h) { return h.lowrank; }
};

}}  // namespace dmod::trust_model

#endif  // DMOD_TRUST_MODEL_H
