// Shared subproblem helpers for the trust-region optimisers in
// trust_kernel.cpp and trustL1_kernel.cpp.
//
// `trust_sub` is a faithful port of the Moré-Sorensen subproblem solver in
// the historical R `trust()` (R/trust.R). The algorithm parametrises the
// Lagrange multiplier as `beep = sigma + lam_min`, so the smallest shifted
// eigenvalue `beta_j = vals_j - lam_min` is always non-negative and the
// degenerate index set is `{j : beta_j <= eig_tol}`, eig_tol being the noise
// floor of the eigen decomposition. `q` is zeroed on that set below the same
// floor, so `C2 == 0` means "g orthogonal to the min-eigenspace" as decided in
// double precision rather than as the LAPACK build happened to round it.
// Branch selection then follows R's classical Moré-Sorensen formulation:
//   Newton              all(vals > eig_tol) and ||H^{-1} g|| <= r
//   easy (incl hard-easy)   C2 > 0  OR  C1 > r^2     (sigma found via root)
//   hard-hard           C2 == 0 AND C1 <= r^2        (step lands on the
//                                                     min-eigenspace boundary)
// where C1 = sum_{j !in imin}(q_j / beta_j)^2, C2 = sum_{j in imin} q_j^2.
//
// BLAS: `q = vecs^T g` and the back-projection `p = +/- vecs * y` go through
// dgemv. Eigen decomposition uses LAPACK dsyevr.
//
// The Coleman-Li (1994, 1996) layer below poses that same subproblem in a
// scaled frame, so a box constraint acts as a change of metric rather than a
// projection and iterates stay strictly interior:
//
//     ghat = D g      Bhat = D H D + C      s = D shat
//
// with D = diag(|v|^{1/2}) and C = diag(|g| * jv) >= 0, so Bhat is positive
// semi-definite whenever H is. Candidate steps are scored with `model_value`,
// the model the subproblem itself minimises, which keeps the predicted and the
// achieved reduction referring to one and the same step.

#ifndef DMOD_TRUST_SUBPROBLEM_H
#define DMOD_TRUST_SUBPROBLEM_H

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

#ifndef FCONE
#define FCONE
#endif


namespace dmod { namespace trust_internal {

// Symmetric K x K eigen via LAPACK dsyevr. `vals` (K, ascending),
// `vecs` (K*K, column-major). A is read-only (copied internally).
inline void eigen_sym_local(const double* A, int K, double* vals, double* vecs) {
  std::vector<double> A_copy(A, A + (std::size_t) K * K);
  std::vector<int>    isuppz(2 * K);
  int info = 0, m_out = 0;
  double abstol = 0.0;
  double wkopt;
  int    iwkopt;
  int    lwork  = -1;
  int    liwork = -1;
  F77_CALL(dsyevr)("V", "A", "U", &K, A_copy.data(), &K,
                   NULL, NULL, NULL, NULL, &abstol, &m_out,
                   vals, vecs, &K, isuppz.data(),
                   &wkopt, &lwork, &iwkopt, &liwork, &info FCONE FCONE FCONE);
  if (info != 0) throw std::runtime_error("dsyevr workspace query failed");
  lwork  = static_cast<int>(wkopt);
  liwork = iwkopt;
  std::vector<double> work(lwork);
  std::vector<int>    iwork(liwork);
  std::copy(A, A + (std::size_t) K * K, A_copy.begin());
  F77_CALL(dsyevr)("V", "A", "U", &K, A_copy.data(), &K,
                   NULL, NULL, NULL, NULL, &abstol, &m_out,
                   vals, vecs, &K, isuppz.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE FCONE);
  if (info != 0) throw std::runtime_error("dsyevr failed");
}

// Moré-Sorensen trust-region subproblem.
// Minimize g^T p + 1/2 p^T H p subject to ||p|| <= r, where H is given
// implicitly by its eigendecomposition (vals ascending, vecs column-major).
inline void trust_sub(int K, const double* g,
                      const double* vals, const double* vecs,
                      double r, double* p_out, double* predicted_red,
                      bool* is_newton, bool* is_hard, bool* is_easy) {
  *is_newton = false;
  *is_hard   = false;
  *is_easy   = false;

  const int    i_one  = 1;
  const double d_one  = 1.0;
  const double d_zero = 0.0;
  const double d_mone = -1.0;

  // q = vecs^T * g     (BLAS dgemv, T mode)
  std::vector<double> q(K);
  F77_CALL(dgemv)("T", &K, &K, &d_one, vecs, &K,
                  g, &i_one, &d_zero, q.data(), &i_one FCONE);

  // Noise floor of the eigen decomposition. Eigenvalues that agree to within it
  // are the same eigenvalue, and one that sits inside it is a zero of
  // unresolvable sign -- both differ between LAPACK builds, which is what made
  // the branch selection below platform-dependent.
  double lam_absmax = 0.0;
  for (int j = 0; j < K; ++j)
    if (std::fabs(vals[j]) > lam_absmax) lam_absmax = std::fabs(vals[j]);
  const double eig_tol =
    8.0 * K * std::numeric_limits<double>::epsilon() * lam_absmax;

  // Newton: take if H is positive definite in double precision AND
  // ||H^{-1} g|| <= r. A singular H whose zero eigenvalue happens to come back
  // just above zero would otherwise pass, and q_j / vals_j then amplifies the
  // roundoff in q_j into a full-size step along the null direction.
  bool all_pos = true;
  for (int j = 0; j < K; ++j) {
    if (!(vals[j] > eig_tol)) { all_pos = false; break; }
  }
  if (all_pos) {
    std::vector<double> y(K);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      y[j] = -q[j] / vals[j];
      pn2 += y[j] * y[j];
    }
    if (std::sqrt(pn2) <= r) {
      F77_CALL(dgemv)("N", &K, &K, &d_one, vecs, &K,
                      y.data(), &i_one, &d_zero, p_out, &i_one FCONE);
      double pred = 0.0;
      for (int j = 0; j < K; ++j) pred += 0.5 * q[j] * q[j] / vals[j];
      *predicted_red = pred;
      *is_newton = true;
      return;
    }
  }

  // Non-Newton path: shift eigenvalues so the smallest becomes 0.
  double lam_min = vals[0];
  for (int j = 1; j < K; ++j) if (vals[j] < lam_min) lam_min = vals[j];

  // The min-eigenspace is a cluster, not one index: a repeated smallest
  // eigenvalue comes back split by O(eps * ||H||), and dividing q_j by that
  // split would turn roundoff into a large C1 term.
  std::vector<double> beta(K);
  std::vector<unsigned char> imin(K);
  for (int j = 0; j < K; ++j) {
    beta[j] = vals[j] - lam_min;
    imin[j] = (beta[j] <= eig_tol) ? 1u : 0u;
    if (imin[j]) beta[j] = 0.0;
  }

  // C3: ||q||^2 = ||g||^2 (orthonormal eigenbasis)
  double C3 = 0.0;
  for (int j = 0; j < K; ++j) C3 += q[j] * q[j];

  // A gradient orthogonal to the min-eigenspace still projects onto it as
  // O(eps * ||g||), because the eigenvector carrying it is only that accurate.
  // Left in, the residue makes C2 > 0, routes the hard case through the easy
  // branch, and the root it then chases is of the size of the residue itself --
  // an arbitrary step along the null direction.
  const double q_tol =
    8.0 * K * std::numeric_limits<double>::epsilon() * std::sqrt(C3);
  for (int j = 0; j < K; ++j)
    if (imin[j] && std::fabs(q[j]) <= q_tol) q[j] = 0.0;

  // C1: contribution from non-degenerate eigenvectors at beep -> 0
  // C2: gradient mass in the min-eigenspace
  double C1 = 0.0, C2 = 0.0;
  for (int j = 0; j < K; ++j) {
    if (imin[j]) {
      C2 += q[j] * q[j];
    } else {
      double t = q[j] / beta[j];
      C1 += t * t;
    }
  }

  std::vector<double> w(K, 0.0);

  if (C2 > 0.0 || C1 > r * r) {
    // Easy / hard-easy: solve for beep on (0, infty) such that ||p|| == r.
    *is_easy = true;
    *is_hard = (C2 == 0.0);

    // fred(beep) = sqrt(1 / sum_j (q_j / (beta_j + beep))^2) - 1/r
    // monotonically increasing in beep on (0, infty).
    auto fred = [&](double beep) -> double {
      if (beep == 0.0) {
        if (C2 > 0.0) return -1.0 / r;
        return std::sqrt(1.0 / C1) - 1.0 / r;
      }
      double s = 0.0;
      for (int j = 0; j < K; ++j) {
        double d = beta[j] + beep;
        double t = q[j] / d;
        s += t * t;
      }
      return std::sqrt(1.0 / s) - 1.0 / r;
    };

    // Bracket [beta_dn, beta_up] from R's trust(): conservative outer bounds
    // derived from C2 and C3 = ||g||^2.
    double beta_dn = std::sqrt(C2) / r;
    double beta_up = std::sqrt(C3) / r;

    double root;
    double f_up = fred(beta_up);
    double f_dn = fred(beta_dn);
    if (f_up <= 0.0) {
      root = beta_up;
    } else if (f_dn >= 0.0) {
      root = beta_dn;
    } else {
      // Bisection. fred is monotone increasing, so this converges robustly
      // to within ~50 iterations even when beta_dn / beta_up span many
      // orders of magnitude.
      double a = beta_dn, b = beta_up;
      for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (a + b);
        if (!(mid > a && mid < b)) { a = mid; b = mid; break; }
        double fm = fred(mid);
        if (fm > 0.0) b = mid; else a = mid;
        // Relative: the root can sit orders of magnitude below beta_up, and an
        // absolute floor would stop the bisection far short of it.
        if ((b - a) <= 1e-14 * std::fabs(b)) break;
      }
      root = 0.5 * (a + b);
    }

    for (int j = 0; j < K; ++j) w[j] = q[j] / (beta[j] + root);

    // p_out = -vecs * w
    F77_CALL(dgemv)("N", &K, &K, &d_mone, vecs, &K,
                    w.data(), &i_one, &d_zero, p_out, &i_one FCONE);

    double m_val = 0.0;
    for (int j = 0; j < K; ++j) {
      double y_j = -w[j];
      m_val += q[j] * y_j + 0.5 * vals[j] * y_j * y_j;
    }
    *predicted_red = -m_val;
    return;
  }

  // Hard-hard: gradient orthogonal to min-eigenspace AND the off-min
  // pseudo-inverse step already fits inside the trust region. Take that
  // step and extend along the min-eigenspace direction to the boundary.
  *is_hard = true;
  *is_easy = false;

  for (int j = 0; j < K; ++j) w[j] = imin[j] ? 0.0 : (q[j] / beta[j]);

  // p_out = -vecs * w
  F77_CALL(dgemv)("N", &K, &K, &d_mone, vecs, &K,
                  w.data(), &i_one, &d_zero, p_out, &i_one FCONE);

  // Extending along the min-eigendirection only pays against negative
  // curvature. For lam_min >= -eig_tol that direction carries no model information
  // (q_jmin == 0 by C2 == 0), so the minimum-norm step is the answer.
  double utry = 0.0;
  if (lam_min < -eig_tol) {
    double pn2 = 0.0;
    for (int i = 0; i < K; ++i) pn2 += p_out[i] * p_out[i];
    utry = std::sqrt(std::max(0.0, r * r - pn2));
    int jmin = -1;
    for (int j = 0; j < K; ++j) if (imin[j]) { jmin = j; break; }
    if (utry > 0.0 && jmin >= 0) {
      for (int i = 0; i < K; ++i) p_out[i] += utry * vecs[i + (std::size_t) jmin * K];
    } else {
      utry = 0.0;
    }
  }

  double m_val = 0.0;
  for (int j = 0; j < K; ++j) {
    // y_j eigen-frame coords: -w[j] for non-min, +utry on the chosen min idx.
    double y_j = -w[j];
    m_val += q[j] * y_j + 0.5 * vals[j] * y_j * y_j;
  }
  // Extension contribution: q_j == 0 on the min-eigenspace by construction, so
  // only the curvature term survives.
  m_val += 0.5 * lam_min * utry * utry;
  *predicted_red = -m_val;
}

// ---------------------------------------------------------------------------
// Coleman-Li interior trust-region-reflective boundary layer
// ---------------------------------------------------------------------------

// |v_i| is the distance to the bound the gradient pushes toward, or 1 when that
// direction is unbounded; jv_i = |d|v_i|/dtheta_i| is 1 exactly when that bound
// is finite. Mirrors fides' get_affine_scaling. Note jv_i * g_i = |g_i|, which
// is why the curvature correction C = diag(|g| * jv) is positive semi-definite.
inline void affine_scaling(int K, const double* theta, const double* g,
                           const double* lb, const double* ub,
                           double* absv, double* jv) {
  for (int i = 0; i < K; ++i) {
    const double b = (g[i] < 0.0) ? ub[i] : lb[i];
    if (std::isfinite(b)) {
      absv[i] = std::fabs(theta[i] - b);
      jv[i]   = 1.0;
    } else {
      absv[i] = 1.0;
      jv[i]   = 0.0;
    }
  }
}

// m(shat) = ghat^T shat + 1/2 shat^T Bhat shat. Bhat is column-major K x K.
inline double model_value(int K, const double* ghat, const double* Bhat,
                          const double* shat) {
  double gs = 0.0, sBs = 0.0;
  for (int i = 0; i < K; ++i) {
    gs += ghat[i] * shat[i];
    double row = 0.0;
    for (int j = 0; j < K; ++j)
      row += Bhat[i + (std::size_t) j * K] * shat[j];
    sBs += shat[i] * row;
  }
  return gs + 0.5 * sBs;
}

// Largest alpha >= 0 with lb <= base + alpha * dir <= ub, +Inf if nothing
// blocks. When `blocked` is non-null it is filled with the coordinates
// attaining that alpha.
inline double max_feasible_fraction(int K, const double* base, const double* dir,
                                    const double* lb, const double* ub,
                                    unsigned char* blocked) {
  double alpha = std::numeric_limits<double>::infinity();
  for (int i = 0; i < K; ++i) {
    if (dir[i] == 0.0) continue;
    const double bound = (dir[i] > 0.0) ? ub[i] : lb[i];
    if (!std::isfinite(bound)) continue;
    const double a = (bound - base[i]) / dir[i];
    if (a < alpha) alpha = a;
  }
  if (alpha < 0.0) alpha = 0.0;
  if (blocked != nullptr) {
    const double tol = 1e-10 * (1.0 + std::fabs(alpha));
    for (int i = 0; i < K; ++i) {
      blocked[i] = 0;
      if (!std::isfinite(alpha) || dir[i] == 0.0) continue;
      const double bound = (dir[i] > 0.0) ? ub[i] : lb[i];
      if (!std::isfinite(bound)) continue;
      if ((bound - base[i]) / dir[i] <= alpha + tol) blocked[i] = 1;
    }
  }
  return alpha;
}

// Largest t >= 0 with ||shat0 + t * d|| <= r.
inline double tr_fraction(int K, const double* shat0, const double* d, double r) {
  double a = 0.0, b = 0.0, c = 0.0;
  for (int i = 0; i < K; ++i) {
    a += d[i] * d[i];
    b += shat0[i] * d[i];
    c += shat0[i] * shat0[i];
  }
  c -= r * r;
  if (!(a > 0.0) || c > 0.0) return 0.0;
  return (-b + std::sqrt(std::max(0.0, b * b - a * c))) / a;
}

// argmin over t in [0, tmax] of m(shat0 + t * d).
inline double line_min(int K, const double* ghat, const double* Bhat,
                       const double* shat0, const double* d, double tmax) {
  if (!(tmax > 0.0)) return 0.0;
  double a = 0.0, b = 0.0;   // a = d^T Bhat d,  b = (ghat + Bhat shat0)^T d
  for (int i = 0; i < K; ++i) {
    double Bd = 0.0, Bs = 0.0;
    for (int j = 0; j < K; ++j) {
      const double Bij = Bhat[i + (std::size_t) j * K];
      Bd += Bij * d[j];
      Bs += Bij * shat0[j];
    }
    a += d[i] * Bd;
    b += (ghat[i] + Bs) * d[i];
  }
  if (a > 0.0) {
    const double t = -b / a;
    return (t < 0.0) ? 0.0 : ((t > tmax) ? tmax : t);
  }
  return (b < 0.0) ? tmax : 0.0;
}

// Choose the best in-box step. `shat` is the trust-region solution in the
// scaled frame (||shat|| <= r); `sqrt_absv` is diag(D). Candidates -- the
// truncated step, its single reflection off the blocking faces, and the scaled
// steepest-descent step -- are all scored with `model_value`, so `*m_out` is
// the predicted change for exactly the step returned in `shat_out` / `s_out`.
// `theta_frac` in (0, 1] is the fraction of the distance to a face that may be
// used, and is what keeps the iterate strictly interior.
inline void stepback(int K, const double* theta,
                     const double* lb, const double* ub,
                     const double* sqrt_absv,
                     const double* ghat, const double* Bhat,
                     const double* shat, double r, double theta_frac,
                     double* shat_out, double* s_out, double* m_out,
                     const char** label) {
  std::vector<double> s(K), cand_shat(K), cand_s(K), dir(K), dir_par(K), base(K);
  std::vector<unsigned char> blocked(K, 0);

  auto to_par = [&](const double* sh, double* out) {
    for (int i = 0; i < K; ++i) out[i] = sqrt_absv[i] * sh[i];
  };
  auto consider = [&](const char* name) {
    const double m = model_value(K, ghat, Bhat, cand_shat.data());
    if (m < *m_out) {
      *m_out = m;
      for (int i = 0; i < K; ++i) { shat_out[i] = cand_shat[i]; s_out[i] = cand_s[i]; }
      *label = name;
    }
  };

  to_par(shat, s.data());
  const double alpha = max_feasible_fraction(K, theta, s.data(), lb, ub, blocked.data());

  // The trust-region step itself, truncated to stay strictly interior.
  const double f = (alpha > 1.0) ? 1.0 : theta_frac * alpha;
  for (int i = 0; i < K; ++i) shat_out[i] = f * shat[i];
  to_par(shat_out, s_out);
  *m_out = model_value(K, ghat, Bhat, shat_out);
  *label = (f < 1.0) ? "truncated" : "full";

  // shat already minimises m over the whole ball; no candidate can beat it.
  if (f >= 1.0) return;

  // Single reflection: walk to the blocking face, flip the blocked components,
  // then minimise along the reflected direction within box and trust region.
  if (std::isfinite(alpha)) {
    for (int i = 0; i < K; ++i) {
      cand_shat[i] = alpha * shat[i];
      dir[i] = blocked[i] ? -shat[i] : shat[i];
    }
    to_par(cand_shat.data(), base.data());
    for (int i = 0; i < K; ++i) base[i] += theta[i];
    to_par(dir.data(), dir_par.data());
    const double t_max = std::min(
        theta_frac * max_feasible_fraction(K, base.data(), dir_par.data(), lb, ub, nullptr),
        tr_fraction(K, cand_shat.data(), dir.data(), r));
    const double t = line_min(K, ghat, Bhat, cand_shat.data(), dir.data(), t_max);
    if (t > 0.0) {
      for (int i = 0; i < K; ++i) cand_shat[i] += t * dir[i];
      to_par(cand_shat.data(), cand_s.data());
      consider("reflected");
    }
  }

  // Cauchy point: the global-convergence guarantee when both steps above block
  // early.
  for (int i = 0; i < K; ++i) { cand_shat[i] = 0.0; dir[i] = -ghat[i]; }
  to_par(dir.data(), dir_par.data());
  const double t_max = std::min(
      theta_frac * max_feasible_fraction(K, theta, dir_par.data(), lb, ub, nullptr),
      tr_fraction(K, cand_shat.data(), dir.data(), r));
  const double t = line_min(K, ghat, Bhat, cand_shat.data(), dir.data(), t_max);
  if (t > 0.0) {
    for (int i = 0; i < K; ++i) cand_shat[i] = t * dir[i];
    to_par(cand_shat.data(), cand_s.data());
    consider("gradient");
  }
}

}}  // namespace dmod::trust_internal

#endif  // DMOD_TRUST_SUBPROBLEM_H
