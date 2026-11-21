// The group kernel of the NLME methods: per evaluation of the marginal, for all groups at once, the
// inner modes, the local model of every group at the mode and the integral of every group by an
// Integrator.
//
// Group i has etas eta in R^K and the inner objective
//   Phi_i(eta) = -2 log p(y_i | eta, psi) + eta' P eta + log|Omega|,
// P the prior precision of the search (Omega^-1 for normal effects), without the constant of the
// normal density. The kernel forms, in this order:
//
//   mode       m_i = argmin Phi_i, trust region on the Gauss-Newton curvature with an SR1 correction;
//              the trial points of all groups in one batched prediction per round (predict_cb).
//   local      Phi_i(m_i + d) ~ Phi_i(m_i) + g_i' d + d' H_i d / 2, g_i the gradient the search
//   model      leaves, H_i the Gauss-Newton curvature with 2 P. The coordinates split into normal
//              ones N and others L (`normal`, all of them unless the caller names a subset).
//   normal     H_NN by the eigendecomposition, eigenvalues floored at eigen_floor_relative
//   block      tr(H_NN) / |N|: log|H_NN| and the inverse H_NN^-1 on the floored eigenvalues.
//   Schur      eliminate_normal(): T - C H_NN^-1 C' and t - C H_NN^-1 g_N for a block C coupled to N.
//              With C = H_LN, T = H_LL, t = g_L it gives S_i and r_i, the curvature and gradient on
//              L once N is integrated out; with C = H_psi,N of the joint Hessian the outer model.
//   integral   -2 log L_i ~ Phi_i(m_i) + term_i, plus constant() over all groups, by an Integrator
//              from the local model alone.
//
// A subject may have fewer coordinates than the problem: `coords` in the subject description names them
// (1-based among the K), the others rest where `pars` holds them; P is then restricted to them.
//
// Interface kernel to integrator. An Integrator reads a LocalModel and returns a GroupIntegral. The
// term is a function of the data parts of the local model, H - 2 P and g - 2 P m, of the mode m and of
// the hyperparameters of the integrator:
//   term     the share of -2 log L_i beyond Phi_i(m_i), without the constant;
//   weight   W_i (K x K) = d term_i / d H_i: the derivative of the term in psi is tr(W_i dH_i/dpsi),
//            dH_i/dpsi the total derivative, the mode moving with psi (the volume callbacks in R);
//   grad     d term_i / d g_i, empty where the term does not depend on g;
//   eta      d term_i / d m_i at fixed data parts, empty likewise;
//   hyper    d term_i / d h for the hyperparameters h of the integrator, in the order it names;
//   mean,    the posterior mean and covariance on the coordinates L the integrator integrates,
//   cov      where it forms them.
// The derivative of -2 log L_i in psi is then that of Phi_i at the mode plus tr(W dH) + grad' dg +
// eta' dm, with dm the mode shift. LaplaceNormal integrates the normal block: term log|H_NN|,
// weight H_NN^-1, constant -|groups| |N| log 2, which is FOCEI. An integrator over L reads S_i, g_L
// and m_L; H includes 2 P on L, which such an integrator, holding a separate prior on L, takes off S_i.

#ifndef DMOD_NLME_KERNEL_H
#define DMOD_NLME_KERNEL_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace dmod {
namespace nlme {

// Tolerance `key` from a control list, `dflt` where the list lacks it.
inline double ctrl_tol(const Rcpp::List& ctrl, const char* key, double dflt) {
  if (ctrl.containsElementNamed(key)) return Rcpp::as<double>(ctrl[key]);
  return dflt;
}

// The radius after a rejected step: a quarter of the radius or of the step, whichever is
// shorter, so that a rejected Newton step well inside the region is not tried again.
inline double shrunk_radius(double r, const double* step, int n) {
  double ss = 0.0;
  for (int k = 0; k < n; ++k) ss += step[k] * step[k];
  return 0.25 * std::min(r, std::sqrt(ss));
}

// Settings of the inner search.
struct InnerControl {
  double rinit = 1.0, rmax = 10.0, steptol = 1e-8, mtol = 1e-12, eigen_floor_relative = 1e-10;
  int iterlim = 30;
};
InnerControl inner_control(const Rcpp::List& control);

// The local model of one group at the mode; matrices K x K column-major, K the coordinates of the
// subject, at the positions `coords` among those of the problem.
struct LocalModel {
  int K = 0, group = 0;
  std::vector<int> coords;
  std::vector<double> eta;               // m_i
  std::vector<double> grad;              // g_i
  std::vector<double> H;                 // H_i
  std::vector<double> eigvals, eigvecs;  // of H_i, ascending, from the last step of the search
  std::vector<int> normal, other;        // N and L, 0-based
  std::vector<double> H_inv;             // H_NN^-1 on the floored eigenvalues, |N| x |N|
  double log_det = 0.0;                  // log|H_NN| on the floored eigenvalues
  std::vector<double> S, r;              // the curvature and gradient on L, empty without L
  double value = 0.0;                    // Phi_i(m_i)
  int iterations = 0;
  bool converged = false;
};

// T (n x n) -= C H_NN^-1 C' and t (n) -= C H_NN^-1 g_N for C (n x |N|, leading dimension n), on the
// floored inverse of `m`; `t` may be null, `g` holds g_N.
void eliminate_normal(int n, const LocalModel& m, const double* C, double* T, double* t, const double* g);

// What an integrator gives for one group; see the header comment.
struct GroupIntegral {
  bool ok = true;
  double term = 0.0;
  std::vector<double> weight, grad, eta, hyper, mean, cov;
};

class Integrator {
 public:
  virtual ~Integrator() = default;
  virtual void integrate(const LocalModel& m, GroupIntegral& out) const = 0;
  virtual double constant(const std::vector<LocalModel>& models) const = 0;
};

// Laplace on the normal block, FOCEI.
class LaplaceNormal : public Integrator {
 public:
  void integrate(const LocalModel& m, GroupIntegral& out) const override;
  double constant(const std::vector<LocalModel>& models) const override;
};

// The inner problem of one evaluation: model_cb and err_cb predict one condition block, predict_cb
// (when not null) all blocks of a round in one call; `pars` the full parameter vector, the etas of
// a group at the rows of `warm`, `meta` the group description of .foceiMeta(), P = Omega_inv (K x K).
struct GroupProblem {
  Rcpp::Function* model_cb = nullptr;
  Rcpp::Function* err_cb = nullptr;
  Rcpp::Function* predict_cb = nullptr;
  Rcpp::Nullable<Rcpp::NumericVector> fixed = R_NilValue;
  Rcpp::List meta;
  std::vector<double> Omega_inv;
  double Omega_log_det = 0.0;
  std::vector<int> normal;  // every coordinate unless normal_given
  bool normal_given = false;
  InnerControl control;
};

// The groups at one evaluation: local models, integrals and the constant.
struct Groups {
  std::vector<LocalModel> models;
  std::vector<GroupIntegral> integrals;
  double constant = 0.0;
};

// The coordinates of a subject among the K of the problem, 0-based.
std::vector<int> subject_coords(const Rcpp::List& meta_i, int K);

// The modes from `warm` (N x K), the local models and their integrals; the modes are written into
// `pars` at the eta positions of the groups.
Groups evaluate_groups(const GroupProblem& prob, Rcpp::NumericVector& pars, const Rcpp::NumericMatrix& warm,
                       const Integrator& integrator);

// The groups as the list the R side reads: eta_modes, gradient, H_GN, H_inv, eigvals, log_det_H,
// value, iterations, converged per group, the integral as term, weight and constant, S and r with L.
Rcpp::List groups_list(const Groups& groups, const Rcpp::List& meta);

}  // namespace nlme
}  // namespace dmod

#endif
