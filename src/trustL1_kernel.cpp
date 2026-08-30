// trust-region optimiser with L1 penalty
//
// Minimises  objfun(theta) + sum_i lambda_i * |theta_i - mu_i|
// (or one-sided  lambda_i * max(0, mu_i - theta_i)) over a subset of
// parameters named in mu, on top of the Moré-Sorensen subproblem in
// trust_subproblem.h. Two L1-specific interventions sit in the inner loop:
//
//  1. L1 active set: a penalised coordinate sitting exactly on its kink
//     theta_i = mu_i is dropped from the reduced subproblem whenever
//     |grad_obj_i| <= lambda_i (two-sided) resp. -grad_obj_i <= lambda_i
//     (one-sided). Below that threshold the L1 force dominates and the
//     coordinate stays pinned.
//
//  2. Kink clamping: after the trust step is added to theta, any penalised
//     coordinate that crossed its kink is snapped back to mu_i, so the next
//     iteration can re-examine the active set. The assignment is verbatim --
//     downstream sparsity tests read an exact equality.
//
// Box bounds are handled by `boundary`, exactly as in trust_kernel.cpp:
// "reflective" applies the Coleman-Li scaling to the coordinates that survive
// the L1 active set, "clip" is the frozen historical scheme. The kink active
// set is orthogonal to that choice and is used by both -- L1 sparsity needs
// coordinates to land exactly on mu, which an interior method cannot deliver.

#include <Rcpp.h>
#include "trust_subproblem.h"
#include "trust_driver.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace Rcpp;
using dmod::trust_internal::affine_scaling;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::model_value;
using dmod::trust_internal::stepback;
using dmod::trust_internal::trust_sub;
using dmod::trust_driver::Blather;
using dmod::trust_driver::Reporter;
using dmod::trust_driver::eval_objfun;
using dmod::trust_driver::fill_bound;
using dmod::trust_driver::fill_parscale;
using dmod::trust_driver::kInf;
using dmod::trust_driver::kStallLimit;
using dmod::trust_driver::push_interior;
using dmod::trust_driver::subproblem_label;

namespace {

// Per-parameter L1 metadata, resolved from the named (mu, lambda) pair.
struct L1Spec {
  std::vector<unsigned char> has;
  std::vector<double>        mu;
  std::vector<double>        lambda;
  bool                       one_sided = false;

  void build(const NumericVector& mu_in, const NumericVector& lambda_in,
             const CharacterVector& parnames, int K, bool one_sided_) {
    has.assign(K, 0);
    mu.assign(K, 0.0);
    lambda.assign(K, 0.0);
    one_sided = one_sided_;
    if (mu_in.size() == 0) return;
    if (!mu_in.hasAttribute("names"))
      stop("trustL1: mu must be a named numeric vector");
    if (lambda_in.size() != mu_in.size())
      stop("trustL1: lambda must have the same length as mu");
    CharacterVector mnames = mu_in.names();
    for (int j = 0; j < mu_in.size(); ++j) {
      std::string nm = as<std::string>(mnames[j]);
      for (int i = 0; i < K; ++i)
        if (as<std::string>(parnames[i]) == nm) {
          has[i] = 1; mu[i] = mu_in[j]; lambda[i] = lambda_in[j];
          break;
        }
    }
  }

  double value(const std::vector<double>& th) const {
    double s = 0.0;
    for (std::size_t i = 0; i < th.size(); ++i) {
      if (!has[i]) continue;
      const double d = th[i] - mu[i];
      if (one_sided) { if (d < 0.0) s += lambda[i] * (-d); }
      else           { s += lambda[i] * std::fabs(d); }
    }
    return s;
  }
  // Subgradient of the penalty. At the kink the coordinate is pinned by the
  // active set and dropped, so the zero returned there is never used.
  double grad(int i, double th_i) const {
    if (!has[i]) return 0.0;
    const double d = th_i - mu[i];
    if (one_sided) return (d < 0.0) ? -lambda[i] : 0.0;
    if (d > 0.0) return  lambda[i];
    if (d < 0.0) return -lambda[i];
    return 0.0;
  }
  bool pinned(int i, double th_i, double grad_obj_i) const {
    if (!has[i] || th_i != mu[i]) return false;
    if (one_sided) return (-grad_obj_i) <= lambda[i];
    return std::fabs(grad_obj_i) <= lambda[i];
  }
  // Snap any coordinate that crossed its kink back onto mu, verbatim.
  void clamp(const std::vector<double>& th, std::vector<double>& th_try) const {
    for (std::size_t i = 0; i < th.size(); ++i) {
      if (!has[i]) continue;
      if ((th[i] - mu[i]) * (th_try[i] - mu[i]) < 0.0) th_try[i] = mu[i];
      if (one_sided && th_try[i] < mu[i])              th_try[i] = mu[i];
    }
  }
};

// -------------------------------------------------------------------------
// Coleman-Li interior trust-region-reflective on the L1-active coordinates
//
// One iteration is `propose` (build the reduced subproblem, take a step, clamp
// the kinks) then `accept` (fold in the objective at the trial point). Between
// them sits the single R callback of the whole loop, which is what lets N solves
// share one batched call: see trustL1_lockstep_impl below. Both drivers run this
// same code, so the lockstep is bit-identical to N separate solves by
// construction rather than by agreement.
// -------------------------------------------------------------------------

// Settings shared by every solve in a lockstep round.
struct RefTune {
  double rmax = 0.0, ftol = 0.0, mtol = 0.0, gtol = 0.0, xtol = 0.0;
  double rmin = 0.0, thetamax = 0.0;
  bool   minimize = true, blather_on = false;
};

struct RefState {
  int K = 0;
  L1Spec l1;
  std::vector<double> ps, lbz, ubz;
  std::vector<double> theta, z, grad_obj, H_full;
  double val = 0.0, f_used = 0.0, r = 0.0;
  double opt_measure = kInf, m_value = 0.0, stepnorm = 0.0;
  bool accept = true, converged = false, bail = false, done = false;
  int iter = 0, n_iter = 0, n_fail = 0, n_stall = 0;
  std::string stop_reason = "iterlim";

  std::vector<int>    active;
  std::vector<double> zr, lbr, ubr, gr, Hr, absv, jv, sqrtv, ghat, Bhat;
  std::vector<double> eigvals, eigvecs, shat, shat_step, s_step, shat_real;
  std::vector<double> z_try, theta_try;
  std::vector<unsigned char> at_bound;
  bool is_newton = false, is_hard = false, is_easy = false;
  const char* sb_label = "full";
  Blather trace;
};

// Bounds, scaling and the starting point. `parinit` is read through `get(i)` so
// a matrix row can be fed in without materialising a vector.
template <typename Get>
void ref_init(RefState& s, int K, Get get, const std::vector<double>& pl,
              const std::vector<double>& pu, const std::vector<double>& ps) {
  s.K = K;
  s.ps = ps;
  s.lbz.assign(K, 0.0); s.ubz.assign(K, 0.0); s.z.assign(K, 0.0);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    s.lbz[i] = ps[i] * pl[i];
    s.ubz[i] = ps[i] * pu[i];
    double zi = ps[i] * get(i);
    if (zi > s.ubz[i]) { any_above = true; zi = s.ubz[i]; }
    if (zi < s.lbz[i]) { any_below = true; zi = s.lbz[i]; }
    s.z[i] = zi;
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");
  push_interior(K, s.z, s.lbz, s.ubz);
  s.theta.assign(K, 0.0);
  for (int i = 0; i < K; ++i) s.theta[i] = s.z[i] / ps[i];
  s.at_bound.assign(K, 0);
  s.z_try.assign(K, 0.0);
  s.theta_try.assign(K, 0.0);
}

// Seed the iteration from the objective at parinit.
void ref_seed(RefState& s, double val_obj, const double* grad, const double* H,
              double rinit) {
  const int K = s.K;
  s.grad_obj.assign(grad, grad + K);
  s.H_full.assign(H, H + (std::size_t) K * K);
  s.val = val_obj + s.l1.value(s.theta);
  s.r = rinit;
  s.f_used = kInf;
}

// Build the reduced subproblem and take a step. Returns false when the state
// stopped before proposing anything (gradient convergence).
bool ref_propose(RefState& s, const RefTune& t) {
  const int K = s.K;
  const double sgn = t.minimize ? 1.0 : -1.0;
  int Kred = 0;

  if (s.accept) {
    s.f_used = t.minimize ? s.val : -s.val;

    // Reduced space: drop only the coordinates pinned at their L1 kink. Box
    // bounds are handled by the scaling, not by dropping.
    s.active.clear();
    for (int i = 0; i < K; ++i)
      if (!s.l1.pinned(i, s.theta[i], s.grad_obj[i])) s.active.push_back(i);
    Kred = static_cast<int>(s.active.size());

    s.zr.assign(Kred, 0.0); s.lbr.assign(Kred, 0.0); s.ubr.assign(Kred, 0.0);
    s.gr.assign(Kred, 0.0); s.Hr.assign((std::size_t) Kred * Kred, 0.0);
    for (int ii = 0; ii < Kred; ++ii) {
      const int i = s.active[ii];
      s.zr[ii]  = s.z[i];
      s.lbr[ii] = s.lbz[i];
      s.ubr[ii] = s.ubz[i];
      s.gr[ii]  = sgn * (s.grad_obj[i] + s.l1.grad(i, s.theta[i])) / s.ps[i];
      for (int jj = 0; jj < Kred; ++jj) {
        const int j = s.active[jj];
        s.Hr[ii + (std::size_t) jj * Kred] =
            sgn * s.H_full[i + (std::size_t) j * K] / (s.ps[i] * s.ps[j]);
      }
    }

    s.absv.assign(Kred, 0.0); s.jv.assign(Kred, 0.0); s.sqrtv.assign(Kred, 0.0);
    affine_scaling(Kred, s.zr.data(), s.gr.data(), s.lbr.data(), s.ubr.data(),
                   s.absv.data(), s.jv.data());
    for (int ii = 0; ii < Kred; ++ii) s.sqrtv[ii] = std::sqrt(s.absv[ii]);

    // Pinned coordinates are stationary by the active-set test, so the
    // optimality measure only has to cover the reduced space.
    s.opt_measure = 0.0;
    for (int ii = 0; ii < Kred; ++ii)
      s.opt_measure = std::max(s.opt_measure, std::fabs(s.absv[ii] * s.gr[ii]));
    // Set before the convergence break, which is exactly when it matters.
    const double btol = std::max(t.gtol, 1e-10);
    std::fill(s.at_bound.begin(), s.at_bound.end(), 0);
    for (int ii = 0; ii < Kred; ++ii)
      s.at_bound[s.active[ii]] = (s.jv[ii] > 0.0 &&
                                  std::fabs(s.absv[ii] * s.gr[ii]) <= btol &&
                                  std::fabs(s.gr[ii]) > btol) ? 1 : 0;

    if (s.opt_measure <= t.gtol) {
      s.converged = true; s.stop_reason = "gradient"; s.done = true; return false;
    }

    s.ghat.assign(Kred, 0.0);
    s.Bhat.assign((std::size_t) Kred * Kred, 0.0);
    for (int ii = 0; ii < Kred; ++ii) s.ghat[ii] = s.sqrtv[ii] * s.gr[ii];
    for (int jj = 0; jj < Kred; ++jj)
      for (int ii = 0; ii < Kred; ++ii)
        s.Bhat[ii + (std::size_t) jj * Kred] =
            s.sqrtv[ii] * s.sqrtv[jj] * s.Hr[ii + (std::size_t) jj * Kred];
    for (int ii = 0; ii < Kred; ++ii)
      s.Bhat[ii + (std::size_t) ii * Kred] += std::fabs(s.gr[ii]) * s.jv[ii];

    s.eigvals.assign(Kred, 0.0);
    s.eigvecs.assign((std::size_t) Kred * Kred, 0.0);
    if (Kred > 0)
      eigen_sym_local(s.Bhat.data(), Kred, s.eigvals.data(), s.eigvecs.data());
  }

  Kred = static_cast<int>(s.active.size());

  if (t.blather_on) {
    s.trace.argpath.insert(s.trace.argpath.end(), s.theta.begin(), s.theta.end());
    s.trace.r.push_back(s.r);
    s.trace.valpath.push_back(s.val);
  }

  s.is_newton = false; s.is_hard = false; s.is_easy = false;
  s.sb_label = "full";
  s.m_value = 0.0;
  s.shat.assign(Kred, 0.0);
  s.shat_step.assign(Kred, 0.0);
  s.s_step.assign(Kred, 0.0);
  if (Kred > 0) {
    double pred_unused = 0.0;
    trust_sub(Kred, s.ghat.data(), s.eigvals.data(), s.eigvecs.data(), s.r,
              s.shat.data(), &pred_unused, &s.is_newton, &s.is_hard, &s.is_easy);
    const double theta_frac =
        std::min(std::max(t.thetamax, 1.0 - s.opt_measure), 1.0 - 1e-12);
    stepback(Kred, s.zr.data(), s.lbr.data(), s.ubr.data(), s.sqrtv.data(),
             s.ghat.data(), s.Bhat.data(), s.shat.data(), s.r, theta_frac,
             s.shat_step.data(), s.s_step.data(), &s.m_value, &s.sb_label);
  }

  s.z_try = s.z;
  for (int ii = 0; ii < Kred; ++ii) s.z_try[s.active[ii]] += s.s_step[ii];
  // See trust_kernel.cpp. mu is validated strictly inside the box, so this
  // never disturbs a pinned kink.
  push_interior(K, s.z_try, s.lbz, s.ubz);
  for (int i = 0; i < K; ++i) s.theta_try[i] = s.z_try[i] / s.ps[i];
  s.l1.clamp(s.theta, s.theta_try);
  for (int i = 0; i < K; ++i) s.z_try[i] = s.ps[i] * s.theta_try[i];

  // Rescore the model at the step actually taken -- the kink clamp shortens
  // individual coordinates after the stepback has chosen a candidate.
  s.shat_real.assign(Kred, 0.0);
  bool rescore = true;
  for (int ii = 0; ii < Kred; ++ii) {
    if (!(s.sqrtv[ii] > 0.0)) { rescore = false; break; }
    s.shat_real[ii] = (s.z_try[s.active[ii]] - s.zr[ii]) / s.sqrtv[ii];
  }
  if (rescore && Kred > 0)
    s.m_value = model_value(Kred, s.ghat.data(), s.Bhat.data(), s.shat_real.data());
  const std::vector<double>& shat_used = rescore ? s.shat_real : s.shat_step;

  s.stepnorm = 0.0;
  for (int ii = 0; ii < Kred; ++ii) s.stepnorm += shat_used[ii] * shat_used[ii];
  s.stepnorm = std::sqrt(s.stepnorm);
  return true;
}

// Fold in the objective at theta_try. Returns false when the loop must stop.
bool ref_accept(RefState& s, const RefTune& t, bool eval_ok, double val_obj_try,
                const double* grad_try, const double* H_try_colmajor) {
  const int K = s.K;
  const double val_try = eval_ok ? val_obj_try + s.l1.value(s.theta_try) : kInf;

  const double pred_pos  = -s.m_value;
  const double ftry_used = t.minimize ? val_try : -val_try;
  const double dval      = std::fabs(ftry_used - s.f_used);
  const double rho = (eval_ok && pred_pos > 0.0)
                       ? (s.f_used - ftry_used) / pred_pos : -kInf;

  if (!eval_ok) {
    s.n_fail++;
    s.accept = false;
    s.r *= 0.25;
    if (s.n_fail >= 3) { s.bail = true; s.stop_reason = "objfun"; }
  } else {
    s.n_fail = 0;
    if (rho < 0.25) {
      s.accept = false;
      s.r = std::min(0.25 * s.r, 0.25 * s.stepnorm);
    } else {
      s.accept = true;
      if (rho > 0.75 && s.stepnorm >= 0.9 * s.r) s.r = std::min(2.0 * s.r, t.rmax);
    }
  }

  // Count rejected steps that also left the objective flat.
  if (s.accept || !eval_ok || dval >= t.ftol) s.n_stall = 0;
  else                                        s.n_stall++;

  if (s.accept && eval_ok) {
    s.theta = s.theta_try;
    s.z     = s.z_try;
    s.val   = val_try;
    s.grad_obj.assign(grad_try, grad_try + K);
    s.H_full.assign(H_try_colmajor, H_try_colmajor + (std::size_t) K * K);
  }

  if (t.blather_on) {
    s.trace.argtry.insert(s.trace.argtry.end(), s.theta_try.begin(), s.theta_try.end());
    s.trace.valtry.push_back(val_try);
    s.trace.accept.push_back(s.accept ? 1 : 0);
    s.trace.preddiff.push_back(s.m_value);
    s.trace.stepnorm.push_back(s.stepnorm);
    s.trace.rho.push_back(rho);
    s.trace.steptype.push_back(subproblem_label(s.is_newton, s.is_hard, s.is_easy));
    s.trace.stepback.push_back(s.sb_label);
  }
  s.n_iter = s.iter;

  if (s.bail) { s.done = true; return false; }
  if (s.accept && eval_ok) {
    if (dval < t.ftol)                { s.converged = true; s.stop_reason = "fvalue";   s.done = true; return false; }
    if (std::fabs(s.m_value) < t.mtol){ s.converged = true; s.stop_reason = "preddiff"; s.done = true; return false; }
    if (t.xtol > 0.0 && s.stepnorm < t.xtol)
                                      { s.converged = true; s.stop_reason = "step";     s.done = true; return false; }
  }
  if (s.r < t.rmin) { s.stop_reason = "radius"; s.done = true; return false; }
  if (s.n_stall >= kStallLimit) {
    s.converged = true; s.stop_reason = "stagnation"; s.done = true; return false;
  }
  return true;
}

// Assemble one solve's return value.
List ref_result(RefState& s, const RefTune& t, const CharacterVector& parnames) {
  const int K = s.K;
  NumericVector arg_out(s.theta.begin(), s.theta.end());
  arg_out.names() = parnames;
  // The combined gradient, so the result is self-consistent with `value`.
  NumericVector grad_out(K);
  for (int i = 0; i < K; ++i)
    grad_out[i] = s.grad_obj[i] + s.l1.grad(i, s.theta[i]);
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) Hess_out(i, j) = s.H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (s.at_bound[i] != 0);
  at_bound_out.names() = parnames;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = s.val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = s.n_iter,
      Named("converged")  = s.converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = s.stop_reason);

  if (t.blather_on) s.trace.attach(result, s.n_iter, K, parnames, t.minimize);
  return result;
}

List trustL1_reflective(Function objfun, NumericVector parinit,
                        const L1Spec& l1,
                        double rinit, double rmax,
                        Nullable<NumericVector> parscale,
                        int iterlim,
                        double ftol, double mtol,
                        double gtol, double xtol,
                        double rmin, double thetamax,
                        bool minimize, bool blather_on,
                        Nullable<NumericVector> parupper,
                        Nullable<NumericVector> parlower,
                        bool printIter, Nullable<CharacterVector> traceFile) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  for (int i = 0; i < K; ++i) {
    if (!(pl[i] < pu[i]))
      stop("trustL1: parlower must lie strictly below parupper");
    // The kink is a target the iterate must be able to sit on exactly, so it
    // has to be reachable from inside the box.
    if (l1.has[i] && !(pl[i] < l1.mu[i] && l1.mu[i] < pu[i]))
      stop("trustL1: mu must lie strictly inside [parlower, parupper]");
  }
  fill_parscale(parscale, ps, K, "trustL1");

  RefState s;
  s.l1 = l1;
  ref_init(s, K, [&](int i) { return parinit[i]; }, pl, pu, ps);

  RefTune t;
  t.rmax = rmax; t.ftol = ftol; t.mtol = mtol; t.gtol = gtol; t.xtol = xtol;
  t.rmin = rmin; t.thetamax = thetamax;
  t.minimize = minimize; t.blather_on = blather_on;

  Reporter report;
  report.init(printIter, iterlim, traceFile, K, parnames);

  NumericVector x_named(s.theta.begin(), s.theta.end());
  x_named.names() = parnames;
  List out_init;
  if (!eval_objfun(objfun, x_named, out_init))
    stop("parinit not feasible: objfun failed");
  double val_obj = as<double>(out_init["value"]);
  if (!std::isfinite(val_obj)) stop("parinit not feasible: value is not finite");
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);
  std::vector<double> H0((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) H0[i + (std::size_t) j * K] = Hmat0(i, j);
  ref_seed(s, val_obj, grad0.begin(), H0.data(), rinit);

  int neval = 1;
  report(neval, s.val, x_named, /*head=*/true);

  std::vector<double> H_try((std::size_t) K * K);
  for (int iter = 1; iter <= iterlim; ++iter) {
    R_CheckUserInterrupt();
    s.iter = iter;
    if (!ref_propose(s, t)) break;

    NumericVector x_try(s.theta_try.begin(), s.theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try);
    double val_obj_try = kInf, val_try = kInf;
    NumericVector grad_try;
    if (eval_ok) {
      val_obj_try = as<double>(out_try["value"]);
      if (!std::isfinite(val_obj_try)) eval_ok = false;
      else {
        grad_try = as<NumericVector>(out_try["gradient"]);
        NumericMatrix Htry_mat = as<NumericMatrix>(out_try["hessian"]);
        for (int j = 0; j < K; ++j)
          for (int i = 0; i < K; ++i)
            H_try[i + (std::size_t) j * K] = Htry_mat(i, j);
        val_try = val_obj_try + l1.value(s.theta_try);
      }
    }
    neval++;
    report(neval, val_try, x_try, /*head=*/false);

    if (!ref_accept(s, t, eval_ok, val_obj_try,
                    eval_ok ? grad_try.begin() : nullptr,
                    eval_ok ? H_try.data() : nullptr)) break;
  }

  if (s.stop_reason == "objfun")
    Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
  else if (s.stop_reason == "radius")
    Rf_warning("Trust radius fell below rmin. Fit is not converged.");
  else if (!s.converged && s.n_iter >= iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  return ref_result(s, t, parnames);
}

// -------------------------------------------------------------------------
// N reflective solves in lock-step
//
// Each subject runs its own trustL1 through the very same ref_propose /
// ref_accept as a single solve; the only difference is that one round collects
// the trial points of every subject that has not stopped and sends them out as
// ONE R call. That call is the ODE solve, and the condition axis inside it is
// what cppDE parallelises. Converged subjects drop out of the round, so the
// round count is the maximum over subjects rather than their sum.
// -------------------------------------------------------------------------
List trustL1_lockstep(Function objfun_many, NumericMatrix parinit,
                      NumericMatrix mu, NumericMatrix lambda,
                      bool one_sided, CharacterVector parnames,
                      double rinit, const RefTune& t,
                      Nullable<NumericVector> parscale, int iterlim,
                      Nullable<NumericVector> parupper,
                      Nullable<NumericVector> parlower) {

  const int N = parinit.nrow(), K = parinit.ncol();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  for (int i = 0; i < K; ++i)
    if (!(pl[i] < pu[i]))
      stop("trustL1: parlower must lie strictly below parupper");
  fill_parscale(parscale, ps, K, "trustL1");

  std::vector<RefState> st(N);
  for (int n = 0; n < N; ++n) {
    RefState& s = st[n];
    s.l1.has.assign(K, 0); s.l1.mu.assign(K, 0.0); s.l1.lambda.assign(K, 0.0);
    s.l1.one_sided = one_sided;
    for (int i = 0; i < K; ++i) {
      if (ISNAN(mu(n, i))) continue;   // NA or NaN: no penalty on this coordinate
      s.l1.has[i] = 1; s.l1.mu[i] = mu(n, i); s.l1.lambda[i] = lambda(n, i);
      if (!(pl[i] < mu(n, i) && mu(n, i) < pu[i]))
        stop("trustL1: mu must lie strictly inside [parlower, parupper]");
    }
    ref_init(s, K, [&](int i) { return parinit(n, i); }, pl, pu, ps);
  }

  // One R call for a set of subjects at their current trial points. Returns a
  // list of objlists, NULL where the evaluation failed.
  auto run_round = [&](const std::vector<int>& idx,
                       const std::vector<const double*>& pts) {
    const int M = static_cast<int>(idx.size());
    List parsL(M);
    IntegerVector which(M);
    for (int j = 0; j < M; ++j) {
      NumericVector v(pts[j], pts[j] + K);
      v.names() = parnames;
      parsL[j]  = v;
      which[j]  = idx[j] + 1;
    }
    return as<List>(objfun_many(Rcpp::Named("parsList") = parsL,
                                Rcpp::Named("which")    = which));
  };

  // Read one objlist. Mirrors the single-solve reader, including the "value not
  // finite counts as a failed evaluation" rule.
  auto read_one = [&](SEXP el, double* val_obj, std::vector<double>& grad,
                      std::vector<double>& H) {
    if (Rf_isNull(el)) return false;
    List o(el);
    double v = as<double>(o["value"]);
    if (!std::isfinite(v)) return false;
    NumericVector g = as<NumericVector>(o["gradient"]);
    NumericMatrix Hm = as<NumericMatrix>(o["hessian"]);
    if (g.size() != K || Hm.nrow() != K || Hm.ncol() != K)
      stop("trustL1 lockstep: objfun returned a gradient or hessian of the wrong size");
    *val_obj = v;
    grad.assign(g.begin(), g.end());
    H.assign((std::size_t) K * K, 0.0);
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) H[i + (std::size_t) j * K] = Hm(i, j);
    return true;
  };

  {
    std::vector<int> idx(N);
    std::vector<const double*> pts(N);
    for (int n = 0; n < N; ++n) { idx[n] = n; pts[n] = st[n].theta.data(); }
    List got = run_round(idx, pts);
    std::vector<double> g, H;
    for (int n = 0; n < N; ++n) {
      double v = 0.0;
      if (!read_one(got[n], &v, g, H))
        stop("parinit not feasible: objfun failed");
      ref_seed(st[n], v, g.data(), H.data(), rinit);
    }
  }

  for (int round = 1; round <= iterlim; ++round) {
    R_CheckUserInterrupt();
    std::vector<int> idx;
    std::vector<const double*> pts;
    for (int n = 0; n < N; ++n) {
      if (st[n].done) continue;
      st[n].iter = round;
      if (!ref_propose(st[n], t)) continue;   // stopped on the gradient test
      idx.push_back(n);
      pts.push_back(st[n].theta_try.data());
    }
    if (idx.empty()) break;

    List got = run_round(idx, pts);
    std::vector<double> g, H;
    for (std::size_t j = 0; j < idx.size(); ++j) {
      const int n = idx[j];
      double v = 0.0;
      const bool ok = read_one(got[j], &v, g, H);
      ref_accept(st[n], t, ok, v, ok ? g.data() : nullptr, ok ? H.data() : nullptr);
    }
  }

  List out(N);
  for (int n = 0; n < N; ++n) {
    RefState& s = st[n];
    if (s.stop_reason == "objfun")
      Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
    else if (s.stop_reason == "radius")
      Rf_warning("Trust radius fell below rmin. Fit is not converged.");
    else if (!s.converged && s.n_iter >= iterlim)
      Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");
    out[n] = ref_result(s, t, parnames);
  }
  return out;
}

// -------------------------------------------------------------------------
// Legacy: active-set reduction plus componentwise clipping
// -------------------------------------------------------------------------
List trustL1_clip(Function objfun, NumericVector parinit,
                  const L1Spec& l1,
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
  fill_parscale(parscale, ps, K, "trustL1");

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
  List out_init = as<List>(objfun(x_named));
  double val_obj = as<double>(out_init["value"]);
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);
  if (!std::isfinite(val_obj)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_obj(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  double val = val_obj + l1.value(theta);
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
        double gi = grad_obj[i];
        bool drop_bound;
        if (minimize) drop_bound = (at_upper[i] && gi < 0.0) || (at_lower[i] && gi > 0.0);
        else          drop_bound = (at_upper[i] && gi > 0.0) || (at_lower[i] && gi < 0.0);
        if (!drop_bound && !l1.pinned(i, theta[i], gi)) active.push_back(i);
      }
      int Kred = static_cast<int>(active.size());
      g_red.assign(Kred, 0.0);
      H_red.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        int i = active[ii];
        g_red[ii] = grad_obj[i] + l1.grad(i, theta[i]);
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
    l1.clamp(theta, theta_try);

    for (int i = 0; i < K; ++i) {
      at_upper[i] = !(theta_try[i] < pu[i]) ? 1 : 0;
      at_lower[i] = !(theta_try[i] > pl[i]) ? 1 : 0;
      if (at_upper[i]) theta_try[i] = pu[i];
      if (at_lower[i]) theta_try[i] = pl[i];
    }

    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try);
    double val_obj_try = kInf, val_try = kInf;
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_obj_try = as<double>(out_try["value"]);
      grad_try    = as<NumericVector>(out_try["gradient"]);
      Htry_mat    = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_obj_try)) eval_ok = false;
      else val_try = val_obj_try + l1.value(theta_try);
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
        Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
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
      theta = theta_try;
      val   = val_try;
      grad_obj.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
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
  arg_out.names() = parnames;
  NumericVector grad_out(K);
  for (int i = 0; i < K; ++i) grad_out[i] = grad_obj[i] + l1.grad(i, theta[i]);
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
List trustL1_impl(Function objfun,
                  NumericVector parinit,
                  NumericVector mu,
                  NumericVector lambda,
                  bool   one_sided,
                  double rinit,
                  double rmax,
                  Nullable<NumericVector> parscale = R_NilValue,
                  int    iterlim   = 100,
                  double ftol      = 1e-6,
                  double mtol      = 1e-6,
                  double gtol      = 1e-6,
                  double xtol      = 0.0,
                  double rmin      = 0.0,
                  double thetamax  = 0.99995,
                  std::string boundary = "reflective",
                  bool   minimize  = true,
                  bool   blather   = false,
                  Nullable<NumericVector>  parupper  = R_NilValue,
                  Nullable<NumericVector>  parlower  = R_NilValue,
                  bool   printIter = false,
                  Nullable<CharacterVector> traceFile = R_NilValue) {

  const int K = parinit.size();
  if (K == 0) stop("trustL1: parinit must be non-empty");
  if (!parinit.hasAttribute("names"))
    stop("trustL1: parinit must be a named numeric vector");
  for (int i = 0; i < K; ++i)
    if (!std::isfinite(parinit[i])) stop("trustL1: parinit not all finite");

  L1Spec l1;
  l1.build(mu, lambda, parinit.names(), K, one_sided);

  if (boundary == "clip")
    return trustL1_clip(objfun, parinit, l1, rinit, rmax, parscale, iterlim,
                        ftol, mtol, minimize, blather,
                        parupper, parlower, printIter, traceFile);
  if (boundary != "reflective")
    stop("trustL1: boundary must be one of \"reflective\", \"clip\"");

  return trustL1_reflective(objfun, parinit, l1, rinit, rmax, parscale, iterlim,
                            ftol, mtol, gtol, xtol, rmin, thetamax,
                            minimize, blather, parupper, parlower,
                            printIter, traceFile);
}


// N reflective trustL1 solves stepping in unison. `objfun_many(parsList, which)`
// returns one objlist per requested subject, NULL where the evaluation failed.
// `mu` and `lambda` are N x K; NA in `mu` marks an unpenalised coordinate.
// [[Rcpp::export]]
List trustL1_lockstep_impl(Function objfun_many,
                           NumericMatrix parinit,
                           NumericMatrix mu,
                           NumericMatrix lambda,
                           bool   one_sided,
                           double rinit,
                           double rmax,
                           Nullable<NumericVector> parscale = R_NilValue,
                           int    iterlim   = 100,
                           double ftol      = 1e-6,
                           double mtol      = 1e-6,
                           double gtol      = 1e-6,
                           double xtol      = 0.0,
                           double rmin      = 0.0,
                           double thetamax  = 0.99995,
                           bool   minimize  = true,
                           bool   blather   = false,
                           Nullable<NumericVector> parupper = R_NilValue,
                           Nullable<NumericVector> parlower = R_NilValue) {

  const int N = parinit.nrow(), K = parinit.ncol();
  if (N == 0 || K == 0) stop("trustL1: parinit must be a non-empty N x K matrix");
  if (mu.nrow() != N || mu.ncol() != K || lambda.nrow() != N || lambda.ncol() != K)
    stop("trustL1: mu and lambda must have the same shape as parinit");
  List dn = parinit.attr("dimnames");
  if (dn.size() < 2 || Rf_isNull(dn[1]))
    stop("trustL1: parinit must carry column names");
  CharacterVector parnames = dn[1];
  for (int n = 0; n < N; ++n)
    for (int i = 0; i < K; ++i)
      if (!std::isfinite(parinit(n, i))) stop("trustL1: parinit not all finite");

  RefTune t;
  t.rmax = rmax; t.ftol = ftol; t.mtol = mtol; t.gtol = gtol; t.xtol = xtol;
  t.rmin = rmin; t.thetamax = thetamax;
  t.minimize = minimize; t.blather_on = blather;

  return trustL1_lockstep(objfun_many, parinit, mu, lambda, one_sided, parnames,
                          rinit, t, parscale, iterlim, parupper, parlower);
}
