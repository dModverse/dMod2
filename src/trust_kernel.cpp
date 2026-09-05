// trust-region optimiser.
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

// Interchangeable Hessian source (Stage 1 of dev/reverseAD-plan.md). The
// quasi-Newton approximation lives in the working Hessian `H_full`, so the
// scaling, eigendecomposition and subproblem downstream are untouched.
enum HessianMethod { HM_GN = 0, HM_BFGS = 1, HM_SR1 = 2, HM_HYBRID = 3 };

// Dense quasi-Newton update of `H` in place. `H` holds the Hessian of `val`; the
// run minimises phi = sgn*val, so the update acts on B = sgn*H with
// y = sgn*(grad_new-grad_old) and writes back H = sgn*B. `s` is the accepted
// step (x frame). BFGS is Powell-damped (stays positive definite); SR1 may be
// indefinite, which the subproblem solves exactly.
void qn_update(int kind, int K, std::vector<double>& H,
               const std::vector<double>& s, const std::vector<double>& y_val,
               double sgn) {
  const std::size_t KK = (std::size_t) K * K;
  std::vector<double> B(KK), y(K);
  for (std::size_t t = 0; t < KK; ++t) B[t] = sgn * H[t];
  for (int i = 0; i < K; ++i) y[i] = sgn * y_val[i];

  std::vector<double> Bs(K, 0.0);
  for (int j = 0; j < K; ++j) {
    const double sj = s[j];
    if (sj != 0.0)
      for (int i = 0; i < K; ++i) Bs[i] += B[i + (std::size_t) j * K] * sj;
  }
  double sy = 0.0, sBs = 0.0;
  for (int i = 0; i < K; ++i) { sy += s[i] * y[i]; sBs += s[i] * Bs[i]; }

  if (kind == HM_BFGS) {
    if (!(sBs > 0.0)) return;                 // no curvature reference; keep B
    double theta = 1.0;
    if (sy < 0.2 * sBs) theta = (0.8 * sBs) / (sBs - sy);
    std::vector<double> rv(K);
    for (int i = 0; i < K; ++i) rv[i] = theta * y[i] + (1.0 - theta) * Bs[i];
    double sr = 0.0;
    for (int i = 0; i < K; ++i) sr += s[i] * rv[i];
    if (!(sr > 0.0)) return;
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        B[i + (std::size_t) j * K] += rv[i] * rv[j] / sr - Bs[i] * Bs[j] / sBs;
  } else {  // HM_SR1
    std::vector<double> w(K);
    for (int i = 0; i < K; ++i) w[i] = y[i] - Bs[i];
    double ws = 0.0, wn = 0.0, sn = 0.0;
    for (int i = 0; i < K; ++i) { ws += w[i] * s[i]; wn += w[i] * w[i]; sn += s[i] * s[i]; }
    if (std::fabs(ws) <= 1e-8 * std::sqrt(sn * wn) || !(wn > 0.0)) return;
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        B[i + (std::size_t) j * K] += w[i] * w[j] / ws;
  }
  for (std::size_t t = 0; t < KK; ++t) H[t] = sgn * B[t];
}

inline const char* hessian_source_label(bool qn_active, int qn_kind) {
  if (!qn_active)        return "gn";
  if (qn_kind == HM_SR1) return "sr1";
  return "bfgs";
}

// -------------------------------------------------------------------------
// Coleman-Li interior trust-region-reflective
// -------------------------------------------------------------------------
List trust_reflective(Function objfun, NumericVector parinit,
                      double rinit, double rmax,
                      Nullable<NumericVector> parscale,
                      int iterlim,
                      double ftol, double mtol,
                      double gtol, double xtol,
                      double rmin, double thetamax,
                      int hessianMethod,
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

  List out_init;
  if (!eval_objfun(objfun, x_named, out_init))
    stop("parinit not feasible: objfun failed");
  double val = as<double>(out_init["value"]);
  if (!std::isfinite(val)) stop("parinit not feasible: value is not finite");
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);

  std::vector<double> grad_full(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  int neval = 1;
  report(neval, val, x_named, /*head=*/true);

  Blather trace;
  std::vector<double> absv(K), jv(K), sqrtv(K), g_z(K), ghat(K);
  std::vector<double> Bhat((std::size_t) K * K), eigvals(K), eigvecs((std::size_t) K * K);
  std::vector<double> shat(K), shat_step(K), s_step(K), z_try(K);
  std::vector<unsigned char> at_bound(K, 0);

  double r = rinit;
  double f_used = minimize ? val : -val;
  double opt_measure = kInf;
  bool   accept = true, converged = false, bail = false;
  int    n_iter = 0, n_fail = 0, n_stall = 0;
  std::string stop_reason = "iterlim";

  // Hessian source: gn passes the objective's J^T J through; bfgs/sr1 maintain
  // their own update seeded from it; hybrid runs gn until the first stagnation,
  // then switches to bfgs once. Hmat0 already seeds H_full.
  const bool is_hybrid = (hessianMethod == HM_HYBRID);
  bool qn_active = (hessianMethod == HM_BFGS || hessianMethod == HM_SR1);
  int  qn_kind  = (hessianMethod == HM_SR1) ? HM_SR1 : HM_BFGS;
  bool switched = false;
  int  qn_neval = 0;
  const double qn_sgn = minimize ? 1.0 : -1.0;
  std::vector<double> s_qn(K), y_qn(K);

  for (int iter = 1; iter <= iterlim; ++iter) {
    R_CheckUserInterrupt();

    if (accept) {
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
      // stationary -- its own gradient is still well away from zero. Set before
      // the convergence break, which is exactly when it matters.
      const double btol = std::max(gtol, 1e-10);
      for (int i = 0; i < K; ++i)
        at_bound[i] = (jv[i] > 0.0 && std::fabs(absv[i] * g_z[i]) <= btol &&
                       std::fabs(g_z[i]) > btol) ? 1 : 0;

      if (opt_measure <= gtol) { converged = true; stop_reason = "gradient"; break; }

      for (int i = 0; i < K; ++i) ghat[i] = sqrtv[i] * g_z[i];
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          Bhat[i + (std::size_t) j * K] =
              sqrtv[i] * sqrtv[j] * sgn * H_full[i + (std::size_t) j * K] / (ps[i] * ps[j]);
      for (int i = 0; i < K; ++i)
        Bhat[i + (std::size_t) i * K] += std::fabs(g_z[i]) * jv[i];

      eigen_sym_local(Bhat.data(), K, eigvals.data(), eigvecs.data());
    }

    if (blather_on) {
      for (int i = 0; i < K; ++i) trace.argpath.push_back(z[i] / ps[i]);
      trace.r.push_back(r);
      trace.valpath.push_back(val);
    }

    bool is_newton = false, is_hard = false, is_easy = false;
    double pred_unused = 0.0;
    trust_sub(K, ghat.data(), eigvals.data(), eigvecs.data(), r,
              shat.data(), &pred_unused, &is_newton, &is_hard, &is_easy);

    // Let the iterate approach a face as the optimality measure falls, but
    // never reach it -- |v| = 0 would freeze that coordinate for good.
    const double theta_frac =
        std::min(std::max(thetamax, 1.0 - opt_measure), 1.0 - 1e-12);

    const char* sb_label = "full";
    double m_value = 0.0;
    stepback(K, z.data(), lbz.data(), ubz.data(), sqrtv.data(),
             ghat.data(), Bhat.data(), shat.data(), r, theta_frac,
             shat_step.data(), s_step.data(), &m_value, &sb_label);

    double stepnorm = 0.0;
    for (int i = 0; i < K; ++i) stepnorm += shat_step[i] * shat_step[i];
    stepnorm = std::sqrt(stepnorm);

    for (int i = 0; i < K; ++i) z_try[i] = z[i] + s_step[i];
    // theta_frac keeps the step interior in exact arithmetic, but once |v| is
    // small the remaining gap underflows and z_try rounds onto the bound,
    // setting |v| = 0 and freezing the coordinate.
    push_interior(K, z_try, lbz, ubz);

    NumericVector x_try(K);
    x_try.names() = parnames;
    for (int i = 0; i < K; ++i) x_try[i] = z_try[i] / ps[i];

    // The quasi-Newton phase maintains its own Hessian, so it asks the
    // objective not to build J^T J; only the gn pass-through reads it back.
    const bool want_h = !qn_active;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try, want_h);
    double val_try = kInf;
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_try  = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      if (want_h) Htry_mat = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_try)) eval_ok = false;
    }
    neval++;
    if (qn_active) qn_neval++;
    report(neval, val_try, x_try, /*head=*/false);

    const double pred_pos  = -m_value;
    const double ftry_used = minimize ? val_try : -val_try;
    const double dval      = std::fabs(ftry_used - f_used);
    const double rho = (eval_ok && pred_pos > 0.0)
                         ? (f_used - ftry_used) / pred_pos : -kInf;

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

    if (accept && eval_ok) {
      // (s, y) from the accepted step, before z/grad_full are overwritten.
      if (qn_active)
        for (int i = 0; i < K; ++i) {
          s_qn[i] = (z_try[i] - z[i]) / ps[i];
          y_qn[i] = grad_try[i] - grad_full[i];
        }
      z = z_try;
      val = val_try;
      grad_full.assign(grad_try.begin(), grad_try.end());
      if (qn_active)
        qn_update(qn_kind, K, H_full, s_qn, y_qn, qn_sgn);
      else
        for (int j = 0; j < K; ++j)
          for (int i = 0; i < K; ++i)
            H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
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
      trace.hsource.push_back(hessian_source_label(qn_active, qn_kind));
    }
    n_iter = iter;

    if (bail) break;
    if (accept && eval_ok) {
      if (dval < ftol) { converged = true; stop_reason = "fvalue"; break; }
      if (std::fabs(m_value) < mtol)             { converged = true; stop_reason = "preddiff";  break; }
      if (xtol > 0.0 && stepnorm < xtol)         { converged = true; stop_reason = "step";   break; }
    }
    if (r < rmin) { stop_reason = "radius"; break; }
    if (n_stall >= kStallLimit) {
      // Hybrid: first stagnation switches gn -> bfgs and restarts the trust
      // region from the last accepted J^T J seed; the second stagnation ends it.
      if (is_hybrid && !switched) {
        switched = true; qn_active = true; qn_kind = HM_BFGS;
        n_stall = 0; r = rinit;
      } else {
        converged = true; stop_reason = "stagnation"; break;
      }
    }
  }

  if (stop_reason == "objfun")
    Rf_warning("trust: objfun evaluation failed 3 times in a row");
  else if (stop_reason == "radius")
    Rf_warning("Trust radius fell below rmin. Fit is not converged.");
  else if (!converged && n_iter >= iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  NumericVector arg_out(K), grad_out(grad_full.begin(), grad_full.end());
  for (int i = 0; i < K; ++i) arg_out[i] = z[i] / ps[i];
  arg_out.names()  = parnames;
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) Hess_out(i, j) = H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (at_bound[i] != 0);
  at_bound_out.names() = parnames;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = n_iter,
      Named("neval")      = neval,
      Named("qnEval")     = qn_neval,
      Named("converged")  = converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = stop_reason);

  if (blather_on) trace.attach(result, n_iter, K, parnames, minimize);
  return result;
}

// -------------------------------------------------------------------------
// Legacy: active-set reduction plus componentwise clipping
// -------------------------------------------------------------------------
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
  List out_init = as<List>(objfun(x_named));
  double val = as<double>(out_init["value"]);
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);
  if (!std::isfinite(val)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_full(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

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
    bool eval_ok = eval_objfun(objfun, x_try, out_try);
    double val_try = kInf;
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_try  = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      Htry_mat = as<NumericMatrix>(out_try["hessian"]);
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
                Nullable<NumericVector> parscale = R_NilValue,
                int    iterlim   = 100,
                double ftol      = 1e-6,
                double mtol      = 1e-6,
                double gtol      = 1e-6,
                double xtol      = 0.0,
                double rmin      = 0.0,
                double thetamax  = 0.99995,
                std::string boundary = "reflective",
                std::string hessianMethod = "gn",
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

  int hm = HM_GN;
  if      (hessianMethod == "gn")     hm = HM_GN;
  else if (hessianMethod == "bfgs")   hm = HM_BFGS;
  else if (hessianMethod == "sr1")    hm = HM_SR1;
  else if (hessianMethod == "hybrid") hm = HM_HYBRID;
  else stop("trust: hessianMethod must be one of \"gn\", \"bfgs\", \"sr1\", \"hybrid\"");

  if (boundary == "clip") {
    if (hm != HM_GN)
      stop("trust: hessianMethod other than \"gn\" requires boundary=\"reflective\"");
    return trust_clip(objfun, parinit, rinit, rmax, parscale, iterlim,
                      ftol, mtol, minimize, blather,
                      parupper, parlower, printIter, traceFile);
  }
  if (boundary != "reflective")
    stop("trust: boundary must be one of \"reflective\", \"clip\"");

  return trust_reflective(objfun, parinit, rinit, rmax, parscale, iterlim,
                          ftol, mtol, gtol, xtol, rmin, thetamax,
                          hm, minimize, blather, parupper, parlower,
                          printIter, traceFile);
}
