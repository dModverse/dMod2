// dMod FOCEI kernel: the outer problem on the group kernel of nlme_kernel.h.
//
// Math anchors (load-bearing, do not change without re-deriving):
//   - Envelope theorem: at eta = eta*(theta), d/dtheta f(theta, eta*) =
//     partial_theta f, since partial_eta f = 0 at the inner stationary point.
//   - Schur complement (block Hessian of the outer marginal at modes):
//       H_outer = H_thetatheta - H_thetaeta H_etaeta^{-1} H_etatheta
//     where H_etaeta is the per-subject GN H_i, summed over subjects.
//   - Laplace approximation:
//       -2 log p(y | theta) ~= sum_i 2 f_i(eta_i*, theta) + log|H_i / (2 pi)|
//     dMod uses H_GN_i in the log|H| term (not Newton).

#include <Rcpp.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <functional>

#include "nlme_kernel.h"
#include "trust_model.h"
#include "trust_driver.h"

using namespace Rcpp;
using dmod::nlme::ctrl_tol;
using dmod::nlme::shrunk_radius;



// Build 0-based indices into `haystack` for each name in `needles`; throws if
// a needle is missing.
static std::vector<int> name_indices(const CharacterVector& haystack,
                                     const CharacterVector& needles,
                                     const char* ctx) {
  std::vector<int> idx(needles.size(), -1);
  for (int i = 0; i < needles.size(); ++i) {
    const std::string nn = as<std::string>(needles[i]);
    for (int j = 0; j < haystack.size(); ++j) {
      if (as<std::string>(haystack[j]) == nn) { idx[i] = j; break; }
    }
    if (idx[i] < 0) {
      std::string msg = std::string(ctx) + ": name not found: " + nn;
      throw std::runtime_error(msg);
    }
  }
  return idx;
}


// Adds the volume derivatives of `correction_cb`, a list of the gradient named by
// outer parameter and an optional Hessian with outer dimnames, to the outer model;
// the callback gets the outer gradient without them as `data_gradient`. Returns the
// `eta_shift`, d eta*/d theta per subject (K x n_outer), or an empty list.
static List add_correction(Function& correction_cb, const CharacterVector& outer_names,
                           const NumericVector& pars_full, const NumericMatrix& H_full,
                           const List& H_inv_list, const List& weight_list, NumericVector& grad_outer,
                           NumericMatrix& Hess) {
  List corr = as<List>(correction_cb(
      Rcpp::Named("full_pars")     = pars_full,
      Rcpp::Named("joint_hessian") = H_full,
      Rcpp::Named("H_inv_list")    = H_inv_list,
      Rcpp::Named("weight_list")   = weight_list,
      Rcpp::Named("data_gradient") = clone(grad_outer)));
  NumericVector g = as<NumericVector>(corr["gradient"]);
  std::vector<int> at = name_indices(outer_names, CharacterVector(g.names()), "correction/gradient");
  for (int j = 0; j < (int)at.size(); ++j) grad_outer[at[j]] += g[j];
  if (!Rf_isNull(corr["hessian"])) {
    NumericMatrix h = as<NumericMatrix>(corr["hessian"]);
    List hd = h.attr("dimnames");
    std::vector<int> hi = name_indices(outer_names, as<CharacterVector>(hd[0]), "correction/hessian");
    for (int c = 0; c < (int)hi.size(); ++c)
      for (int r = 0; r < (int)hi.size(); ++r) Hess(hi[r], hi[c]) += h(r, c);
  }
  return corr.containsElementNamed("eta_shift") && !Rf_isNull(corr["eta_shift"]) ?
           as<List>(corr["eta_shift"]) : List();
}

// Moves the modes in `eta` (N x K) to first order along the step, eta += S_i step, the columns of
// S_i at the outer positions in its attribute "idx" where it has one; a subject whose move exceeds
// `radius`, the region the inner search starts with, keeps its modes.
static void shift_modes(NumericMatrix& eta, const List& shift, const std::vector<double>& step, double radius) {
  const int N = eta.nrow(), K = eta.ncol(), n = step.size();
  if (shift.size() != N) return;
  std::vector<double> d(K);
  for (int i = 0; i < N; ++i) {
    if (Rf_isNull(shift[i])) continue;
    NumericMatrix S = as<NumericMatrix>(shift[i]);
    SEXP idx = Rf_getAttrib(S, Rf_install("idx"));
    const bool local = !Rf_isNull(idx);
    if (S.nrow() != K || (!local && S.ncol() != n) || (local && Rf_length(idx) != S.ncol())) continue;
    double norm2 = 0.0;
    for (int k = 0; k < K; ++k) {
      d[k] = 0.0;
      for (int j = 0; j < S.ncol(); ++j) d[k] += S(k, j) * step[local ? INTEGER(idx)[j] - 1 : j];
      norm2 += d[k] * d[k];
    }
    if (std::isfinite(norm2) && std::sqrt(norm2) <= radius) for (int k = 0; k < K; ++k) eta(i, k) += d[k];
  }
}

// Outer FOCEI objective: OFV from the per-subject inner values, and one joint_cb
// call at the modes for the outer gradient and the cross block H[outer, eta]
// that the Schur complement needs.
List focei_outer_objfn(Function model_cb,
                            Function err_cb,
                            Function joint_cb,
                            NumericVector outer_pars,
                            NumericMatrix eta_warmstart,
                            List subject_meta,
                            NumericMatrix Omega_inv_mat,
                            double Omega_log_det,
                            Nullable<NumericVector> fixed,
                            List inner_ctrl,
                            std::string correction_mode = "none",
                            Nullable<Function> correction_cb_opt = R_NilValue,
                       Nullable<Function> predict_cb = R_NilValue,
                       bool apply_correction = true,
                       bool newton_value = false) {
  CharacterVector subjects        = subject_meta["subjects"];
  List            eta_names_list  = subject_meta["eta_names"];
  CharacterVector outer_names     = subject_meta["outer_names"];
  IntegerVector   outer_idx_full  = subject_meta["outer_idx_in_full"];
  NumericVector   other_etas_init = subject_meta["other_etas_init"];
  CharacterVector pars_full_names = subject_meta["pars_full_names"];
  const int N = subjects.size();
  const int n_outer = outer_names.size();

  // Build pars_full
  const int n_full = pars_full_names.size();
  NumericVector pars_full(n_full);
  pars_full.attr("names") = pars_full_names;
  for (int i = 0; i < n_full; ++i) pars_full[i] = 0.0;
  {
    CharacterVector init_names = other_etas_init.names();
    std::vector<int> idx = name_indices(pars_full_names, init_names,
                                        "fast_outer/other_etas_init");
    for (int i = 0; i < (int)idx.size(); ++i)
      pars_full[idx[i]] = other_etas_init[i];
  }
  for (int i = 0; i < n_outer; ++i)
    pars_full[outer_idx_full[i] - 1] = outer_pars[i];

  // the group kernel with Laplace on the normal block; the modes go into pars_full
  dmod::nlme::GroupProblem prob;
  prob.model_cb = &model_cb;
  prob.err_cb = &err_cb;
  Function pcb = predict_cb.isNotNull() ? Function(predict_cb) : model_cb;
  if (predict_cb.isNotNull()) prob.predict_cb = &pcb;
  prob.fixed = fixed;
  prob.meta = subject_meta;
  prob.Omega_inv.assign(Omega_inv_mat.begin(), Omega_inv_mat.end());
  prob.Omega_log_det = Omega_log_det;
  prob.control = dmod::nlme::inner_control(inner_ctrl);
  const dmod::nlme::LaplaceNormal laplace;
  const dmod::nlme::Groups groups = dmod::nlme::evaluate_groups(prob, pars_full, eta_warmstart, laplace);
  List inner_res = dmod::nlme::groups_list(groups, subject_meta);

  double sum_terms = 0.0;
  for (int i = 0; i < N; ++i) sum_terms += groups.integrals[i].term;

  // Gauss-Newton joint Hessian: the local models use the GN inner block,
  // and the stage-2 implicit chain must differentiate the same matrix it inverts.
  const bool want_deriv2 = false;
  List out_mode = as<List>(joint_cb(
      Rcpp::Named("pars")       = pars_full,
      Rcpp::Named("fixed")      = fixed,
      Rcpp::Named("deriv")      = true,
      Rcpp::Named("deriv2")     = want_deriv2,
      Rcpp::Named("conditions") = subject_meta["conditions"]));
  // the joint objective at the modes: the inner values of the subjects plus the data of
  // the conditions that belong to no subject
  const double OFV = as<double>(out_mode["value"]) + sum_terms + groups.constant;
  NumericVector grad_full = as<NumericVector>(out_mode["gradient"]);
  NumericMatrix H_full    = as<NumericMatrix>(out_mode["hessian"]);
  CharacterVector grad_names = grad_full.names();
  List Hdim = H_full.attr("dimnames");
  CharacterVector H_rownames = as<CharacterVector>(Hdim[0]);

  std::vector<int> outer_idx_grad = name_indices(grad_names, outer_names,
                                                 "fast_outer/grad");
  std::vector<int> outer_idx_hess = name_indices(H_rownames, outer_names,
                                                 "fast_outer/hess");
  NumericVector grad_outer(n_outer);
  for (int i = 0; i < n_outer; ++i) grad_outer[i] = grad_full[outer_idx_grad[i]];
  grad_outer.attr("names") = outer_names;

  // the outer model: the normal block of every group eliminated from the joint Gauss-Newton Hessian,
  // and the eta gradient g_i the inner search leaves taken through the coupling, grad -= H_oi H_i^-1 g_i,
  // the gradient at the exact modes to first order
  std::vector<double> Hess_col(n_outer * n_outer);
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r) Hess_col[r + c * n_outer] = H_full(outer_idx_hess[r], outer_idx_hess[c]);
  std::vector<double> H_oi, g_N;
  // with `newton_value` the value at the exact modes to second order in the eta gradient of the joint
  // objective, OFV - g_N' H_NN^-1 g_N / 2; the volume term stays that at the modes found
  double decrement = 0.0;
  for (int i = 0; i < N; ++i) {
    const dmod::nlme::LocalModel& m = groups.models[i];
    CharacterVector eta_nm_i = as<CharacterVector>(eta_names_list[i]);
    std::vector<int> eh = name_indices(H_rownames, eta_nm_i, "fast_outer/eta_hess");
    std::vector<int> eg = name_indices(grad_names, eta_nm_i, "fast_outer/eta_grad");
    const int n = m.normal.size();
    H_oi.assign((std::size_t) n_outer * n, 0.0);
    g_N.assign(n, 0.0);
    for (int k = 0; k < n; ++k) {
      g_N[k] = grad_full[eg[m.normal[k]]];
      for (int r = 0; r < n_outer; ++r) H_oi[r + k * n_outer] = H_full(outer_idx_hess[r], eh[m.normal[k]]);
    }
    dmod::nlme::eliminate_normal(n_outer, m, H_oi.data(), Hess_col.data(), grad_outer.begin(), g_N.data());
    if (newton_value) dmod::nlme::eliminate_normal(1, m, g_N.data(), &decrement, nullptr, g_N.data());
  }
  NumericMatrix Hess(n_outer, n_outer);
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r) Hess(r, c) = Hess_col[r + c * n_outer];
  Hess.attr("dimnames") = List::create(outer_names, outer_names);
  List H_inv_list = inner_res["H_inv"], weight_list = inner_res["weight"];

  List eta_shift;
  if (apply_correction && correction_mode != "none" && correction_cb_opt.isNotNull()) {
    Function correction_cb(correction_cb_opt);
    eta_shift = add_correction(correction_cb, outer_names, pars_full, H_full, H_inv_list, weight_list,
                               grad_outer, Hess);
  }

  return List::create(
    Rcpp::Named("value")       = OFV + decrement / 2.0,
    Rcpp::Named("gradient")    = grad_outer,
    Rcpp::Named("hessian")     = Hess,
    Rcpp::Named("eta_modes")   = inner_res["eta_modes"],
    Rcpp::Named("H_GN")        = inner_res["H_GN"],
    Rcpp::Named("H_inv")       = H_inv_list,
    Rcpp::Named("weight")      = weight_list,
    Rcpp::Named("log_det_H")   = inner_res["log_det_H"],
    Rcpp::Named("sum_logdetH") = sum_terms,
    Rcpp::Named("iterations")  = inner_res["iterations"],
    Rcpp::Named("converged")   = inner_res["converged"],
    Rcpp::Named("pars_full")   = pars_full,
    Rcpp::Named("H_full")      = H_full,
    Rcpp::Named("eta_shift")   = eta_shift);
}


// The two halves of an outer evaluation: `trial` the value at theta from the modes started at
// `warm`, `complete` the gradient, the Hessian of the model and the mode shifts at an accepted
// trial; a failure of either is a rejected trial.
template <class Model>
struct OuterEval {
  std::function<double(NumericVector&, NumericMatrix&, RObject&)> trial;
  std::function<void(RObject&, std::vector<double>&, typename Model::Hess&, List&)> complete;
  std::function<NumericMatrix(const RObject&)> modes;
};

// The outer trust region of FOCEI on a curvature model of trust_model.h, one loop for dense and
// sparse: Moré-Sorensen steps accepted by the ratio of actual to predicted reduction, until the model
// gains less than mtol in the region rmax or five rejected trials leave the value flat within ftol.
template <class Model>
List focei_outer_loop(OuterEval<Model>& ev, NumericVector init, CharacterVector outer_names,
                      List outer_ctrl, double inner_rinit, int N, int K, NumericMatrix eta_init) {
  const int n_outer = init.size();
  double rinit   = outer_ctrl.containsElementNamed("rinit") ? as<double>(outer_ctrl["rinit"]) : 1.0;
  double rmax    = outer_ctrl.containsElementNamed("rmax") ? as<double>(outer_ctrl["rmax"]) : 10.0;
  int    iterlim = outer_ctrl.containsElementNamed("iterlim") ? as<int>(outer_ctrl["iterlim"]) : 100;
  double ftol    = ctrl_tol(outer_ctrl, "ftol", 1e-6);
  double mtol    = ctrl_tol(outer_ctrl, "mtol", 1e-6);

  NumericVector theta(clone(init));
  theta.attr("names") = outer_names;
  NumericMatrix eta_warm = eta_init.nrow() == N && eta_init.ncol() == K ? clone(eta_init) : NumericMatrix(N, K);
  Model model(n_outer);
  // the inner iterations of all trials, summed over the groups
  int inner_iter = 0;
  auto count = [&](const RObject& T) {
    if (T.isNULL()) return;
    List out(T);
    if (out.containsElementNamed("iterations")) for (int it : as<IntegerVector>(out["iterations"])) inner_iter += it;
  };
  const std::vector<double> ones(n_outer, 1.0), zeros(n_outer, 0.0);

  RObject keep;
  double f_curr = ev.trial(theta, eta_warm, keep);
  count(keep);
  std::vector<double> g_curr;
  typename Model::Hess H_curr;
  List shift_curr;
  ev.complete(keep, g_curr, H_curr, shift_curr);
  // a copy: shift_modes() moves the warm start in place, the modes of the evaluation stay
  eta_warm = clone(ev.modes(keep));
  model.build(H_curr, ones.data(), 1.0, zeros.data());

  double r = rinit;
  bool   converged = false;
  int    iter = 0, n_stall = 0;
  NumericVector theta_try(n_outer);
  theta_try.attr("names") = outer_names;
  std::vector<double> p_step(n_outer), p_max(n_outer);
  std::vector<double> trace_value, trace_radius, trace_rho;
  std::vector<int>    trace_accept;
  trace_value .push_back(f_curr);
  trace_radius.push_back(r);
  trace_rho   .push_back(0.0);
  trace_accept.push_back(1);

  for (iter = 1; iter <= iterlim; ++iter) {
    R_CheckUserInterrupt();
    bool is_newton = false, is_hard = false, is_easy = false, newton_max = false;
    double predicted_red = 0.0;
    model.solve(g_curr.data(), r, p_step.data(), &is_newton, &is_hard, &is_easy, &predicted_red);
    // what the model can gain at all, in the largest region the fit allows: a region that rejections
    // shrank makes the step's own prediction small without the fit being near the optimum
    double reachable = predicted_red;
    if (r < rmax) model.solve(g_curr.data(), rmax, p_max.data(), &newton_max, &is_hard, &is_easy, &reachable);

    for (int j = 0; j < n_outer; ++j) theta_try[j] = theta[j] + p_step[j];

    NumericMatrix eta_warm_saved(clone(eta_warm));
    shift_modes(eta_warm, shift_curr, p_step, inner_rinit);
    RObject T;
    double f_try = std::numeric_limits<double>::infinity();
    bool eval_ok = true;
    try {
      f_try = ev.trial(theta_try, eta_warm, T);
      count(T);
    } catch (std::exception&) { eval_ok = false; }

    double actual_red = f_curr - f_try;
    double rho = (predicted_red > 0.0 && std::isfinite(f_try))
                   ? actual_red / predicted_red
                   : -std::numeric_limits<double>::infinity();

    bool is_terminate = std::isfinite(f_try) && std::fabs(reachable) < mtol &&
                        (std::fabs(actual_red) < ftol || std::fabs(predicted_red) < mtol);
    bool accept;
    if (is_terminate) {
      accept = (f_try < f_curr);
    } else if (rho < 0.25) {
      accept = false;
      r = shrunk_radius(r, p_step.data(), n_outer);
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    // the derivatives of an accepted trial that the loop goes on from; a failure rejects it
    std::vector<double> g_try;
    typename Model::Hess H_try;
    List shift_try;
    if (accept && eval_ok) {
      try {
        ev.complete(T, g_try, H_try, shift_try);
      } catch (std::exception&) { accept = false; r = r * 0.25; }
    }

    if (accept && eval_ok) {
      for (int j = 0; j < n_outer; ++j) theta[j] = theta_try[j];
      f_curr = f_try;
      g_curr = std::move(g_try);
      H_curr = std::move(H_try);
      shift_curr = shift_try;
      eta_warm = clone(ev.modes(T));
      keep = T;
      model.build(H_curr, ones.data(), 1.0, zeros.data());
    } else {
      eta_warm = eta_warm_saved;
    }

    trace_value .push_back(f_curr);
    trace_radius.push_back(r);
    trace_rho   .push_back(rho);
    trace_accept.push_back(accept ? 1 : 0);

    if (is_terminate) { converged = true; break; }
    // as trust(): consecutive rejected trials that leave the objective flat within ftol, the noise
    // of the evaluation, end the fit with nothing left to resolve
    if (accept || !eval_ok || std::fabs(actual_red) >= ftol) n_stall = 0;
    else if (++n_stall >= dmod::trust_driver::kStallLimit) { converged = true; break; }
  }
  // a loop that reaches iterlim leaves iter one past it
  iter = std::min(iter, iterlim);

  List trace = List::create(
      Rcpp::Named("value")  = NumericVector(trace_value.begin(), trace_value.end()),
      Rcpp::Named("radius") = NumericVector(trace_radius.begin(), trace_radius.end()),
      Rcpp::Named("rho")    = NumericVector(trace_rho.begin(), trace_rho.end()),
      Rcpp::Named("accept") = IntegerVector(trace_accept.begin(), trace_accept.end()));
  return List::create(
      Rcpp::Named("argument")   = theta,
      Rcpp::Named("value")      = f_curr,
      Rcpp::Named("keep")       = keep,
      Rcpp::Named("iterations") = iter,
      Rcpp::Named("innerIterations") = inner_iter,
      Rcpp::Named("converged")  = converged,
      Rcpp::Named("trace")      = trace);
}


// Outer trust over focei_outer_objfn. Omega_inv and log|Omega| are rebuilt each
// outer iteration from the Cholesky parameters in the outer vector.
// [[Rcpp::export]]
List focei_run(Function model_cb,
                    Function err_cb,
                    Function joint_cb,
                    NumericVector init,
                    List subject_meta,
                    Nullable<NumericVector> fixed,
                    List control,
                    std::string correction_mode = "none",
                    Nullable<Function> correction_cb = R_NilValue,
                       Nullable<Function> predict_cb = R_NilValue,
                       Nullable<NumericMatrix> eta_init = R_NilValue) {
  CharacterVector subjects   = subject_meta["subjects"];
  CharacterVector outer_names = subject_meta["outer_names"];
  const int N = subjects.size();
  const int K = as<int>(subject_meta["K"]);
  const int n_outer = outer_names.size();
  if ((int)init.size() != n_outer)
    throw std::runtime_error("focei_run: init length mismatch.");

  List          om_meta    = as<List>(subject_meta["omega_meta"]);
  CharacterVector chol_pars  = as<CharacterVector>(om_meta["chol_pars"]);
  IntegerMatrix  chol_loc   = as<IntegerMatrix>(om_meta["chol_loc"]);
  LogicalVector  is_diag    = as<LogicalVector>(om_meta["is_diag"]);
  const int n_chol = chol_pars.size();

  // Build Omega from chol pars using exp-on-diag, value-on-off-diag, matching
  // omegaSpec$buildL.
  auto build_Omega_inv_and_logdet = [&](const NumericVector& theta)
      -> std::pair<NumericMatrix, double> {
    // Map chol par values from theta (named NumericVector)
    // a Cholesky parameter is estimated, in theta, or fixed
    CharacterVector tn = theta.names();
    NumericVector fx = fixed.isNotNull() ? NumericVector(fixed) : NumericVector(0);
    CharacterVector fn = fx.size() ? CharacterVector(fx.names()) : CharacterVector(0);
    NumericMatrix L(K, K);
    for (int m = 0; m < n_chol; ++m) {
      std::string nm = as<std::string>(chol_pars[m]);
      double v = NA_REAL;
      bool found = false;
      for (int j = 0; j < tn.size() && !found; ++j)
        if (as<std::string>(tn[j]) == nm) { v = theta[j]; found = true; }
      for (int j = 0; j < fn.size() && !found; ++j)
        if (as<std::string>(fn[j]) == nm) { v = fx[j]; found = true; }
      if (!found)
        throw std::runtime_error(std::string("focei_run: Cholesky parameter neither estimated nor fixed: ") + nm);
      int r = chol_loc(m, 0) - 1;
      int c = chol_loc(m, 1) - 1;
      L(r, c) = is_diag[m] ? std::exp(v) : v;
    }
    // Omega = L L^T, so Omega_inv = L^-T L^-1. Solve L X = I via forward sub.
    NumericMatrix L_inv(K, K);
    for (int j = 0; j < K; ++j) {
      L_inv(j, j) = 1.0 / L(j, j);
      for (int i = j + 1; i < K; ++i) {
        double s = 0.0;
        for (int q = j; q < i; ++q) s += L(i, q) * L_inv(q, j);
        L_inv(i, j) = -s / L(i, i);
      }
    }
    NumericMatrix Omega_inv(K, K);
    for (int r = 0; r < K; ++r)
      for (int c = 0; c < K; ++c) {
        double s = 0.0;
        for (int q = 0; q < K; ++q) s += L_inv(q, r) * L_inv(q, c);
        Omega_inv(r, c) = s;
      }
    double log_det = 0.0;
    for (int j = 0; j < K; ++j) log_det += 2.0 * std::log(L(j, j));
    return std::make_pair(Omega_inv, log_det);
  };

  List inner_ctrl = control.containsElementNamed("inner") ?
                      as<List>(control["inner"]) : List();
  const double inner_rinit = inner_ctrl.containsElementNamed("rinit") ? as<double>(inner_ctrl["rinit"]) : 1.0;
  List outer_ctrl = control.containsElementNamed("outer") ?
                      as<List>(control["outer"]) : List();
  const bool newton_value = control.containsElementNamed("newtonValue") && as<bool>(control["newtonValue"]);

  auto eval_outer = [&](NumericVector& th, NumericMatrix& warm, bool correct) -> List {
    auto OL = build_Omega_inv_and_logdet(th);
    return focei_outer_objfn(
        model_cb, err_cb, joint_cb, th, warm,
        subject_meta, OL.first, OL.second, fixed, inner_ctrl,
        correction_mode, correction_cb, predict_cb, correct, newton_value);
  };
  const bool corrected = correction_mode != "none" && correction_cb.isNotNull();
  // the trial point by the inner search with the outer gradient and Schur Hessian; at an accepted
  // one the volume correction
  OuterEval<dmod::trust_model::DenseModel> ev;
  ev.trial = [&](NumericVector& th, NumericMatrix& warm, RObject& T) {
    List out = eval_outer(th, warm, false);
    T = out;
    return as<double>(out["value"]);
  };
  ev.complete = [&](RObject& T, std::vector<double>& g, dmod::trust_model::DenseModel::Hess& H, List& shift) {
    List out(T);
    NumericVector gt = clone(as<NumericVector>(out["gradient"]));
    NumericMatrix Ht = clone(as<NumericMatrix>(out["hessian"]));
    if (corrected) {
      Function cb(correction_cb);
      shift = add_correction(cb, outer_names, as<NumericVector>(out["pars_full"]),
                             as<NumericMatrix>(out["H_full"]), as<List>(out["H_inv"]), as<List>(out["weight"]), gt, Ht);
    }
    out["gradient"] = gt;
    out["hessian"] = Ht;
    g.assign(gt.begin(), gt.end());
    H.H.assign(Ht.begin(), Ht.end());
    T = out;
  };
  ev.modes = [](const RObject& T) { return NumericMatrix(as<List>(T)["eta_modes"]); };
  List res = focei_outer_loop(ev, init, outer_names, outer_ctrl, inner_rinit, N, K,
                              eta_init.isNotNull() ? NumericMatrix(eta_init) : NumericMatrix(0, 0));

  List out_final = as<List>(res["keep"]);
  return List::create(
      Rcpp::Named("argument")   = res["argument"],
      Rcpp::Named("value")      = as<double>(out_final["value"]),
      Rcpp::Named("gradient")   = as<NumericVector>(out_final["gradient"]),
      Rcpp::Named("hessian")    = as<NumericMatrix>(out_final["hessian"]),
      Rcpp::Named("etaModes")   = as<NumericMatrix>(out_final["eta_modes"]),
      Rcpp::Named("H_GN")       = as<List>(out_final["H_GN"]),
      Rcpp::Named("H_inv")      = as<List>(out_final["H_inv"]),
      Rcpp::Named("weight")     = as<List>(out_final["weight"]),
      Rcpp::Named("log_det_H")  = as<NumericVector>(out_final["log_det_H"]),
      Rcpp::Named("sum_logdetH")= as<double>(out_final["sum_logdetH"]),
      Rcpp::Named("iterations") = res["iterations"],
      Rcpp::Named("innerIterations") = res["innerIterations"],
      Rcpp::Named("converged")  = res["converged"],
      Rcpp::Named("trace")      = res["trace"]);
}

// FOCEI on the sparse model by focei_outer_loop(), trial and derivatives from R: `trial(theta, warm)`
// a list of `value` and `eta_modes`, `complete(trial)` one of `gradient`, `kernel`, `hessian`,
// `lowrank` and `shift`, the mode shifts with their outer columns in "idx".
// [[Rcpp::export]]
List focei_run_sparse(Function trial, Function complete, NumericVector init, int N, int K,
                      List outer_ctrl, double inner_rinit, Nullable<NumericMatrix> eta_init = R_NilValue) {
  CharacterVector outer_names = init.names();
  const int n_outer = init.size();
  OuterEval<dmod::trust_model::SparseModel> ev;
  ev.trial = [&](NumericVector& th, NumericMatrix& warm, RObject& T) {
    List out = trial(th, warm);
    T = out;
    return as<double>(out["value"]);
  };
  ev.complete = [&](RObject& T, std::vector<double>& g, dmod::trust_model::SparseModel::Hess& H, List& shift) {
    List out = complete(T);
    NumericVector gt = out["gradient"];
    if (gt.size() != n_outer) throw std::runtime_error("focei_run_sparse: the gradient has the wrong length");
    if (!dmod::trust_model::SparseModel::read(out, n_outer, H))
      throw std::runtime_error("focei_run_sparse: no sparse Hessian");
    g.assign(gt.begin(), gt.end());
    shift = out.containsElementNamed("shift") ? List(out["shift"]) : List();
    T = out;
  };
  ev.modes = [](const RObject& T) { return NumericMatrix(as<List>(T)["eta_modes"]); };
  return focei_outer_loop(ev, init, outer_names, outer_ctrl, inner_rinit, N, K,
                          eta_init.isNotNull() ? NumericMatrix(eta_init) : NumericMatrix(0, 0));
}
