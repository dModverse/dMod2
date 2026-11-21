// NLME group kernel: inner modes, local models and integrals of all groups of one evaluation; the
// interface and the formulas are in nlme_kernel.h.
//
// Math anchors (load-bearing, do not change without re-deriving):
//   - Inner objective per group, Gauss-Newton in eta: second derivatives of prediction and sigma
//     dropped, first-order cross terms of an eta-dependent sigma kept; the prior eta' P eta adds 2 P.
//   - Laplace on the normal block: int exp(-Phi/2) d eta ~ exp(-Phi(m)/2) (2 pi)^{|N|/2} |H_NN / 2|^-1/2,
//     so -2 log L_i ~ Phi_i(m_i) + log|H_NN| - |N| log 2 up to the constant of the normal density,
//     which the prior omits as well.

#include <Rcpp.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "residual_kernel.h"
#include "trust_subproblem.h"
#include "nlme_kernel.h"

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif

namespace dmod {
namespace nlme {

// Whether a trial is accepted: by the ratio of the actual to the predicted decrease, and near the
// modes, where a Newton step predicts less than 1e-6 of the value, also by a smaller gradient: the
// solver's error in the value there exceeds what the step gains, the gradient is still measured.
static bool accept_trial(double rho, bool is_newton, double predicted_red, double val, const std::vector<double>& grad,
                         const std::vector<double>& grad_try) {
  if (rho >= 0.25) return true;
  if (!is_newton || predicted_red > 1e-6 * std::max(1.0, std::fabs(val))) return false;
  double g0 = 0.0, g1 = 0.0;
  for (std::size_t k = 0; k < grad.size(); ++k) { g0 += grad[k] * grad[k]; g1 += grad_try[k] * grad_try[k]; }
  return g1 < g0;
}

// The modes are found where the Gauss-Newton step H^-1 g, on H and not on the secant model, is below
// `steptol`, as the volume term is linear in their error; or where a step on the boundary gains less
// than `mtol` of the value `val`, the region shrunk below the distance to the modes.
static bool modes_found(int K, const std::vector<double>& grad, const std::vector<double>& eigvals,
                        const std::vector<double>& eigvecs, double predicted_red, bool is_newton, double steptol,
                        double mtol, double val) {
  if (!is_newton) return std::fabs(predicted_red) <= mtol * std::max(1.0, std::fabs(val));
  double tr = 0.0;
  for (int k = 0; k < K; ++k) tr += std::fabs(eigvals[k]);
  const double floor = 1e-10 * tr / std::max(K, 1);
  double s2 = 0.0;
  for (int k = 0; k < K; ++k) {
    double c = 0.0;
    for (int j = 0; j < K; ++j) c += eigvecs[j + (std::size_t) K * k] * grad[j];
    const double step = c / std::max(eigvals[k], floor);
    s2 += step * step;
  }
  return std::sqrt(s2) <= steptol;
}


// Thin wrappers over the Moré-Sorensen subproblem and eigendecomposition in
// trust_subproblem.h, shared with trust_kernel.cpp and covering the hard case.
// A Hessian that is not finite stops the evaluation: LAPACK does not return on NaN input.
static inline void eigen_sym(double* A, int K, double* vals, double* vecs) {
  for (std::size_t i = 0; i < (std::size_t) K * K; ++i)
    if (!std::isfinite(A[i])) throw std::runtime_error("FOCEI: the Hessian is not finite.");
  dmod::trust_internal::eigen_sym_local(A, K, vals, vecs);
}

// The inner trust steps on the Gauss-Newton Hessian H plus S, the curvature of the prediction
// in the etas that Gauss-Newton drops, learned by SR1 from the gradients of accepted steps:
// (H_new + S) step = g_new - g_old after the update.
static void secant_update(int K, std::vector<double>& S, const std::vector<double>& step,
                          const std::vector<double>& g_old, const std::vector<double>& g_new,
                          const std::vector<double>& H_new) {
  std::vector<double> v(K);
  double vs = 0.0, vv = 0.0, ss = 0.0;
  for (int r = 0; r < K; ++r) {
    double hs = 0.0;
    for (int c = 0; c < K; ++c) hs += (H_new[r + c * K] + S[r + c * K]) * step[c];
    v[r] = g_new[r] - g_old[r] - hs;
  }
  for (int r = 0; r < K; ++r) { vs += v[r] * step[r]; vv += v[r] * v[r]; ss += step[r] * step[r]; }
  if (!(std::fabs(vs) > 1e-8 * std::sqrt(vv * ss))) return;
  for (int c = 0; c < K; ++c)
    for (int r = 0; r < K; ++r) S[r + c * K] += v[r] * v[c] / vs;
}

// Eigen decomposition of the step model H + S, or of H where H + S is not positive definite.
static void model_eigen(int K, const std::vector<double>& H, const std::vector<double>& S,
                        std::vector<double>& vals, std::vector<double>& vecs) {
  std::vector<double> B(H);
  for (std::size_t q = 0; q < B.size(); ++q) B[q] += S[q];
  eigen_sym(B.data(), K, vals.data(), vecs.data());
  if (vals[0] <= 0.0) {
    B = H;
    eigen_sym(B.data(), K, vals.data(), vecs.data());
  }
}

static inline void trust_subproblem(int K, const double* g,
                                    const double* vals, const double* vecs,
                                    double r, double* p, double* predicted_red,
                                    bool* is_newton) {
  bool is_hard = false, is_easy = false;
  dmod::trust_internal::trust_sub(K, g, vals, vecs, r, p, predicted_red,
                                  is_newton, &is_hard, &is_easy);
}


// Parameter names of a [time, observable, parameter] deriv array.
static CharacterVector deriv_par_names(const NumericVector& d) {
  const List dn = d.attr("dimnames");
  return as<CharacterVector>(dn[2]);
}

// Inner objective for one subject from one model_cb and one err_cb call, without
// the normL2 / constraintL2 stack. Gauss-Newton: second derivatives of pred and
// sigma are dropped, first-order cross terms of an eta-dependent sigma kept.
struct EvalResult {
  double value;
  std::vector<double> grad;  // K
  std::vector<double> hess;  // K*K column-major
};

static inline bool all_finite(const EvalResult& e) {
  if (!std::isfinite(e.value)) return false;
  for (double v : e.grad) if (!std::isfinite(v)) return false;
  for (double v : e.hess) if (!std::isfinite(v)) return false;
  return true;
}

// The prediction and error-model frames of a subject, one per condition block;
// `failed` where the prediction of one of the subject blocks failed.
struct SubjectFrames {
  std::vector<NumericMatrix> pred, err;
  bool failed = false;
};

// Split in two so the lock-step driver can batch the callbacks: this half
// runs them for the condition blocks of one subject, `eval_one_subject_from`
// is the pure math.
static void predict_one_subject(
    Function& model_cb, Function& err_cb, NumericVector& pars,
    Nullable<NumericVector> fixed, const List& meta_i, const double* eta_block,
    SubjectFrames* out) {
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  const int K = static_cast<int>(eta_idx_in_pars.size());
  for (int k = 0; k < K; ++k) pars[eta_idx_in_pars[k] - 1] = eta_block[k];

  List blocks = meta_i["blocks"];
  out->pred.assign(blocks.size(), NumericMatrix());
  out->err.assign(blocks.size(), NumericMatrix());
  for (int b = 0; b < blocks.size(); ++b) {
    List block = blocks[b];
    NumericVector model_times = block["times"];
    CharacterVector condition = block["condition"];
    List pred_list = as<List>(model_cb(
        Rcpp::Named("times")      = model_times,
        Rcpp::Named("pars")       = pars,
        Rcpp::Named("fixed")      = fixed,
        Rcpp::Named("deriv")      = true,
        Rcpp::Named("conditions") = condition));
    out->pred[b] = as<NumericMatrix>(pred_list[0]);
    NumericVector pinner = out->pred[b].attr("parameters");
    List err_list = as<List>(err_cb(
        Rcpp::Named("out")        = out->pred[b],
        Rcpp::Named("pars")       = pinner,
        Rcpp::Named("conditions") = condition));
    out->err[b] = as<NumericMatrix>(err_list[0]);
  }
}

// Data term of one condition block of a subject in the K subject etas, Gauss-Newton,
// added to value, grad and hess.
static void accumulate_block(
    NumericMatrix& pred_i,
    NumericMatrix& err_i,
    const List& meta_i,
    const CharacterVector& eta_names,
    int K,
    double& value,
    std::vector<double>& grad,
    std::vector<double>& hess) {
  IntegerVector t_idx_in_pred        = meta_i["t_idx_in_pred"];
  NumericVector y_data               = meta_i["y_data"];
  IntegerVector o_idx_in_pred        = meta_i["o_idx_in_pred"];
  IntegerVector o_idx_in_deriv       = meta_i["o_idx_in_deriv"];
  IntegerVector t_idx_in_err         = meta_i["t_idx_in_err"];
  IntegerVector o_idx_in_err         = meta_i["o_idx_in_err"];
  IntegerVector o_idx_in_err_deriv   = meta_i["o_idx_in_err_deriv"];
  // a recorded sigma replaces the error model of the row, NA where the model applies
  const bool have_sigma_data = meta_i.containsElementNamed("sigma_data");
  NumericVector sigma_data = have_sigma_data ? as<NumericVector>(meta_i["sigma_data"]) : NumericVector(0);

  const int T = y_data.size();

  NumericVector deriv_attr = as<NumericVector>(pred_i.attr("deriv"));
  IntegerVector deriv_dim  = deriv_attr.attr("dim");
  const int Dp0 = deriv_dim[0];  // model deriv: time
  const int Dp1 = deriv_dim[1];  // model deriv: observable
  // eta columns by name: the parameter axis depends on `fixed` and on the call; an eta
  // without a column enters only the error model
  std::vector<int> eta_col(K, -1);
  {
    const CharacterVector dpn = deriv_par_names(deriv_attr);
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < dpn.size(); ++j)
        if (as<std::string>(dpn[j]) == as<std::string>(eta_names[k])) { eta_col[k] = j; break; }
  }

  // The err deriv attribute is read only where it exists, the row has an
  // err-deriv observable and the eta has a column in it (indices > 0).
  Nullable<NumericVector> err_deriv_attr_opt = err_i.attr("deriv");
  const bool         have_err_deriv = err_deriv_attr_opt.isNotNull();
  NumericVector      err_deriv_attr;
  int                De0 = 0, De1 = 0;
  if (have_err_deriv) {
    err_deriv_attr = as<NumericVector>(err_deriv_attr_opt);
    IntegerVector dim_e = err_deriv_attr.attr("dim");
    De0 = dim_e[0];  // err deriv: time
    De1 = dim_e[1];  // err deriv: observable
  }
  // an eta without a column in the err deriv does not enter sigma
  std::vector<int> eta_col_err(K, -1);
  if (have_err_deriv) {
    const CharacterVector epn = deriv_par_names(err_deriv_attr);
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < epn.size(); ++j)
        if (as<std::string>(epn[j]) == as<std::string>(eta_names[k])) { eta_col_err[k] = j; break; }
  }

  // Per-row predictions, sigmas and Jacobian rows in packed row-major buffers
  // of stride K, the layout the shared residual kernel reads.
  std::vector<double> pred_row(T);
  std::vector<double> sigma_row(T);
  std::vector<double> y_row(T);
  std::vector<double> dpred_row((std::size_t) T * (std::size_t) K, 0.0);
  std::vector<double> dsigma_row;
  bool sigma_has_eta = false;  // any row with non-trivial err-deriv slice?

  for (int row = 0; row < T; ++row) {
    const int tp = t_idx_in_pred [row] - 1;
    const int op = o_idx_in_pred [row] - 1;
    const int od = o_idx_in_deriv[row] - 1;
    const int te = t_idx_in_err  [row] - 1;
    const int oe = o_idx_in_err  [row] - 1;

    const bool recorded = have_sigma_data && !ISNAN(sigma_data[row]);
    pred_row [row] = pred_i(tp, op);
    sigma_row[row] = recorded ? sigma_data[row] : err_i(te, oe);
    y_row    [row] = y_data[row];

    double* dpred_i = dpred_row.data() + (std::size_t) row * (std::size_t) K;
    for (int k = 0; k < K; ++k) {
      const int p = eta_col[k];
      dpred_i[k] = p < 0 ? 0.0 : deriv_attr[tp + od * Dp0 + p * Dp0 * Dp1];
    }

    // dsigma row: populated lazily on first eta-dependent row.
    const int oed_one = o_idx_in_err_deriv[row];
    if (!recorded && have_err_deriv && oed_one > 0) {
      if (!sigma_has_eta) {
        dsigma_row.assign((std::size_t) T * (std::size_t) K, 0.0);
        sigma_has_eta = true;
      }
      double* dsigma_i = dsigma_row.data() + (std::size_t) row * (std::size_t) K;
      const int oed = oed_one - 1;
      for (int k = 0; k < K; ++k) {
        const int pe = eta_col_err[k];
        if (pe >= 0) {
          dsigma_i[k] = err_deriv_attr[te + oed * De0 + pe * De0 * De1];
        }
      }
    }
  }

  // Shared residual kernel, Gauss-Newton in eta: the rows above the limit, then the
  // censored rows (M3), which the meta orders last. The plain sum of squares of these rows is not reported.
  const int n_aloq = meta_i.containsElementNamed("n_aloq") ? as<int>(meta_i["n_aloq"]) : T;
  double chi2_unused = 0.0;
  dmod::AccumOpts opts;
  opts.use_deriv2_exact    = false;
  opts.bloq_mode           = dmod::BloqMode::M3;
  if (meta_i.containsElementNamed("t_df")) opts.t_df = as<double>(meta_i["t_df"]);
  opts.sigma_depends_on_par = sigma_has_eta;
  // All first-order Hessian parts.
  opts.aloq_part1 = opts.aloq_part2 = opts.aloq_part3 = true;

  dmod::accumulate_aloq_residual(
      /*n_obs =*/ n_aloq,
      /*n_par =*/ K,
      pred_row.data(),
      dpred_row.data(),
      /*d2pred =*/ nullptr,
      y_row.data(),
      sigma_row.data(),
      sigma_has_eta ? dsigma_row.data() : nullptr,
      /*d2sigma=*/ nullptr,
      /*lloq   =*/ nullptr,
      opts,
      value,
      chi2_unused,
      grad.data(),
      hess.data());
  if (n_aloq < T) {
    NumericVector lloq = meta_i["lloq"];
    const std::size_t off = (std::size_t) n_aloq * (std::size_t) K;
    dmod::accumulate_bloq_residual(
        T - n_aloq, K,
        pred_row.data() + n_aloq,
        dpred_row.data() + off,
        /*d2pred =*/ nullptr,
        y_row.data() + n_aloq,
        sigma_row.data() + n_aloq,
        sigma_has_eta ? dsigma_row.data() + off : nullptr,
        /*d2sigma=*/ nullptr,
        lloq.begin() + n_aloq,
        opts,
        value,
        grad.data(),
        hess.data());
  }
}

// Inner objective of one subject: the data terms of the subject condition blocks plus
// the MVN prior of the subject etas.
static EvalResult eval_one_subject_from(
    SubjectFrames& frames,
    NumericVector& pars,
    const List& meta_i,
    const double* eta_block,
    const double* Omega_inv,
    double Omega_log_det) {
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  CharacterVector eta_names     = meta_i["eta_names"];
  const int K = static_cast<int>(eta_idx_in_pars.size());
  for (int k = 0; k < K; ++k) pars[eta_idx_in_pars[k] - 1] = eta_block[k];

  double value = 0.0;
  std::vector<double> grad(K, 0.0);
  std::vector<double> hess(K * K, 0.0);
  List blocks = meta_i["blocks"];
  for (int b = 0; b < blocks.size(); ++b)
    accumulate_block(frames.pred[b], frames.err[b], as<List>(blocks[b]), eta_names, K,
                     value, grad, hess);

  // MVN prior: quadratic + log|Omega| per subject; the constraintL2_mvn
  // convention drops the 2*pi factor (joint sums log|Omega| N times overall).
  std::vector<double> Oeta(K, 0.0);
  for (int j = 0; j < K; ++j)
    for (int k = 0; k < K; ++k)
      Oeta[k] += Omega_inv[k + j * K] * eta_block[j];
  for (int k = 0; k < K; ++k) value += eta_block[k] * Oeta[k];
  value += Omega_log_det;
  for (int k = 0; k < K; ++k) grad[k] += 2.0 * Oeta[k];
  for (int j = 0; j < K; ++j)
    for (int k = 0; k < K; ++k)
      hess[k + j * K] += 2.0 * Omega_inv[k + j * K];

  EvalResult out;
  out.value = value;
  out.grad  = std::move(grad);
  out.hess  = std::move(hess);
  return out;
}


static EvalResult eval_one_subject(
    Function& model_cb, Function& err_cb, NumericVector& pars,
    Nullable<NumericVector> fixed, const List& meta_i, const double* eta_block,
    const double* Omega_inv, double Omega_log_det) {
  SubjectFrames frames;
  predict_one_subject(model_cb, err_cb, pars, fixed, meta_i, eta_block, &frames);
  return eval_one_subject_from(frames, pars, meta_i, eta_block, Omega_inv, Omega_log_det);
}


// The local model from the state the search ends in, before the normal block is factored.
static LocalModel end_state(int K, const std::vector<double>& eta, const std::vector<double>& grad,
                            const std::vector<double>& hess, const std::vector<double>& eigvals,
                            const std::vector<double>& eigvecs, double value, int iterations, bool converged) {
  LocalModel m;
  m.K = K;
  m.eta = eta;
  m.grad = grad;
  m.H = hess;
  m.eigvals = eigvals;
  m.eigvecs = eigvecs;
  m.value = value;
  m.iterations = iterations;
  m.converged = converged;
  return m;
}


// Inner trust for one subject on eval_one_subject; writes the final eta into `pars`.
static LocalModel inner_trust_one_subject(
    Function& model_cb,
    Function& err_cb,
    NumericVector& pars,
    const Nullable<NumericVector>& fixed,
    const List& meta_i,
    const NumericVector& eta_init,
    const double* Omega_inv,
    double Omega_log_det,
    double rinit, double rmax,
    int    iterlim,
    double steptol, double mtol) {
  const int K = eta_init.size();

  std::vector<double> eta_curr(eta_init.begin(), eta_init.end());
  EvalResult fr = eval_one_subject(
      model_cb, err_cb, pars, fixed, meta_i,
      eta_curr.data(), Omega_inv, Omega_log_det);
  double val_curr = fr.value;
  std::vector<double> grad(fr.grad), hess(fr.hess);

  // Eigen on current hessian
  std::vector<double> eigvals(K), eigvecs(K * K);
  eigen_sym(hess.data(), K, eigvals.data(), eigvecs.data());

  double r = rinit;
  bool   converged = false;
  int    iter = 0;
  std::vector<double> eta_try(K), p_step(K);
  std::vector<double> S((std::size_t) K * K, 0.0), mvals(eigvals), mvecs(eigvecs);

  for (iter = 1; iter <= iterlim; ++iter) {
    bool   is_newton    = false;
    double predicted_red = 0.0;
    trust_subproblem(K, grad.data(), mvals.data(), mvecs.data(),
                     r, p_step.data(), &predicted_red, &is_newton);
    if (modes_found(K, grad, eigvals, eigvecs, predicted_red, is_newton, steptol, mtol, val_curr)) {
      converged = true;
      break;
    }
    for (int k = 0; k < K; ++k) eta_try[k] = eta_curr[k] + p_step[k];

    EvalResult ft = eval_one_subject(
        model_cb, err_cb, pars, fixed, meta_i,
        eta_try.data(), Omega_inv, Omega_log_det);
    // a trial whose value or derivatives are not finite is rejected like one that does not decrease
    if (!all_finite(ft)) {
      r = shrunk_radius(r, p_step.data(), K);
      std::fill(S.begin(), S.end(), 0.0);
      model_eigen(K, hess, S, mvals, mvecs);
      continue;
    }
    double val_try = ft.value;
    std::vector<double> grad_try(ft.grad), hess_try(ft.hess);
    std::vector<double> eigvals_try(K), eigvecs_try(K * K);
    eigen_sym(hess_try.data(), K, eigvals_try.data(), eigvecs_try.data());

    double actual_red = val_curr - val_try;
    double rho = (predicted_red > 0.0 && std::isfinite(val_try))
                   ? actual_red / predicted_red
                   : -std::numeric_limits<double>::infinity();

    const bool accept = std::isfinite(val_try) && accept_trial(rho, is_newton, predicted_red, val_curr, grad, grad_try);
    if (!accept) r = shrunk_radius(r, p_step.data(), K);
    else if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);

    if (accept) {
      secant_update(K, S, p_step, grad, grad_try, hess_try);
      eta_curr = eta_try;
      val_curr = val_try;
      grad     = grad_try;
      hess     = hess_try;
      eigvals  = eigvals_try;
      eigvecs  = eigvecs_try;
    } else {
      std::fill(S.begin(), S.end(), 0.0);
    }
    model_eigen(K, hess, S, mvals, mvecs);
  }

  // Final eta into pars for the joint call at the modes.
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  for (int k = 0; k < K; ++k) {
    pars[eta_idx_in_pars[k] - 1] = eta_curr[k];
  }

  // a loop that reaches iterlim leaves iter one past it
  return end_state(K, eta_curr, grad, hess, eigvals, eigvecs, val_curr, std::min(iter, iterlim), converged);
}
// Lock-step inner trust: every active subject takes one step per round and the
// round's predictions go out as one batched callback. Subjects are independent,
// so this is exact; converged ones drop out.
static std::vector<LocalModel> inner_trust_lockstep(
    Function& predict_cb, NumericVector& pars_full, NumericMatrix& eta_warmstart,
    List& subject_meta_fast, const std::vector<std::vector<int>>& coords,
    const std::vector<std::vector<double>>& P_sub, double Omega_log_det,
    Nullable<NumericVector> fixed, int N,
    double rinit, double rmax, int iterlim, double steptol, double mtol) {

  struct LS {
    std::vector<double> eta, grad, hess, eigvals, eigvecs, eta_try, p_step, S, mvals, mvecs;
    double val = 0.0, r = 0.0, predicted_red = 0.0;
    int iter = 0, K = 0;
    bool converged = false, done = false, is_newton = false;
    NumericVector pars;
  };
  std::vector<LS> st(N);

  // One batched callback for the condition blocks of a set of subjects at given
  // etas; the frames come back per subject.
  auto run_round = [&](const std::vector<int>& idx,
                       const std::vector<const double*>& etas) {
    const int M = static_cast<int>(idx.size());
    int n_req = 0;
    for (int j = 0; j < M; ++j)
      n_req += as<List>(as<List>(subject_meta_fast[idx[j]])["blocks"]).size();
    List timesL(n_req), parsL(n_req);
    CharacterVector conds(n_req);
    int q = 0;
    for (int j = 0; j < M; ++j) {
      const int i = idx[j];
      List meta_i = as<List>(subject_meta_fast[i]);
      IntegerVector eip = meta_i["eta_idx_in_pars"];
      NumericVector pj = clone(st[i].pars);
      for (int k = 0; k < st[i].K; ++k) pj[eip[k] - 1] = etas[j][k];
      st[i].pars = pj;
      List blocks = meta_i["blocks"];
      for (int b = 0; b < blocks.size(); ++b, ++q) {
        List block = blocks[b];
        parsL[q]  = pj;
        timesL[q] = block["times"];
        conds[q]  = as<CharacterVector>(block["condition"])[0];
      }
    }
    List got = as<List>(predict_cb(Rcpp::Named("timesList")  = timesL,
                                   Rcpp::Named("parsList")   = parsL,
                                   Rcpp::Named("conditions") = conds,
                                   Rcpp::Named("fixed")      = fixed));
    List preds = as<List>(got["pred"]), errs = as<List>(got["err"]);
    std::vector<SubjectFrames> frames(M);
    q = 0;
    for (int j = 0; j < M; ++j) {
      const int nb = as<List>(as<List>(subject_meta_fast[idx[j]])["blocks"]).size();
      for (int b = 0; b < nb; ++b, ++q) {
        SEXP pq = preds[q], eq = errs[q];
        if (TYPEOF(pq) != REALSXP || TYPEOF(eq) != REALSXP) { frames[j].failed = true; continue; }
        frames[j].pred.push_back(as<NumericMatrix>(pq));
        frames[j].err.push_back(as<NumericMatrix>(eq));
      }
    }
    return frames;
  };

  for (int i = 0; i < N; ++i) {
    const int K = st[i].K = coords[i].size();
    st[i].pars = clone(pars_full);
    st[i].eta.assign(K, 0.0);
    for (int k = 0; k < K; ++k) st[i].eta[k] = eta_warmstart(i, coords[i][k]);
    st[i].eta_try.assign(K, 0.0);
    st[i].p_step.assign(K, 0.0);
    st[i].r = rinit;
  }

  {
    std::vector<int> idx(N);
    std::vector<const double*> etas(N);
    for (int i = 0; i < N; ++i) { idx[i] = i; etas[i] = st[i].eta.data(); }
    std::vector<SubjectFrames> frames = run_round(idx, etas);
    for (int i = 0; i < N; ++i) {
      if (frames[i].failed)
        Rcpp::stop("FOCEI: the prediction fails at the modes the inner search starts from.");
      List meta_i = as<List>(subject_meta_fast[i]);
      const int K = st[i].K;
      EvalResult fr = eval_one_subject_from(frames[i], st[i].pars, meta_i,
                                            st[i].eta.data(), P_sub[i].data(),
                                            Omega_log_det);
      st[i].val = fr.value; st[i].grad = fr.grad; st[i].hess = fr.hess;
      st[i].eigvals.assign(K, 0.0); st[i].eigvecs.assign((size_t)K * K, 0.0);
      eigen_sym(st[i].hess.data(), K, st[i].eigvals.data(), st[i].eigvecs.data());
      st[i].S.assign((size_t)K * K, 0.0);
      st[i].mvals = st[i].eigvals; st[i].mvecs = st[i].eigvecs;
    }
  }

  for (int round = 1; round <= iterlim; ++round) {
    std::vector<int> idx;
    std::vector<const double*> etas;
    for (int i = 0; i < N; ++i) {
      if (st[i].done) continue;
      const int K = st[i].K;
      st[i].iter = round;
      trust_subproblem(K, st[i].grad.data(), st[i].mvals.data(),
                       st[i].mvecs.data(), st[i].r, st[i].p_step.data(),
                       &st[i].predicted_red, &st[i].is_newton);
      if (modes_found(K, st[i].grad, st[i].eigvals, st[i].eigvecs, st[i].predicted_red, st[i].is_newton, steptol,
                      mtol, st[i].val)) {
        st[i].converged = true; st[i].done = true; continue;
      }
      for (int k = 0; k < K; ++k) st[i].eta_try[k] = st[i].eta[k] + st[i].p_step[k];
      idx.push_back(i);
      etas.push_back(st[i].eta_try.data());
    }
    if (idx.empty()) break;

    std::vector<SubjectFrames> frames = run_round(idx, etas);

    for (size_t j = 0; j < idx.size(); ++j) {
      const int i = idx[j];
      const int K = st[i].K;
      // a trial where the prediction fails, or whose value or derivatives are not finite, is
      // rejected like one that does not decrease
      auto reject = [&]() {
        st[i].r = shrunk_radius(st[i].r, st[i].p_step.data(), K);
        std::fill(st[i].S.begin(), st[i].S.end(), 0.0);
        model_eigen(K, st[i].hess, st[i].S, st[i].mvals, st[i].mvecs);
      };
      if (frames[j].failed) { reject(); continue; }
      List meta_i = as<List>(subject_meta_fast[i]);
      EvalResult ft = eval_one_subject_from(frames[j], st[i].pars, meta_i,
                                            st[i].eta_try.data(),
                                            P_sub[i].data(), Omega_log_det);
      if (!all_finite(ft)) { reject(); continue; }
      std::vector<double> eigvals_try(K), eigvecs_try((size_t)K * K);
      eigen_sym(ft.hess.data(), K, eigvals_try.data(), eigvecs_try.data());

      const double actual_red = st[i].val - ft.value;
      const double rho = (st[i].predicted_red > 0.0 && std::isfinite(ft.value))
                           ? actual_red / st[i].predicted_red
                           : -std::numeric_limits<double>::infinity();
      const bool accept = std::isfinite(ft.value) &&
                          accept_trial(rho, st[i].is_newton, st[i].predicted_red, st[i].val, st[i].grad, ft.grad);
      if (!accept) st[i].r = shrunk_radius(st[i].r, st[i].p_step.data(), K);
      else if (rho > 0.75 && !st[i].is_newton) st[i].r = std::min(2.0 * st[i].r, rmax);
      if (accept) {
        secant_update(K, st[i].S, st[i].p_step, st[i].grad, ft.grad, ft.hess);
        st[i].eta = st[i].eta_try; st[i].val = ft.value;
        st[i].grad = ft.grad; st[i].hess = ft.hess;
        st[i].eigvals = eigvals_try; st[i].eigvecs = eigvecs_try;
      } else {
        std::fill(st[i].S.begin(), st[i].S.end(), 0.0);
      }
      model_eigen(K, st[i].hess, st[i].S, st[i].mvals, st[i].mvecs);
    }
  }

  std::vector<LocalModel> out(N);
  for (int i = 0; i < N; ++i) {
    const int K = st[i].K;
    IntegerVector eip = as<List>(subject_meta_fast[i])["eta_idx_in_pars"];
    for (int k = 0; k < K; ++k) pars_full[eip[k] - 1] = st[i].eta[k];
    out[i] = end_state(K, st[i].eta, st[i].grad, st[i].hess, st[i].eigvals, st[i].eigvecs, st[i].val,
                       st[i].iter, st[i].converged);
  }
  return out;
}


// log|A| and A^-1 of a symmetric n x n matrix from the eigendecomposition, eigenvalues floored at
// eflr |tr(A)| / n.
static void floored_inverse(int n, const std::vector<double>& eigvals, const std::vector<double>& eigvecs,
                            double eflr, std::vector<double>& A_inv, double& log_det) {
  double tr_H = 0.0;
  for (int j = 0; j < n; ++j) tr_H += eigvals[j];
  double eps_floor = eflr * std::fabs(tr_H) / n;
  std::vector<double> lam_floor(n);
  for (int j = 0; j < n; ++j) lam_floor[j] = std::max(eigvals[j], eps_floor);
  log_det = 0.0;
  for (int j = 0; j < n; ++j) log_det += std::log(lam_floor[j]);
  A_inv.assign((std::size_t) n * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int q = 0; q < n; ++q)
        s += eigvecs[i + q * n] * eigvecs[j + q * n] / lam_floor[q];
      A_inv[i + j * n] = s;
    }
}

void eliminate_normal(int n, const LocalModel& m, const double* C, double* T, double* t, const double* g) {
  int k = static_cast<int>(m.normal.size());
  if (n == 0 || k == 0) return;
  std::vector<double> tmp((std::size_t) n * k);
  const double one = 1.0, m_one = -1.0, zero = 0.0;
  // tmp = C H_NN^-1, T -= tmp C'
  F77_CALL(dgemm)("N", "N", &n, &k, &k, &one, C, &n, m.H_inv.data(), &k, &zero, tmp.data(), &n FCONE FCONE);
  F77_CALL(dgemm)("N", "T", &n, &n, &k, &m_one, tmp.data(), &n, C, &n, &one, T, &n FCONE FCONE);
  if (!t) return;
  for (int r = 0; r < n; ++r) {
    double d = 0.0;
    for (int q = 0; q < k; ++q) d += tmp[r + q * n] * g[q];
    t[r] -= d;
  }
}

// The normal block of `m` on the coordinates `isN`: the floored inverse and log-determinant, and
// with other coordinates S and r on them.
static void factor_normal(LocalModel& m, const std::vector<bool>& isN, double eflr) {
  const int K = m.K;
  m.normal.clear();
  m.other.clear();
  for (int j = 0; j < K; ++j) (isN[j] ? m.normal : m.other).push_back(j);
  const int n = m.normal.size(), l = m.other.size();
  if (n == K) {
    floored_inverse(K, m.eigvals, m.eigvecs, eflr, m.H_inv, m.log_det);
    return;
  }
  m.H_inv.clear();
  m.log_det = 0.0;
  if (n > 0) {
    std::vector<double> A((std::size_t) n * n), vals(n), vecs((std::size_t) n * n);
    for (int c = 0; c < n; ++c)
      for (int r = 0; r < n; ++r) A[r + c * n] = m.H[m.normal[r] + (std::size_t) K * m.normal[c]];
    eigen_sym(A.data(), n, vals.data(), vecs.data());
    floored_inverse(n, vals, vecs, eflr, m.H_inv, m.log_det);
  }
  m.S.assign((std::size_t) l * l, 0.0);
  m.r.assign(l, 0.0);
  std::vector<double> C((std::size_t) l * n), gN(n);
  for (int c = 0; c < l; ++c) {
    m.r[c] = m.grad[m.other[c]];
    for (int r = 0; r < l; ++r) m.S[r + c * l] = m.H[m.other[r] + (std::size_t) K * m.other[c]];
  }
  for (int c = 0; c < n; ++c) {
    gN[c] = m.grad[m.normal[c]];
    for (int r = 0; r < l; ++r) C[r + c * l] = m.H[m.other[r] + (std::size_t) K * m.normal[c]];
  }
  eliminate_normal(l, m, C.data(), m.S.data(), m.r.data(), gN.data());
}

void LaplaceNormal::integrate(const LocalModel& m, GroupIntegral& out) const {
  const int K = m.K, n = m.normal.size();
  out.term = m.log_det;
  if (n == K) {
    out.weight = m.H_inv;
    return;
  }
  out.weight.assign((std::size_t) K * K, 0.0);
  for (int c = 0; c < n; ++c)
    for (int r = 0; r < n; ++r) out.weight[m.normal[r] + (std::size_t) K * m.normal[c]] = m.H_inv[r + c * n];
}

double LaplaceNormal::constant(const std::vector<LocalModel>& models) const {
  double n = 0.0;
  for (const LocalModel& m : models) n += m.normal.size();
  return -n * std::log(2.0);
}

InnerControl inner_control(const List& control) {
  InnerControl c;
  if (control.containsElementNamed("rinit")) c.rinit = as<double>(control["rinit"]);
  if (control.containsElementNamed("rmax")) c.rmax = as<double>(control["rmax"]);
  if (control.containsElementNamed("iterlim")) c.iterlim = as<int>(control["iterlim"]);
  c.steptol = ctrl_tol(control, "steptol", "steptol", 1e-8);
  c.mtol = ctrl_tol(control, "mtol", "mterm", 1e-12);
  if (control.containsElementNamed("eigen_floor_relative"))
    c.eigen_floor_relative = as<double>(control["eigen_floor_relative"]);
  return c;
}

std::vector<int> subject_coords(const List& meta_i, int K) {
  std::vector<int> c;
  if (meta_i.containsElementNamed("coords")) {
    IntegerVector ci = meta_i["coords"];
    for (int j : ci) {
      if (j < 1 || j > K) throw std::runtime_error("nlme_group_kernel: coords out of range.");
      c.push_back(j - 1);
    }
    return c;
  }
  for (int j = 0; j < K; ++j) c.push_back(j);
  return c;
}

// The K x K matrix P on the coordinates `c`.
static std::vector<double> restrict_square(const std::vector<double>& P, int K, const std::vector<int>& c) {
  const int n = c.size();
  std::vector<double> out((std::size_t) n * n);
  for (int b = 0; b < n; ++b)
    for (int a = 0; a < n; ++a) out[a + (std::size_t) n * b] = P[c[a] + (std::size_t) K * c[b]];
  return out;
}

Groups evaluate_groups(const GroupProblem& prob, NumericVector& pars, const NumericMatrix& warm_in,
                       const Integrator& integrator) {
  CharacterVector subjects = prob.meta["subjects"];
  List subject_meta_fast = prob.meta["fast_meta"];
  const int K = as<int>(prob.meta["K"]);
  const int N = subjects.size();
  NumericMatrix warm(warm_in);
  if (warm.nrow() != N || warm.ncol() != K)
    throw std::runtime_error("nlme_group_kernel: eta_warmstart shape mismatch.");
  const InnerControl& c = prob.control;
  std::vector<std::vector<int>> coords(N);
  std::vector<std::vector<double>> P_sub(N);
  for (int i = 0; i < N; ++i) {
    coords[i] = subject_coords(as<List>(subject_meta_fast[i]), K);
    P_sub[i] = restrict_square(prob.Omega_inv, K, coords[i]);
  }
  Groups out;
  if (prob.predict_cb) {
    out.models = inner_trust_lockstep(*prob.predict_cb, pars, warm, subject_meta_fast, coords, P_sub,
                                      prob.Omega_log_det, prob.fixed, N, c.rinit, c.rmax, c.iterlim,
                                      c.steptol, c.mtol);
  } else {
    // the etas of the other groups stay at their warm start: a working copy per group
    std::vector<double> base_pars(pars.begin(), pars.end());
    out.models.resize(N);
    for (int i = 0; i < N; ++i) {
      NumericVector eta_init(coords[i].size());
      for (std::size_t k = 0; k < coords[i].size(); ++k) eta_init[k] = warm(i, coords[i][k]);
      NumericVector working_pars(base_pars.begin(), base_pars.end());
      working_pars.attr("names") = pars.attr("names");
      out.models[i] = inner_trust_one_subject(*prob.model_cb, *prob.err_cb, working_pars, prob.fixed,
                                              as<List>(subject_meta_fast[i]), eta_init, P_sub[i].data(),
                                              prob.Omega_log_det, c.rinit, c.rmax, c.iterlim, c.steptol, c.mtol);
    }
    for (int i = 0; i < N; ++i) {
      IntegerVector eip = as<List>(subject_meta_fast[i])["eta_idx_in_pars"];
      for (std::size_t k = 0; k < coords[i].size(); ++k) pars[eip[k] - 1] = out.models[i].eta[k];
    }
  }
  out.integrals.resize(N);
  std::vector<bool> isN(K, !prob.normal_given);
  for (int j : prob.normal) isN[j] = true;
  for (int i = 0; i < N; ++i) {
    LocalModel& m = out.models[i];
    m.group = i;
    m.coords = coords[i];
    std::vector<bool> sub(coords[i].size());
    for (std::size_t k = 0; k < coords[i].size(); ++k) sub[k] = isN[coords[i][k]];
    factor_normal(m, sub, c.eigen_floor_relative);
    integrator.integrate(m, out.integrals[i]);
  }
  out.constant = integrator.constant(out.models);
  return out;
}

// A column-major n x n buffer as a matrix with the names `nm` on both sides.
static NumericMatrix named_square(const std::vector<double>& x, int n, const CharacterVector& nm) {
  NumericMatrix out(n, n);
  std::copy(x.begin(), x.end(), out.begin());
  out.attr("dimnames") = List::create(nm, nm);
  return out;
}

// A vector with the names `nm`.
static NumericVector named_vector(const std::vector<double>& x, const CharacterVector& nm) {
  NumericVector out(x.begin(), x.end());
  out.attr("names") = nm;
  return out;
}

List groups_list(const Groups& groups, const List& meta) {
  CharacterVector subjects = meta["subjects"];
  List fast = meta["fast_meta"];
  const int K = as<int>(meta["K"]);
  const int N = subjects.size();
  NumericMatrix eta_modes(N, K), grad_mat(N, K), eigvals_mat(N, K);
  NumericVector log_det_H(N), value(N), term(N);
  IntegerVector iterations(N);
  LogicalVector converged(N), ok(N);
  List H_GN_list(N), H_inv_list(N), weight_list(N), S_list(N), r_list(N), dgrad(N), deta(N), hyper(N), mean(N),
      cov(N);
  bool other = false, derivs = false;
  for (int i = 0; i < N; ++i) {
    const LocalModel& m = groups.models[i];
    const GroupIntegral& gi = groups.integrals[i];
    std::fill(eigvals_mat.row(i).begin(), eigvals_mat.row(i).end(), NA_REAL);
    for (int k = 0; k < m.K; ++k) {
      eta_modes(i, m.coords[k]) = m.eta[k]; grad_mat(i, m.coords[k]) = m.grad[k]; eigvals_mat(i, k) = m.eigvals[k];
    }
    CharacterVector eta_nm = as<CharacterVector>(as<List>(fast[i])["eta_names"]);
    CharacterVector nN(m.normal.size()), nL(m.other.size());
    for (std::size_t j = 0; j < m.normal.size(); ++j) nN[j] = eta_nm[m.normal[j]];
    for (std::size_t j = 0; j < m.other.size(); ++j) nL[j] = eta_nm[m.other[j]];
    H_GN_list[i] = named_square(m.H, m.K, eta_nm);
    H_inv_list[i] = named_square(m.H_inv, m.normal.size(), nN);
    weight_list[i] = named_square(gi.weight, m.K, eta_nm);
    if (!m.other.empty()) {
      other = true;
      S_list[i] = named_square(m.S, m.other.size(), nL);
      NumericVector r(m.r.begin(), m.r.end());
      r.attr("names") = nL;
      r_list[i] = r;
    }
    if (!gi.grad.empty() || !gi.eta.empty()) {
      derivs = true;
      dgrad[i] = named_vector(gi.grad, eta_nm);
      deta[i] = named_vector(gi.eta, eta_nm);
      hyper[i] = NumericVector(gi.hyper.begin(), gi.hyper.end());
      mean[i] = named_vector(gi.mean, nL);
      if (!gi.cov.empty()) cov[i] = named_square(gi.cov, m.other.size(), nL);
    }
    log_det_H[i] = m.log_det; value[i] = m.value; term[i] = gi.term;
    iterations[i] = m.iterations; converged[i] = m.converged; ok[i] = gi.ok;
  }
  rownames(eta_modes) = subjects;
  rownames(grad_mat) = subjects;
  rownames(eigvals_mat) = subjects;
  log_det_H.names() = subjects;
  value.names() = subjects;
  term.names() = subjects;
  iterations.names() = subjects;
  converged.names() = subjects;
  H_GN_list.names() = subjects;
  H_inv_list.names() = subjects;
  weight_list.names() = subjects;
  List out = List::create(
    Rcpp::Named("eta_modes")  = eta_modes,
    Rcpp::Named("gradient")   = grad_mat,
    Rcpp::Named("H_GN")       = H_GN_list,
    Rcpp::Named("H_inv")      = H_inv_list,
    Rcpp::Named("eigvals")    = eigvals_mat,
    Rcpp::Named("log_det_H")  = log_det_H,
    Rcpp::Named("value")      = value,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("term")       = term,
    Rcpp::Named("weight")     = weight_list,
    Rcpp::Named("constant")   = groups.constant);
  if (other) {
    S_list.names() = subjects;
    r_list.names() = subjects;
    out["S"] = S_list;
    out["r"] = r_list;
  }
  if (derivs) {
    ok.names() = subjects;
    for (List* x : {&dgrad, &deta, &hyper, &mean, &cov}) x->names() = subjects;
    out["ok"] = ok;
    out["dgrad"] = dgrad;
    out["deta"] = deta;
    out["hyper"] = hyper;
    out["mean"] = mean;
    out["cov"] = cov;
  }
  return out;
}

}  // namespace nlme
}  // namespace dmod


// The group kernel with Laplace on the normal block for all groups at the start `eta_warmstart`;
// control$normal (1-based) names the normal coordinates where not all are. The modes are written
// into `pars_full`.
// [[Rcpp::export]]
List nlme_group_kernel(Function model_cb, Function err_cb, NumericVector pars_full, NumericMatrix eta_warmstart,
                       List subject_meta, NumericMatrix Omega_inv_mat, double Omega_log_det,
                       Nullable<NumericVector> fixed, List control, Nullable<Function> predict_cb = R_NilValue) {
  const int K = as<int>(subject_meta["K"]);
  if (Omega_inv_mat.nrow() != K || Omega_inv_mat.ncol() != K)
    throw std::runtime_error("nlme_group_kernel: Omega_inv_mat shape mismatch.");
  dmod::nlme::GroupProblem prob;
  prob.model_cb = &model_cb;
  prob.err_cb = &err_cb;
  Function pcb = predict_cb.isNotNull() ? Function(predict_cb) : model_cb;
  if (predict_cb.isNotNull()) prob.predict_cb = &pcb;
  prob.fixed = fixed;
  prob.meta = subject_meta;
  prob.Omega_inv.assign(Omega_inv_mat.begin(), Omega_inv_mat.end());
  prob.Omega_log_det = Omega_log_det;
  prob.control = dmod::nlme::inner_control(control);
  if (control.containsElementNamed("normal")) {
    prob.normal_given = true;
    IntegerVector nm = control["normal"];
    for (int j : nm) {
      if (j < 1 || j > K) throw std::runtime_error("nlme_group_kernel: control$normal out of range.");
      prob.normal.push_back(j - 1);
    }
  }
  dmod::nlme::LaplaceNormal laplace;
  dmod::nlme::Groups groups = dmod::nlme::evaluate_groups(prob, pars_full, eta_warmstart, laplace);
  return dmod::nlme::groups_list(groups, subject_meta);
}
