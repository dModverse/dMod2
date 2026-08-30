// dMod FOCEI kernel.
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

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "residual_kernel.h"
#include "trust_subproblem.h"

using namespace Rcpp;

// R >= 3.6.2 appends hidden Fortran string lengths to BLAS/LAPACK char args.
#ifndef FCONE
#define FCONE
#endif

// Tolerance from a control list, accepting the legacy trust() spelling.
static double ctrl_tol(const List& ctrl, const char* key, const char* legacy,
                       double dflt) {
  if (ctrl.containsElementNamed(key))    return as<double>(ctrl[key]);
  if (ctrl.containsElementNamed(legacy)) return as<double>(ctrl[legacy]);
  return dflt;
}


// Smoke-test entry: one joint() call, returns OFV + ncalls for boundary
// validation.
// [[Rcpp::export]]
List focei_kernel_ping(Function joint_cb,
                       NumericVector pars,
                       Nullable<NumericVector> fixed,
                       CharacterVector conditions) {
  List out = as<List>(joint_cb(
      Rcpp::Named("pars")       = pars,
      Rcpp::Named("fixed")      = fixed,
      Rcpp::Named("deriv")      = true,
      Rcpp::Named("conditions") = conditions));

  double value = as<double>(out["value"]);
  int    grad_len = 0;
  bool   has_hess = false;
  if (out.containsElementNamed("gradient") && !Rf_isNull(out["gradient"])) {
    NumericVector g = as<NumericVector>(out["gradient"]);
    grad_len = g.size();
  }
  if (out.containsElementNamed("hessian") && !Rf_isNull(out["hessian"])) {
    has_hess = true;
  }
  return List::create(
      Rcpp::Named("value")        = value,
      Rcpp::Named("gradient_len") = grad_len,
      Rcpp::Named("has_hessian")  = has_hess,
      Rcpp::Named("ncalls")       = 1);
}


static inline double vecs_at(const std::vector<double>& V, int K,
                             int row, int col) {
  return V[row + col * K];
}


// The Moré-Sorensen subproblem and its eigendecomposition are shared with
// trust_kernel.cpp; see trust_subproblem.h. Local signatures are kept so the
// call sites below read unchanged. The shared solver additionally handles the
// hard case (gradient orthogonal to the min-eigenspace), which the private
// version did not.
static inline void eigen_sym(double* A, int K, double* vals, double* vecs) {
  dmod::trust_internal::eigen_sym_local(A, K, vals, vecs);
}

static inline void trust_subproblem(int K, const double* g,
                                    const double* vals, const double* vecs,
                                    double r, double* p, double* predicted_red,
                                    bool* is_newton) {
  bool is_hard = false, is_easy = false;
  dmod::trust_internal::trust_sub(K, g, vals, vecs, r, p, predicted_red,
                                  is_newton, &is_hard, &is_easy);
}


// Per-subject inner trust loop. Calls joint(deriv=TRUE) per iteration;
// `pars` is mutated in place at `eta_idx_global` to reflect the current eta_i,
// other entries untouched.
struct InnerResult {
  std::vector<double> eta;        // K
  std::vector<double> grad;       // K (analytical gradient at the final eta)
  std::vector<double> H_GN;       // K*K column-major
  std::vector<double> H_inv;      // K*K column-major (Moore-Penrose with eigen floor)
  std::vector<double> eigvals;    // K (last-iter eigenvalues, ascending)
  double              log_det_H;
  double              value;
  int                 iterations;
  bool                converged;
};


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




// Stage-2 `d log|H| / d theta` correction is computed in R via
// `.computeFoceiCorrection` (R/nlme.R); the kernel calls it as an
// `Rcpp::Function` callback once per outer iter (see `correction_cb`).


// Fast inner objective for one subject. Bypasses normL2 / constraintL2 /
// evalConditionResidual / res / nll. Multi-output and eta-dependent sigma are
// supported: when sigma depends on eta (proportional/combined error models),
// the corresponding columns of attr(err, "deriv") feed Js_k = d sigma/d eta_k
// into the gradient and Hessian. The math is the M3-class likelihood (no BLOQ
// correction) with the GN-without-deriv2 convention: residual second
// derivatives of pred and sigma are dropped, but
// all first-derivative cross terms (Part1/2/3) are retained.
//
// Per call: ONE model_cb + ONE err_cb. `Omega_inv` is precomputed once per
// outer iter (chol pars don't change during inner trust) and reused.
struct EvalResult {
  double value;
  std::vector<double> grad;  // K
  std::vector<double> hess;  // K*K column-major
};

// Split in two so the lock-step driver can batch the callbacks: this half
// runs them for one subject, `eval_one_subject_from` is the pure math.
static void predict_one_subject(
    Function& model_cb, Function& err_cb, NumericVector& pars,
    Nullable<NumericVector> fixed, const List& meta_i, const double* eta_block,
    NumericMatrix* pred_out, NumericMatrix* err_out) {
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  const int K = static_cast<int>(eta_idx_in_pars.size());
  for (int k = 0; k < K; ++k) pars[eta_idx_in_pars[k] - 1] = eta_block[k];

  NumericVector model_times = meta_i["times"];
  CharacterVector condition = meta_i["condition"];
  List pred_list = as<List>(model_cb(
      Rcpp::Named("times")      = model_times,
      Rcpp::Named("pars")       = pars,
      Rcpp::Named("fixed")      = fixed,
      Rcpp::Named("deriv")      = true,
      Rcpp::Named("conditions") = condition));
  *pred_out = as<NumericMatrix>(pred_list[0]);
  NumericVector pinner = pred_out->attr("parameters");
  List err_list = as<List>(err_cb(
      Rcpp::Named("out")        = *pred_out,
      Rcpp::Named("pars")       = pinner,
      Rcpp::Named("conditions") = condition));
  *err_out = as<NumericMatrix>(err_list[0]);
}

static EvalResult eval_one_subject_from(
    NumericMatrix& pred_i,
    NumericMatrix& err_i,
    NumericVector& pars,
    const List& meta_i,
    const double* eta_block,
    const double* Omega_inv,
    double Omega_log_det) {
  const int K = static_cast<int>(as<IntegerVector>(meta_i["eta_idx_in_pars"]).size());
  IntegerVector eta_idx_in_pars      = meta_i["eta_idx_in_pars"];
  IntegerVector t_idx_in_pred        = meta_i["t_idx_in_pred"];
  NumericVector y_data               = meta_i["y_data"];
  IntegerVector o_idx_in_pred        = meta_i["o_idx_in_pred"];
  IntegerVector eta_idx_in_deriv     = meta_i["eta_idx_in_deriv"];
  IntegerVector o_idx_in_deriv       = meta_i["o_idx_in_deriv"];
  IntegerVector t_idx_in_err         = meta_i["t_idx_in_err"];
  IntegerVector o_idx_in_err         = meta_i["o_idx_in_err"];
  IntegerVector eta_idx_in_err_deriv = meta_i["eta_idx_in_err_deriv"];
  IntegerVector o_idx_in_err_deriv   = meta_i["o_idx_in_err_deriv"];

  const int T = y_data.size();

  for (int k = 0; k < K; ++k) {
    pars[eta_idx_in_pars[k] - 1] = eta_block[k];
  }

  NumericVector deriv_attr = as<NumericVector>(pred_i.attr("deriv"));
  IntegerVector deriv_dim  = deriv_attr.attr("dim");
  const int Dp0 = deriv_dim[0];  // model deriv: time
  const int Dp1 = deriv_dim[1];  // model deriv: observable

  // Err deriv attribute: present when sigma is differentiable w.r.t. outer
  // pars (i.e. virtually always under Y()-compiled errmodels). We index it
  // only when (a) the attribute exists, (b) the row's o_idx_in_err_deriv > 0,
  // and (c) the eta of interest has a column in the err deriv (> 0).
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

  double value = 0.0;
  std::vector<double> grad(K, 0.0);
  std::vector<double> hess(K * K, 0.0);

  // Pre-gather per-row predictions, sigmas, and Jacobian rows into packed
  // row-major buffers that the shared residual kernel consumes directly.
  // Stride is K (parameter dim) so each row's full Jacobian is contiguous.
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

    pred_row [row] = pred_i(tp, op);
    sigma_row[row] = err_i (te, oe);
    y_row    [row] = y_data[row];

    double* dpred_i = dpred_row.data() + (std::size_t) row * (std::size_t) K;
    for (int k = 0; k < K; ++k) {
      const int p = eta_idx_in_deriv[k] - 1;
      dpred_i[k] = deriv_attr[tp + od * Dp0 + p * Dp0 * Dp1];
    }

    // dsigma row: populated lazily on first eta-dependent row.
    const int oed_one = o_idx_in_err_deriv[row];
    if (have_err_deriv && oed_one > 0) {
      if (!sigma_has_eta) {
        dsigma_row.assign((std::size_t) T * (std::size_t) K, 0.0);
        sigma_has_eta = true;
      }
      double* dsigma_i = dsigma_row.data() + (std::size_t) row * (std::size_t) K;
      const int oed = oed_one - 1;
      for (int k = 0; k < K; ++k) {
        const int pe_one = eta_idx_in_err_deriv[k];
        if (pe_one > 0) {
          const int pe = pe_one - 1;
          dsigma_i[k] = err_deriv_attr[te + oed * De0 + pe * De0 * De1];
        }
      }
    }
  }

  // Delegate the ALOQ residual math (Part0+1+2+3, sigma(eta), value+grad+hess)
  // to the shared kernel. FOCEI never has BLOQ data, no exact-deriv2, and no
  // Bessel correction, so the opts struct is minimal.
  dmod::AccumOpts opts;
  opts.use_deriv2_exact    = false;
  opts.bloq_mode           = dmod::BloqMode::NONE;
  opts.sigma_depends_on_par = sigma_has_eta;
  // Hessian parts: all on, matching the original eval_one_subject behaviour.
  opts.aloq_part1 = opts.aloq_part2 = opts.aloq_part3 = true;

  dmod::accumulate_aloq_residual(
      /*n_obs =*/ T,
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
      grad.data(),
      hess.data());

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
  NumericMatrix pred_i, err_i;
  predict_one_subject(model_cb, err_cb, pars, fixed, meta_i, eta_block,
                      &pred_i, &err_i);
  return eval_one_subject_from(pred_i, err_i, pars, meta_i, eta_block,
                               Omega_inv, Omega_log_det);
}


// Rcpp-exported wrapper around eval_one_subject. Used by the joint-block
// Bayesian sampler (Pfad B / Particle-Gibbs) to evaluate the per-subject
// conditional posterior p(eta_i | y_i, theta, Omega) as an objlist:
//   value    = data NLL_i(eta_i) + eta_i^T Omega^-1 eta_i + log|Omega|
//   gradient = 2 J^T r / sigma^2 + 2 Omega^-1 eta_i (Fisher info + prior)
//   hessian  = 2 J^T (1/sigma^2) J + 2 Omega^-1     (Gauss-Newton + prior)
// All over the K eta entries of subject i.
//
// `pars_full` must contain placeholders for ALL parameters consumed by
// model_cb (typically prdfn = g * x * p): structural pars + the K eta_i
// for this subject. The eta slots are filled from `eta_block`.
// [[Rcpp::export]]
List focei_eval_one_subject(
    Function model_cb,
    Function err_cb,
    NumericVector pars_full,
    Nullable<NumericVector> fixed,
    List meta_i,
    NumericVector eta_block,
    NumericMatrix Omega_inv,
    double Omega_log_det) {

  NumericVector pars_copy = clone(pars_full);
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  const int K = eta_block.size();
  if ((int)eta_idx_in_pars.size() != K)
    stop("focei_eval_one_subject: eta_block length mismatch.");
  if (Omega_inv.nrow() != K || Omega_inv.ncol() != K)
    stop("focei_eval_one_subject: Omega_inv shape mismatch.");

  EvalResult er = eval_one_subject(
      model_cb, err_cb, pars_copy, fixed, meta_i,
      eta_block.begin(), Omega_inv.begin(), Omega_log_det);

  NumericVector grad(K);
  for (int k = 0; k < K; ++k) grad[k] = er.grad[k];
  NumericMatrix H(K, K);
  for (int c = 0; c < K; ++c) for (int r = 0; r < K; ++r)
    H(r, c) = er.hess[r + K * c];
  CharacterVector eta_names = meta_i["eta_names"];
  grad.attr("names") = eta_names;
  H.attr("dimnames") = List::create(eta_names, eta_names);

  return List::create(_["value"] = er.value,
                      _["gradient"] = grad,
                      _["hessian"]  = H);
}


// Per-subject inner trust using eval_one_subject. Returns the same
// InnerResult shape as inner_trust_one_subject().
// Eigen-floored log|H| and H_inv, shared by the serial and lock-step drivers.
static InnerResult finish_inner(int K, const std::vector<double>& eta_curr,
                                const std::vector<double>& grad,
                                const std::vector<double>& hess,
                                const std::vector<double>& eigvals,
                                const std::vector<double>& eigvecs,
                                double val_curr, int iter, bool converged,
                                double eigen_floor_relative) {
  double tr_H = 0.0;
  for (int j = 0; j < K; ++j) tr_H += eigvals[j];
  double eps_floor = eigen_floor_relative * std::fabs(tr_H) / K;
  std::vector<double> lam_floor(K);
  for (int j = 0; j < K; ++j) lam_floor[j] = std::max(eigvals[j], eps_floor);
  double log_det = 0.0;
  for (int j = 0; j < K; ++j) log_det += std::log(lam_floor[j]);
  std::vector<double> H_inv(K * K, 0.0);
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < K; ++j) {
      double s = 0.0;
      for (int q = 0; q < K; ++q)
        s += eigvecs[i + q * K] * eigvecs[j + q * K] / lam_floor[q];
      H_inv[i + j * K] = s;
    }

  InnerResult res;
  res.eta        = eta_curr;
  res.grad       = grad;
  res.H_GN       = hess;
  res.H_inv      = H_inv;
  res.eigvals    = eigvals;
  res.log_det_H  = log_det;
  res.value      = val_curr;
  res.iterations = iter;
  res.converged  = converged;
  return res;
}


static InnerResult inner_trust_one_subject(
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
    double ftol, double mtol,
    double eigen_floor_relative) {
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

  for (iter = 1; iter <= iterlim; ++iter) {
    bool   is_newton    = false;
    double predicted_red = 0.0;
    trust_subproblem(K, grad.data(), eigvals.data(), eigvecs.data(),
                     r, p_step.data(), &predicted_red, &is_newton);
    for (int k = 0; k < K; ++k) eta_try[k] = eta_curr[k] + p_step[k];

    EvalResult ft = eval_one_subject(
        model_cb, err_cb, pars, fixed, meta_i,
        eta_try.data(), Omega_inv, Omega_log_det);
    double val_try = ft.value;
    std::vector<double> grad_try(ft.grad), hess_try(ft.hess);
    std::vector<double> eigvals_try(K), eigvecs_try(K * K);
    eigen_sym(hess_try.data(), K, eigvals_try.data(), eigvecs_try.data());

    double actual_red = val_curr - val_try;
    double rho = (predicted_red > 0.0 && std::isfinite(val_try))
                   ? actual_red / predicted_red
                   : -std::numeric_limits<double>::infinity();

    bool is_terminate = std::isfinite(val_try) &&
                        (std::fabs(actual_red)   < ftol ||
                         std::fabs(predicted_red) < mtol);
    bool accept;
    if (is_terminate) {
      accept = (val_try < val_curr);
    } else if (rho < 0.25) {
      accept = false;
      r = r * 0.25;
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    if (accept) {
      eta_curr = eta_try;
      val_curr = val_try;
      grad     = grad_try;
      hess     = hess_try;
      eigvals  = eigvals_try;
      eigvecs  = eigvecs_try;
    }
    if (is_terminate) { converged = true; break; }
  }

  // Write final eta into pars (for downstream joint-at-modes call)
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  for (int k = 0; k < K; ++k) {
    pars[eta_idx_in_pars[k] - 1] = eta_curr[k];
  }

  return finish_inner(K, eta_curr, grad, hess, eigvals, eigvecs,
                     val_curr, iter, converged, eigen_floor_relative);
}


// Lock-step inner trust: every active subject takes one step per round, and
// the round's predictions go out as a single batched callback. The subjects'
// trust iterations are independent, so stepping them together is exact;
// converged subjects drop out so the round count follows the slowest one that
// is still moving.
static List inner_trust_lockstep(
    Function& predict_cb, NumericVector& pars_full, NumericMatrix& eta_warmstart,
    List& subject_meta_fast, CharacterVector& subjects,
    const std::vector<double>& Omega_inv, double Omega_log_det,
    Nullable<NumericVector> fixed, int K, int N,
    double rinit, double rmax, int iterlim, double ftol, double mtol,
    double eflr) {

  struct LS {
    std::vector<double> eta, grad, hess, eigvals, eigvecs, eta_try, p_step;
    double val = 0.0, r = 0.0, predicted_red = 0.0;
    int iter = 0;
    bool converged = false, done = false, is_newton = false;
    NumericVector pars;
  };
  std::vector<LS> st(N);

  // One batched callback for a set of subjects at given etas.
  auto run_round = [&](const std::vector<int>& idx,
                       const std::vector<const double*>& etas) {
    const int M = static_cast<int>(idx.size());
    List timesL(M), parsL(M);
    CharacterVector conds(M);
    for (int j = 0; j < M; ++j) {
      const int i = idx[j];
      List meta_i = as<List>(subject_meta_fast[i]);
      IntegerVector eip = meta_i["eta_idx_in_pars"];
      NumericVector pj = clone(st[i].pars);
      for (int k = 0; k < K; ++k) pj[eip[k] - 1] = etas[j][k];
      st[i].pars = pj;
      parsL[j]  = pj;
      timesL[j] = meta_i["times"];
      conds[j]  = as<CharacterVector>(meta_i["condition"])[0];
    }
    return as<List>(predict_cb(Rcpp::Named("timesList")  = timesL,
                               Rcpp::Named("parsList")   = parsL,
                               Rcpp::Named("conditions") = conds,
                               Rcpp::Named("fixed")      = fixed));
  };

  for (int i = 0; i < N; ++i) {
    st[i].pars = clone(pars_full);
    st[i].eta.assign(K, 0.0);
    for (int k = 0; k < K; ++k) st[i].eta[k] = eta_warmstart(i, k);
    st[i].eta_try.assign(K, 0.0);
    st[i].p_step.assign(K, 0.0);
    st[i].r = rinit;
  }

  {
    std::vector<int> idx(N);
    std::vector<const double*> etas(N);
    for (int i = 0; i < N; ++i) { idx[i] = i; etas[i] = st[i].eta.data(); }
    List got = run_round(idx, etas);
    List preds = as<List>(got["pred"]), errs = as<List>(got["err"]);
    for (int i = 0; i < N; ++i) {
      NumericMatrix pi = as<NumericMatrix>(preds[i]);
      NumericMatrix ei = as<NumericMatrix>(errs[i]);
      List meta_i = as<List>(subject_meta_fast[i]);
      EvalResult fr = eval_one_subject_from(pi, ei, st[i].pars, meta_i,
                                            st[i].eta.data(), Omega_inv.data(),
                                            Omega_log_det);
      st[i].val = fr.value; st[i].grad = fr.grad; st[i].hess = fr.hess;
      st[i].eigvals.assign(K, 0.0); st[i].eigvecs.assign((size_t)K * K, 0.0);
      eigen_sym(st[i].hess.data(), K, st[i].eigvals.data(), st[i].eigvecs.data());
    }
  }

  for (int round = 1; round <= iterlim; ++round) {
    std::vector<int> idx;
    std::vector<const double*> etas;
    for (int i = 0; i < N; ++i) {
      if (st[i].done) continue;
      st[i].iter = round;
      trust_subproblem(K, st[i].grad.data(), st[i].eigvals.data(),
                       st[i].eigvecs.data(), st[i].r, st[i].p_step.data(),
                       &st[i].predicted_red, &st[i].is_newton);
      for (int k = 0; k < K; ++k) st[i].eta_try[k] = st[i].eta[k] + st[i].p_step[k];
      idx.push_back(i);
      etas.push_back(st[i].eta_try.data());
    }
    if (idx.empty()) break;

    List got = run_round(idx, etas);
    List preds = as<List>(got["pred"]), errs = as<List>(got["err"]);

    for (size_t j = 0; j < idx.size(); ++j) {
      const int i = idx[j];
      NumericMatrix pi = as<NumericMatrix>(preds[j]);
      NumericMatrix ei = as<NumericMatrix>(errs[j]);
      List meta_i = as<List>(subject_meta_fast[i]);
      EvalResult ft = eval_one_subject_from(pi, ei, st[i].pars, meta_i,
                                            st[i].eta_try.data(),
                                            Omega_inv.data(), Omega_log_det);
      std::vector<double> eigvals_try(K), eigvecs_try((size_t)K * K);
      eigen_sym(ft.hess.data(), K, eigvals_try.data(), eigvecs_try.data());

      const double actual_red = st[i].val - ft.value;
      const double rho = (st[i].predicted_red > 0.0 && std::isfinite(ft.value))
                           ? actual_red / st[i].predicted_red
                           : -std::numeric_limits<double>::infinity();
      const bool is_terminate = std::isfinite(ft.value) &&
                                (std::fabs(actual_red) < ftol ||
                                 std::fabs(st[i].predicted_red) < mtol);
      bool accept;
      if (is_terminate) {
        accept = (ft.value < st[i].val);
      } else if (rho < 0.25) {
        accept = false; st[i].r *= 0.25;
      } else {
        accept = true;
        if (rho > 0.75 && !st[i].is_newton) st[i].r = std::min(2.0 * st[i].r, rmax);
      }
      if (accept) {
        st[i].eta = st[i].eta_try; st[i].val = ft.value;
        st[i].grad = ft.grad; st[i].hess = ft.hess;
        st[i].eigvals = eigvals_try; st[i].eigvecs = eigvecs_try;
      }
      if (is_terminate) { st[i].converged = true; st[i].done = true; }
    }
  }

  NumericMatrix eta_modes(N, K), grad_mat(N, K), eigvals_mat(N, K);
  NumericVector log_det_H(N), value(N);
  IntegerVector iterations(N);
  LogicalVector converged(N);
  List H_GN_list(N), H_inv_list(N);
  for (int i = 0; i < N; ++i) {
    IntegerVector eip = as<List>(subject_meta_fast[i])["eta_idx_in_pars"];
    for (int k = 0; k < K; ++k) pars_full[eip[k] - 1] = st[i].eta[k];
    InnerResult r = finish_inner(K, st[i].eta, st[i].grad, st[i].hess,
                                 st[i].eigvals, st[i].eigvecs, st[i].val,
                                 st[i].iter, st[i].converged, eflr);
    for (int k = 0; k < K; ++k) {
      eta_modes(i, k) = r.eta[k]; grad_mat(i, k) = r.grad[k];
      eigvals_mat(i, k) = r.eigvals[k];
    }
    NumericMatrix hg(K, K), hi(K, K);
    for (int c = 0; c < K; ++c) for (int rr = 0; rr < K; ++rr) {
      hg(rr, c) = r.H_GN[rr + c * K]; hi(rr, c) = r.H_inv[rr + c * K];
    }
    CharacterVector eta_nm =
      as<CharacterVector>(as<List>(subject_meta_fast[i])["eta_names"]);
    List dn = List::create(eta_nm, eta_nm);
    hg.attr("dimnames") = dn;
    hi.attr("dimnames") = dn;
    H_GN_list[i] = hg; H_inv_list[i] = hi;
    log_det_H[i] = r.log_det_H; value[i] = r.value;
    iterations[i] = r.iterations; converged[i] = r.converged;
  }
  rownames(eta_modes)   = subjects;
  rownames(grad_mat)    = subjects;
  rownames(eigvals_mat) = subjects;
  log_det_H .names() = subjects;
  value     .names() = subjects;
  iterations.names() = subjects;
  converged .names() = subjects;
  H_GN_list .names() = subjects;
  H_inv_list.names() = subjects;

  return List::create(
    Rcpp::Named("eta_modes")  = eta_modes,
    Rcpp::Named("gradient")   = grad_mat,
    Rcpp::Named("H_GN")       = H_GN_list,
    Rcpp::Named("H_inv")      = H_inv_list,
    Rcpp::Named("eigvals")    = eigvals_mat,
    Rcpp::Named("log_det_H")  = log_det_H,
    Rcpp::Named("value")      = value,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged")  = converged);
}


// Batched per-subject inner trust using the fast objective. Mirrors
// focei_inner_trust_batched but takes model_cb + err_cb + fast subject_meta
// entries instead of joint_cb. Omega_inv is precomputed once in R per outer
// iter (chol pars don't change inside inner trust).
// [[Rcpp::export]]
List focei_inner_trust(Function model_cb,
                            Function err_cb,
                            NumericVector pars_full,
                            NumericMatrix eta_warmstart,
                            List subject_meta,
                            NumericMatrix Omega_inv_mat,
                            double Omega_log_det,
                            Nullable<NumericVector> fixed,
                            List control,
                            Nullable<Function> predict_cb = R_NilValue) {
  CharacterVector subjects = subject_meta["subjects"];
  List subject_meta_fast   = subject_meta["fast_meta"];
  const int K = as<int>(subject_meta["K"]);
  const int N = subjects.size();
  if (eta_warmstart.nrow() != N || eta_warmstart.ncol() != K)
    throw std::runtime_error("focei_inner_trust: eta_warmstart shape mismatch.");
  if (Omega_inv_mat.nrow() != K || Omega_inv_mat.ncol() != K)
    throw std::runtime_error("focei_inner_trust: Omega_inv_mat shape mismatch.");

  std::vector<double> Omega_inv(K * K);
  for (int c = 0; c < K; ++c)
    for (int r = 0; r < K; ++r)
      Omega_inv[r + c * K] = Omega_inv_mat(r, c);

  double rinit   = control.containsElementNamed("rinit") ?
                     as<double>(control["rinit"]) : 1.0;
  double rmax    = control.containsElementNamed("rmax") ?
                     as<double>(control["rmax"]) : 10.0;
  int    iterlim = control.containsElementNamed("iterlim") ?
                     as<int>(control["iterlim"]) : 30;
  double ftol    = ctrl_tol(control, "ftol", "fterm", 1e-7);
  double mtol    = ctrl_tol(control, "mtol", "mterm", 1e-7);
  double eflr    = control.containsElementNamed("eigen_floor_relative") ?
                     as<double>(control["eigen_floor_relative"]) : 1e-10;

  if (predict_cb.isNotNull()) {
    Function pcb(predict_cb);
    return inner_trust_lockstep(pcb, pars_full, eta_warmstart, subject_meta_fast,
                                subjects, Omega_inv, Omega_log_det, fixed, K, N,
                                rinit, rmax, iterlim, ftol, mtol, eflr);
  }

  // Same "freeze other subjects' etas at warmstart" semantics as the
  // non-fast path: working copy per subject.
  std::vector<double> base_pars(pars_full.begin(), pars_full.end());

  NumericMatrix eta_modes(N, K);
  NumericMatrix grad_mat(N, K);
  NumericMatrix eigvals_mat(N, K);
  NumericVector log_det_H(N);
  NumericVector value(N);
  IntegerVector iterations(N);
  LogicalVector converged(N);
  List H_GN_list(N);
  List H_inv_list(N);

  for (int i = 0; i < N; ++i) {
    NumericVector eta_init(K);
    for (int k = 0; k < K; ++k) eta_init[k] = eta_warmstart(i, k);

    NumericVector working_pars(base_pars.begin(), base_pars.end());
    working_pars.attr("names") = pars_full.attr("names");

    List meta_i = as<List>(subject_meta_fast[i]);

    InnerResult r = inner_trust_one_subject(
        model_cb, err_cb, working_pars, fixed, meta_i, eta_init,
        Omega_inv.data(), Omega_log_det,
        rinit, rmax, iterlim, ftol, mtol, eflr);

    for (int k = 0; k < K; ++k) {
      eta_modes  (i, k) = r.eta[k];
      grad_mat   (i, k) = r.grad[k];
      eigvals_mat(i, k) = r.eigvals[k];
    }
    NumericMatrix H_GN_i (K, K);
    NumericMatrix H_inv_i(K, K);
    for (int c = 0; c < K; ++c)
      for (int rw = 0; rw < K; ++rw) {
        H_GN_i (rw, c) = r.H_GN [rw + c * K];
        H_inv_i(rw, c) = r.H_inv[rw + c * K];
      }
    CharacterVector eta_nm = as<CharacterVector>(as<List>(meta_i)["eta_names"]);
    List dn = List::create(eta_nm, eta_nm);
    H_GN_i .attr("dimnames") = dn;
    H_inv_i.attr("dimnames") = dn;
    H_GN_list [i] = H_GN_i;
    H_inv_list[i] = H_inv_i;
    log_det_H [i] = r.log_det_H;
    value     [i] = r.value;
    iterations[i] = r.iterations;
    converged [i] = r.converged;
  }

  // Write the converged modes into pars_full (caller-visible side effect)
  IntegerMatrix eta_idx_global = subject_meta["eta_idx_global"];
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < K; ++k)
      pars_full[eta_idx_global(i, k) - 1] = eta_modes(i, k);

  rownames(eta_modes)   = subjects;
  rownames(grad_mat)    = subjects;
  rownames(eigvals_mat) = subjects;
  log_det_H .names() = subjects;
  value     .names() = subjects;
  iterations.names() = subjects;
  converged .names() = subjects;
  H_GN_list .names() = subjects;
  H_inv_list.names() = subjects;

  return List::create(
    Rcpp::Named("eta_modes")  = eta_modes,
    Rcpp::Named("gradient")   = grad_mat,
    Rcpp::Named("H_GN")       = H_GN_list,
    Rcpp::Named("H_inv")      = H_inv_list,
    Rcpp::Named("eigvals")    = eigvals_mat,
    Rcpp::Named("log_det_H")  = log_det_H,
    Rcpp::Named("value")      = value,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged")  = converged);
}


// Fast variant of focei_outer_objfn: uses focei_inner_trust per-subject
// and recovers OFV by summing per-subject fast values. One joint_cb call is
// still required at modes for the OUTER structural gradient and Hessian
// (Schur complement needs cross-block H[outer, eta]).
// [[Rcpp::export]]
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
                       Nullable<Function> predict_cb = R_NilValue) {
  CharacterVector subjects        = subject_meta["subjects"];
  IntegerMatrix   eta_idx_global  = subject_meta["eta_idx_global"];
  List            eta_names_list  = subject_meta["eta_names"];
  CharacterVector outer_names     = subject_meta["outer_names"];
  IntegerVector   outer_idx_full  = subject_meta["outer_idx_in_full"];
  NumericVector   other_etas_init = subject_meta["other_etas_init"];
  CharacterVector pars_full_names = subject_meta["pars_full_names"];
  const int K = as<int>(subject_meta["K"]);
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

  // Fast inner trust (mutates pars_full to carry modes)
  List inner_res = focei_inner_trust(
      model_cb, err_cb, pars_full, eta_warmstart,
      subject_meta, Omega_inv_mat, Omega_log_det, fixed, inner_ctrl, predict_cb);
  NumericMatrix eta_modes  = as<NumericMatrix>(inner_res["eta_modes"]);
  List          H_GN_list  = as<List>(inner_res["H_GN"]);
  List          H_inv_list = as<List>(inner_res["H_inv"]);
  NumericVector log_det_H  = as<NumericVector>(inner_res["log_det_H"]);
  IntegerVector iterations = as<IntegerVector>(inner_res["iterations"]);
  LogicalVector converged  = as<LogicalVector>(inner_res["converged"]);
  NumericVector per_subj_value = as<NumericVector>(inner_res["value"]);

  // sum_i fast_value_i equals joint_at_modes.value by construction
  // (per-subject data nLL + quad + log|Omega| sums to normL2 + constraintL2_mvn).
  double joint_value = 0.0;
  for (int i = 0; i < N; ++i) joint_value += per_subj_value[i];
  double sum_logdetH = 0.0;
  for (int i = 0; i < N; ++i) sum_logdetH += log_det_H[i];
  const double OFV = joint_value + sum_logdetH
                    - (double)N * (double)K * std::log(2.0);

  // Only remaining full joint() call per outer iter, for the outer structural
  // gradient + Hessian.
  //
  // Gauss-Newton, deliberately. H_inv_list below is the inverse of the GN
  // inner Hessian, and the Laplace term uses log|H_GN|; the stage-2 implicit
  // chain differentiates the same matrix it inverts, so joint_hessian must be
  // the GN block too. Mixing an exact Newton block into that identity is
  // inconsistent and measurably worsens agreement with the exact marginal
  // (|grad| of bayesNLMEMarginal at the FOCEI mode: 0.024 GN vs 0.068 Newton
  // on tests/testthat/test-bayesNLME.R).
  //
  // This used to be `correction_mode != "none"`, which was silently ignored:
  // "+.objfn" dropped deriv2 before it reached the operands. Now that it does
  // not, the choice has to be stated.
  const bool want_deriv2 = false;
  List out_mode = as<List>(joint_cb(
      Rcpp::Named("pars")       = pars_full,
      Rcpp::Named("fixed")      = fixed,
      Rcpp::Named("deriv")      = true,
      Rcpp::Named("deriv2")     = want_deriv2,
      Rcpp::Named("conditions") = subjects));
  NumericVector grad_full = as<NumericVector>(out_mode["gradient"]);
  NumericMatrix H_full    = as<NumericMatrix>(out_mode["hessian"]);
  CharacterVector grad_names = grad_full.names();
  List Hdim = H_full.attr("dimnames");
  CharacterVector H_rownames = as<CharacterVector>(Hdim[0]);

  std::vector<int> outer_idx_grad = name_indices(grad_names, outer_names,
                                                 "fast_outer/grad");
  std::vector<int> outer_idx_hess = name_indices(H_rownames, outer_names,
                                                 "fast_outer/hess");
  std::vector<std::vector<int> > eta_idx_hess(N);
  for (int i = 0; i < N; ++i) {
    CharacterVector eta_nm_i = as<CharacterVector>(eta_names_list[i]);
    eta_idx_hess[i] = name_indices(H_rownames, eta_nm_i,
                                   "fast_outer/eta_hess");
  }

  NumericVector grad_outer(n_outer);
  for (int i = 0; i < n_outer; ++i) grad_outer[i] = grad_full[outer_idx_grad[i]];
  grad_outer.attr("names") = outer_names;

  // Stage-2 correction is honoured as eager-equivalent on the fast path:
  // lagged-Taylor bookkeeping is not wired through here, so we always apply
  // the current-iter analytical correction when the callback is provided.
  if (correction_mode != "none" && correction_cb_opt.isNotNull()) {
    Function correction_cb(correction_cb_opt);
    NumericVector corr = as<NumericVector>(correction_cb(
        Rcpp::Named("full_pars")     = pars_full,
        Rcpp::Named("joint_hessian") = H_full,
        Rcpp::Named("H_inv_list")    = H_inv_list));
    // Add by name (corr is named by outer params; missing names are 0)
    CharacterVector cnms = corr.names();
    for (int j = 0; j < cnms.size(); ++j) {
      std::string nm = as<std::string>(cnms[j]);
      for (int i = 0; i < n_outer; ++i)
        if (as<std::string>(outer_names[i]) == nm) {
          grad_outer[i] += corr[j];
          break;
        }
    }
  }

  // Schur Hessian. Per subject:
  //   tmp  = H_oi * H_inv_i              (n_outer x K)       via dgemm "N", "N"
  //   Hess = Hess - tmp * H_oi^T          (n_outer x n_outer) via dgemm "N", "T"
  NumericMatrix Hess(n_outer, n_outer);
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r)
      Hess(r, c) = H_full(outer_idx_hess[r], outer_idx_hess[c]);
  std::vector<double> H_oi(n_outer * K), tmp(n_outer * K);
  std::vector<double> H_inv_col(K * K);
  std::vector<double> Hess_col(n_outer * n_outer);
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r) Hess_col[r + c * n_outer] = Hess(r, c);
  const double one = 1.0, m_one = -1.0, zero = 0.0;
  for (int i = 0; i < N; ++i) {
    NumericMatrix H_inv_i = as<NumericMatrix>(H_inv_list[i]);
    // Pack H_inv_i to column-major buffer (it already is, but force layout)
    for (int c = 0; c < K; ++c)
      for (int r = 0; r < K; ++r) H_inv_col[r + c * K] = H_inv_i(r, c);
    // Pack H_oi from H_full
    for (int k = 0; k < K; ++k)
      for (int r = 0; r < n_outer; ++r)
        H_oi[r + k * n_outer] = H_full(outer_idx_hess[r], eta_idx_hess[i][k]);
    // tmp = H_oi * H_inv_i  (n_outer x K)
    F77_CALL(dgemm)("N", "N", &n_outer, &K, &K, &one,
                    H_oi.data(),     &n_outer,
                    H_inv_col.data(), &K,
                    &zero, tmp.data(), &n_outer FCONE FCONE);
    // Hess -= tmp * H_oi^T  (n_outer x n_outer)
    F77_CALL(dgemm)("N", "T", &n_outer, &n_outer, &K, &m_one,
                    tmp.data(),  &n_outer,
                    H_oi.data(), &n_outer,
                    &one, Hess_col.data(), &n_outer FCONE FCONE);
  }
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r) Hess(r, c) = Hess_col[r + c * n_outer];
  Hess.attr("dimnames") = List::create(outer_names, outer_names);

  return List::create(
    Rcpp::Named("value")       = OFV,
    Rcpp::Named("gradient")    = grad_outer,
    Rcpp::Named("hessian")     = Hess,
    Rcpp::Named("eta_modes")   = eta_modes,
    Rcpp::Named("H_GN")        = H_GN_list,
    Rcpp::Named("H_inv")       = H_inv_list,
    Rcpp::Named("log_det_H")   = log_det_H,
    Rcpp::Named("sum_logdetH") = sum_logdetH,
    Rcpp::Named("iterations")  = iterations,
    Rcpp::Named("converged")   = converged);
}


// Fast variant of focei_run: outer trust in C++ over the fast outer objfn.
// joint_cb is called once per outer iter for the structural gradient +
// Hessian; Omega_inv / Omega_log_det are recomputed each outer iter since
// chol pars sit in outer_pars (chol par locations come via subject_meta).
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
                       Nullable<Function> predict_cb = R_NilValue) {
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
    CharacterVector tn = theta.names();
    std::vector<int> idx(n_chol, -1);
    for (int m = 0; m < n_chol; ++m) {
      std::string nm = as<std::string>(chol_pars[m]);
      for (int j = 0; j < tn.size(); ++j)
        if (as<std::string>(tn[j]) == nm) { idx[m] = j; break; }
      if (idx[m] < 0)
        throw std::runtime_error(std::string("focei_run: chol par missing from outer pars: ") + nm);
    }
    NumericMatrix L(K, K);
    for (int m = 0; m < n_chol; ++m) {
      int r = chol_loc(m, 0) - 1;
      int c = chol_loc(m, 1) - 1;
      double v = theta[idx[m]];
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
  List outer_ctrl = control.containsElementNamed("outer") ?
                      as<List>(control["outer"]) : List();
  double rinit   = outer_ctrl.containsElementNamed("rinit") ?
                     as<double>(outer_ctrl["rinit"]) : 1.0;
  double rmax    = outer_ctrl.containsElementNamed("rmax") ?
                     as<double>(outer_ctrl["rmax"]) : 10.0;
  int    iterlim = outer_ctrl.containsElementNamed("iterlim") ?
                     as<int>(outer_ctrl["iterlim"]) : 100;
  double ftol    = ctrl_tol(outer_ctrl, "ftol", "fterm", 1e-6);
  double mtol    = ctrl_tol(outer_ctrl, "mtol", "mterm", 1e-6);

  NumericVector theta(clone(init));
  theta.attr("names") = outer_names;
  NumericMatrix eta_warm(N, K);

  auto eval_outer = [&](NumericVector& th) -> List {
    auto OL = build_Omega_inv_and_logdet(th);
    return focei_outer_objfn(
        model_cb, err_cb, joint_cb, th, eta_warm,
        subject_meta, OL.first, OL.second, fixed, inner_ctrl,
        correction_mode, correction_cb, predict_cb);
  };

  List out_curr = eval_outer(theta);
  double       f_curr = as<double>(out_curr["value"]);
  NumericVector g_curr = as<NumericVector>(out_curr["gradient"]);
  NumericMatrix H_curr = as<NumericMatrix>(out_curr["hessian"]);
  eta_warm = as<NumericMatrix>(out_curr["eta_modes"]);

  double r = rinit;
  bool   converged = false;
  int    iter = 0;
  NumericVector theta_try(n_outer);
  theta_try.attr("names") = outer_names;
  std::vector<double> p_step(n_outer);
  std::vector<double> trace_value, trace_radius, trace_rho;
  std::vector<int>    trace_accept;
  trace_value .push_back(f_curr);
  trace_radius.push_back(r);
  trace_rho   .push_back(0.0);
  trace_accept.push_back(1);

  for (iter = 1; iter <= iterlim; ++iter) {
    std::vector<double> H_buf(n_outer * n_outer);
    for (int c = 0; c < n_outer; ++c)
      for (int r0 = 0; r0 < n_outer; ++r0)
        H_buf[r0 + c * n_outer] = H_curr(r0, c);
    std::vector<double> vals(n_outer), vecs(n_outer * n_outer);
    eigen_sym(H_buf.data(), n_outer, vals.data(), vecs.data());

    std::vector<double> g_vec(g_curr.begin(), g_curr.end());
    bool   is_newton    = false;
    double predicted_red = 0.0;
    trust_subproblem(n_outer, g_vec.data(), vals.data(), vecs.data(),
                     r, p_step.data(), &predicted_red, &is_newton);

    for (int j = 0; j < n_outer; ++j) theta_try[j] = theta[j] + p_step[j];

    NumericMatrix eta_warm_saved(clone(eta_warm));
    List out_try;
    double f_try = std::numeric_limits<double>::infinity();
    bool eval_ok = true;
    try {
      out_try = eval_outer(theta_try);
      f_try = as<double>(out_try["value"]);
    } catch (std::exception&) { eval_ok = false; }

    double actual_red = f_curr - f_try;
    double rho = (predicted_red > 0.0 && std::isfinite(f_try))
                   ? actual_red / predicted_red
                   : -std::numeric_limits<double>::infinity();

    bool is_terminate = std::isfinite(f_try) &&
                        (std::fabs(actual_red) < ftol ||
                         std::fabs(predicted_red) < mtol);
    bool accept;
    if (is_terminate) {
      accept = (f_try < f_curr);
    } else if (rho < 0.25) {
      accept = false;
      r = r * 0.25;
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    if (accept && eval_ok) {
      for (int j = 0; j < n_outer; ++j) theta[j] = theta_try[j];
      f_curr = f_try;
      g_curr = as<NumericVector>(out_try["gradient"]);
      H_curr = as<NumericMatrix>(out_try["hessian"]);
      eta_warm = as<NumericMatrix>(out_try["eta_modes"]);
    } else {
      eta_warm = eta_warm_saved;
    }

    trace_value .push_back(f_curr);
    trace_radius.push_back(r);
    trace_rho   .push_back(rho);
    trace_accept.push_back(accept ? 1 : 0);

    if (is_terminate) { converged = true; break; }
  }

  List out_final = eval_outer(theta);
  NumericMatrix eta_modes_f = as<NumericMatrix>(out_final["eta_modes"]);
  List          H_GN_list_f = as<List>(out_final["H_GN"]);
  List          H_inv_list_f = as<List>(out_final["H_inv"]);
  NumericVector log_det_H_f = as<NumericVector>(out_final["log_det_H"]);

  List trace = List::create(
      Rcpp::Named("value")  = NumericVector(trace_value.begin(), trace_value.end()),
      Rcpp::Named("radius") = NumericVector(trace_radius.begin(), trace_radius.end()),
      Rcpp::Named("rho")    = NumericVector(trace_rho.begin(), trace_rho.end()),
      Rcpp::Named("accept") = IntegerVector(trace_accept.begin(), trace_accept.end()));

  return List::create(
      Rcpp::Named("argument")   = theta,
      Rcpp::Named("value")      = as<double>(out_final["value"]),
      Rcpp::Named("gradient")   = as<NumericVector>(out_final["gradient"]),
      Rcpp::Named("hessian")    = as<NumericMatrix>(out_final["hessian"]),
      Rcpp::Named("etaModes")   = eta_modes_f,
      Rcpp::Named("H_GN")       = H_GN_list_f,
      Rcpp::Named("H_inv")      = H_inv_list_f,
      Rcpp::Named("log_det_H")  = log_det_H_f,
      Rcpp::Named("sum_logdetH")= as<double>(out_final["sum_logdetH"]),
      Rcpp::Named("iterations") = iter,
      Rcpp::Named("converged")  = converged,
      Rcpp::Named("trace")      = trace);
}


// Generic C++ trust-region outer loop around an arbitrary R objective.
// Used by the C++ FOCEI adapter when correction != "none": the R-side
// emObjfn(correction = ...) closure already carries the Stage-2 math, so we
// just drive it with C++ trust step bookkeeping. objfn_cb must accept a
