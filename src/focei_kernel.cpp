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

// Tolerance from a control list, accepting the deprecated trust() alias.
static double ctrl_tol(const List& ctrl, const char* key, const char* legacy,
                       double dflt) {
  if (ctrl.containsElementNamed(key))    return as<double>(ctrl[key]);
  if (ctrl.containsElementNamed(legacy)) return as<double>(ctrl[legacy]);
  return dflt;
}


static inline double vecs_at(const std::vector<double>& V, int K,
                             int row, int col) {
  return V[row + col * K];
}


// Thin wrappers over the Moré-Sorensen subproblem and eigendecomposition in
// trust_subproblem.h, shared with trust_kernel.cpp and covering the hard case.
static inline void eigen_sym(double* A, int K, double* vals, double* vecs) {
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

// The radius after a rejected step: a quarter of the radius or of the step, whichever is
// shorter, so that a rejected Newton step well inside the region is not tried again.
static double shrunk_radius(double r, const double* step, int n) {
  double ss = 0.0;
  for (int k = 0; k < n; ++k) ss += step[k] * step[k];
  return 0.25 * std::min(r, std::sqrt(ss));
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


// Result of one subject's inner trust loop.
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




// Parameter names of a [time, observable, parameter] deriv array.
static CharacterVector deriv_par_names(const NumericVector& d) {
  const List dn = d.attr("dimnames");
  return as<CharacterVector>(dn[2]);
}

// The derivative of sum_i log|H_i| in theta comes from R (`.foceiVolume`) through
// `correction_cb`, at the start, at accepted trials and at the final point; the
// trust ratio of a trial needs its value only.


// Inner objective for one subject from one model_cb and one err_cb call, without
// the normL2 / constraintL2 stack. Gauss-Newton: second derivatives of pred and
// sigma are dropped, first-order cross terms of an eta-dependent sigma kept.
struct EvalResult {
  double value;
  std::vector<double> grad;  // K
  std::vector<double> hess;  // K*K column-major
};

// The prediction and error-model frames of a subject, one per condition block.
struct SubjectFrames {
  std::vector<NumericMatrix> pred, err;
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

// Data term of one condition block of a subject in its K etas, Gauss-Newton,
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
  // a recorded sigma replaces the error model of its row, NA where the model applies
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
  // censored rows (M3), which the meta orders last. Its plain sum of squares is not reported.
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

// Inner objective of one subject: the data terms of its condition blocks plus
// the MVN prior of its etas.
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


// Conditional objective of subject i over its K etas as an objlist: data NLL
// + eta^T Omega^-1 eta + log|Omega|, Gauss-Newton Hessian. `pars_full` holds
// every model_cb parameter; its eta slots are overwritten from `eta_block`.
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


// Inner trust for one subject on eval_one_subject; writes the final eta into `pars`.
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
  std::vector<double> S((std::size_t) K * K, 0.0), mvals(eigvals), mvecs(eigvecs);

  for (iter = 1; iter <= iterlim; ++iter) {
    bool   is_newton    = false;
    double predicted_red = 0.0;
    trust_subproblem(K, grad.data(), mvals.data(), mvecs.data(),
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
      r = shrunk_radius(r, p_step.data(), K);
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

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
    if (is_terminate) { converged = true; break; }
    model_eigen(K, hess, S, mvals, mvecs);
  }

  // Final eta into pars for the joint call at the modes.
  IntegerVector eta_idx_in_pars = meta_i["eta_idx_in_pars"];
  for (int k = 0; k < K; ++k) {
    pars[eta_idx_in_pars[k] - 1] = eta_curr[k];
  }

  // a loop that reaches iterlim leaves iter one past it
  return finish_inner(K, eta_curr, grad, hess, eigvals, eigvecs,
                     val_curr, std::min(iter, iterlim), converged, eigen_floor_relative);
}


// Lock-step inner trust: every active subject takes one step per round and the
// round's predictions go out as one batched callback. Subjects are independent,
// so this is exact; converged ones drop out.
static List inner_trust_lockstep(
    Function& predict_cb, NumericVector& pars_full, NumericMatrix& eta_warmstart,
    List& subject_meta_fast, CharacterVector& subjects,
    const std::vector<double>& Omega_inv, double Omega_log_det,
    Nullable<NumericVector> fixed, int K, int N,
    double rinit, double rmax, int iterlim, double ftol, double mtol,
    double eflr) {

  struct LS {
    std::vector<double> eta, grad, hess, eigvals, eigvecs, eta_try, p_step, S, mvals, mvecs;
    double val = 0.0, r = 0.0, predicted_red = 0.0;
    int iter = 0;
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
      for (int k = 0; k < K; ++k) pj[eip[k] - 1] = etas[j][k];
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
        frames[j].pred.push_back(as<NumericMatrix>(preds[q]));
        frames[j].err.push_back(as<NumericMatrix>(errs[q]));
      }
    }
    return frames;
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
    std::vector<SubjectFrames> frames = run_round(idx, etas);
    for (int i = 0; i < N; ++i) {
      List meta_i = as<List>(subject_meta_fast[i]);
      EvalResult fr = eval_one_subject_from(frames[i], st[i].pars, meta_i,
                                            st[i].eta.data(), Omega_inv.data(),
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
      st[i].iter = round;
      trust_subproblem(K, st[i].grad.data(), st[i].mvals.data(),
                       st[i].mvecs.data(), st[i].r, st[i].p_step.data(),
                       &st[i].predicted_red, &st[i].is_newton);
      for (int k = 0; k < K; ++k) st[i].eta_try[k] = st[i].eta[k] + st[i].p_step[k];
      idx.push_back(i);
      etas.push_back(st[i].eta_try.data());
    }
    if (idx.empty()) break;

    std::vector<SubjectFrames> frames = run_round(idx, etas);

    for (size_t j = 0; j < idx.size(); ++j) {
      const int i = idx[j];
      List meta_i = as<List>(subject_meta_fast[i]);
      EvalResult ft = eval_one_subject_from(frames[j], st[i].pars, meta_i,
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
        accept = false; st[i].r = shrunk_radius(st[i].r, st[i].p_step.data(), K);
      } else {
        accept = true;
        if (rho > 0.75 && !st[i].is_newton) st[i].r = std::min(2.0 * st[i].r, rmax);
      }
      if (accept) {
        secant_update(K, st[i].S, st[i].p_step, st[i].grad, ft.grad, ft.hess);
        st[i].eta = st[i].eta_try; st[i].val = ft.value;
        st[i].grad = ft.grad; st[i].hess = ft.hess;
        st[i].eigvals = eigvals_try; st[i].eigvecs = eigvecs_try;
      } else {
        std::fill(st[i].S.begin(), st[i].S.end(), 0.0);
      }
      if (is_terminate) { st[i].converged = true; st[i].done = true; continue; }
      model_eigen(K, st[i].hess, st[i].S, st[i].mvals, st[i].mvecs);
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


// Inner trust for all subjects on the fast objective. Omega_inv is fixed during
// the inner loop, so R computes it once per outer iteration.
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

  // Other subjects' etas stay at their warmstart: working copy per subject.
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


// Adds the volume derivatives of `correction_cb`, a list of the gradient named by
// outer parameter and an optional Hessian with outer dimnames, to the outer model;
// the callback gets the outer gradient without them as `data_gradient`. Returns its
// `eta_shift`, d eta*/d theta per subject (K x n_outer), or an empty list.
static List add_correction(Function& correction_cb, const CharacterVector& outer_names,
                           const NumericVector& pars_full, const NumericMatrix& H_full,
                           const List& H_inv_list, NumericVector& grad_outer,
                           NumericMatrix& Hess) {
  List corr = as<List>(correction_cb(
      Rcpp::Named("full_pars")     = pars_full,
      Rcpp::Named("joint_hessian") = H_full,
      Rcpp::Named("H_inv_list")    = H_inv_list,
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

// Moves the modes in `eta` (N x K) to first order along the step: eta += S_i step.
static void shift_modes(NumericMatrix& eta, const List& shift, const std::vector<double>& step) {
  const int N = eta.nrow(), K = eta.ncol(), n = step.size();
  if (shift.size() != N) return;
  for (int i = 0; i < N; ++i) {
    NumericMatrix S = as<NumericMatrix>(shift[i]);
    if (S.nrow() != K || S.ncol() != n) continue;
    for (int k = 0; k < K; ++k) {
      double d = 0.0;
      for (int j = 0; j < n; ++j) d += S(k, j) * step[j];
      eta(i, k) += d;
    }
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
                       bool apply_correction = true) {
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

  // Fast inner trust (mutates pars_full to hold the modes)
  List inner_res = focei_inner_trust(
      model_cb, err_cb, pars_full, eta_warmstart,
      subject_meta, Omega_inv_mat, Omega_log_det, fixed, inner_ctrl, predict_cb);
  NumericMatrix eta_modes  = as<NumericMatrix>(inner_res["eta_modes"]);
  List          H_GN_list  = as<List>(inner_res["H_GN"]);
  List          H_inv_list = as<List>(inner_res["H_inv"]);
  NumericVector log_det_H  = as<NumericVector>(inner_res["log_det_H"]);
  IntegerVector iterations = as<IntegerVector>(inner_res["iterations"]);
  LogicalVector converged  = as<LogicalVector>(inner_res["converged"]);

  double sum_logdetH = 0.0;
  for (int i = 0; i < N; ++i) sum_logdetH += log_det_H[i];

  // Gauss-Newton joint Hessian: H_inv_list and log|H| use the GN inner block,
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
  const double OFV = as<double>(out_mode["value"]) + sum_logdetH
                    - (double)N * (double)K * std::log(2.0);
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
  std::vector<std::vector<int> > eta_idx_grad(N);
  for (int i = 0; i < N; ++i)
    eta_idx_grad[i] = name_indices(grad_names, as<CharacterVector>(eta_names_list[i]),
                                   "fast_outer/eta_grad");
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
    // the gradient at the exact modes to first order: the inner trust stops on the value, so
    // the eta gradient g_i left over enters through the coupling, grad -= H_oi H_i^-1 g_i
    for (int r = 0; r < n_outer; ++r) {
      double d = 0.0;
      for (int k = 0; k < K; ++k) d += tmp[r + k * n_outer] * grad_full[eta_idx_grad[i][k]];
      grad_outer[r] -= d;
    }
  }
  for (int c = 0; c < n_outer; ++c)
    for (int r = 0; r < n_outer; ++r) Hess(r, c) = Hess_col[r + c * n_outer];
  Hess.attr("dimnames") = List::create(outer_names, outer_names);

  List eta_shift;
  if (apply_correction && correction_mode != "none" && correction_cb_opt.isNotNull()) {
    Function correction_cb(correction_cb_opt);
    eta_shift = add_correction(correction_cb, outer_names, pars_full, H_full, H_inv_list,
                               grad_outer, Hess);
  }

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
    Rcpp::Named("converged")   = converged,
    Rcpp::Named("pars_full")   = pars_full,
    Rcpp::Named("H_full")      = H_full,
    Rcpp::Named("eta_shift")   = eta_shift);
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

  auto eval_outer = [&](NumericVector& th, bool correct) -> List {
    auto OL = build_Omega_inv_and_logdet(th);
    return focei_outer_objfn(
        model_cb, err_cb, joint_cb, th, eta_warm,
        subject_meta, OL.first, OL.second, fixed, inner_ctrl,
        correction_mode, correction_cb, predict_cb, correct);
  };
  const bool corrected = correction_mode != "none" && correction_cb.isNotNull();

  List out_curr = eval_outer(theta, true);
  // the evaluation of the point the loop stands on, with its volume correction
  List out_keep = out_curr;
  double       f_curr = as<double>(out_curr["value"]);
  NumericVector g_curr = as<NumericVector>(out_curr["gradient"]);
  NumericMatrix H_curr = as<NumericMatrix>(out_curr["hessian"]);
  eta_warm = as<NumericMatrix>(out_curr["eta_modes"]);
  List shift_curr = as<List>(out_curr["eta_shift"]);

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
    shift_modes(eta_warm, shift_curr, p_step);
    List out_try;
    double f_try = std::numeric_limits<double>::infinity();
    bool eval_ok = true;
    try {
      out_try = eval_outer(theta_try, false);
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
      r = shrunk_radius(r, p_step.data(), n_outer);
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    // the volume derivatives of an accepted trial that the loop goes on from;
    // a failure rejects the trial
    NumericVector g_try;
    NumericMatrix H_try;
    List shift_try;
    if (accept && eval_ok) {
      g_try = clone(as<NumericVector>(out_try["gradient"]));
      H_try = clone(as<NumericMatrix>(out_try["hessian"]));
      if (corrected) {
        try {
          Function cb(correction_cb);
          shift_try = add_correction(cb, outer_names, as<NumericVector>(out_try["pars_full"]),
                                     as<NumericMatrix>(out_try["H_full"]), as<List>(out_try["H_inv"]),
                                     g_try, H_try);
        } catch (std::exception&) { accept = false; r = r * 0.25; }
      }
    }

    if (accept && eval_ok) {
      for (int j = 0; j < n_outer; ++j) theta[j] = theta_try[j];
      f_curr = f_try;
      g_curr = g_try;
      H_curr = H_try;
      shift_curr = shift_try;
      eta_warm = as<NumericMatrix>(out_try["eta_modes"]);
      out_keep = out_try;
      out_keep["gradient"] = g_try;
      out_keep["hessian"] = H_try;
    } else {
      eta_warm = eta_warm_saved;
    }

    trace_value .push_back(f_curr);
    trace_radius.push_back(r);
    trace_rho   .push_back(rho);
    trace_accept.push_back(accept ? 1 : 0);

    if (is_terminate) { converged = true; break; }
  }
  // a loop that reaches iterlim leaves iter one past it
  iter = std::min(iter, iterlim);

  List out_final = out_keep;
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
