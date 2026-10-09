// dMod shared residual kernel.
//
// Math reference (load-bearing):
//
//   wr  = (pred - val) / sigma   ("weighted residual")
//   w0  = pred / sigma           ("weighted zero"; only used for M4)
//
//   For ALOQ rows (val > lloq):
//     obj   += wr^2 + log(2 pi sigma^2)
//     chi2  += wr^2                    (the data term without its normaliser)
//     grad  += 2 * (wr * dwr/dtheta + 1/sigma * dsigma/dtheta)
//     hess  += 2 * outer(dwr/dtheta, dwr/dtheta)              -- Part0 (GN)
//            + ALOQ_part1 + ALOQ_part2 + ALOQ_part3            -- toggleable
//            + (if use_deriv2_exact) 2 * wr/sigma * d2pred/dtheta^2
//            +                       (2/sigma) (1 - wr^2) * d2sigma/dtheta^2
//
//   For BLOQ rows (val <= lloq):
//     M3:   obj  += -2 log Phi(-wr)
//           grad += 2 G(-wr) dwr/dtheta
//           hess += 2 (-wr*G(-wr) + G(-wr)^2) outer(dwr, dwr) + parts2+3
//                + (if use_deriv2_exact) 2 G(-wr)/sigma * d2pred/dtheta^2
//                +                       -2 wr * G(-wr) / sigma * d2sigma/dtheta^2
//     M4*:  obj  += -2 log(1 - Phi(wr)/Phi(w0))               (with stability)
//           grad += 2 (c1 dwr - c2 dw0 + c3 dw0)
//           hess += corresponding 3-part GN form
//                + (if use_deriv2_exact) 2 (c1-c2+c3)/sigma * d2pred/dtheta^2
//                +                       -2/sigma (c1*wr + (c3-c2)*w0) * d2sigma/dtheta^2
//
//   M4BEAL additionally adds an exact ALOQ-side correction:
//     obj  += 2 log Phi(w0)
//     grad += 2 G(w0) dw0/dtheta
//     hess += 2 (-w0*G(w0) - G(w0)^2) outer(dw0, dw0) + cross-sigma term
//          + (if use_deriv2_exact) 2 G(w0)/sigma * d2pred/dtheta^2
//          +                       -2 w0 G(w0)/sigma * d2sigma/dtheta^2
//
// Conventions:
//   - dpred is stored row-major: dpred[row*n_par + p] = d pred_row / d par_p.
//   - d2pred is stored row-major over (row, p, q):
//     d2pred[row*n_par*n_par + p*n_par + q] = d^2 pred_row / d par_p d par_q.
//   - dsigma / d2sigma analogous shape (or nullptr if sigma is independent).
//   - hess accumulator is stored column-major (LAPACK convention).
//
// Threading: every function in this header is reentrant. The caller manages
// thread-local accumulators; concurrent calls touching the same hess buffer
// are NOT supported.
//
// All accumulator arrays are INOUT: caller is responsible for zeroing before
// the first call, and may call accumulate_aloq + accumulate_bloq in sequence
// to add contributions from both partitions.

#ifndef DMOD_RESIDUAL_KERNEL_H
#define DMOD_RESIDUAL_KERNEL_H

#include <cstddef>

namespace dmod {

enum class BloqMode {
  NONE,    // no BLOQ data at all (the default)
  M1,      // BLOQ rows exist but are excluded from the objective
  M3,      // -2 log Phi(-wr)
  M4NM,    // M4 method
  M4BEAL   // M4 with ALOQ-side correction term
};

struct AccumOpts {
  // If false, the Hessian accumulation is skipped entirely (value and gradient
  // still computed). hess_acc may then be null. Used by the quasi-Newton
  // Hessian sources, which maintain their own approximation.
  bool build_hessian = true;

  // If true, H += sum_i w_pred_i d^2 pred_i / d theta^2 + w_sig_i d^2 sigma_i /
  // d theta^2 with the per-row weights of the math reference above. The d2sigma
  // term needs sigma_depends_on_par and a non-null d2sigma buffer.
  bool use_deriv2_exact = false;

  // Selects BLOQ likelihood treatment. Has no effect on accumulate_aloq.
  // (M4BEAL however does add a correction to the ALOQ branch.)
  BloqMode bloq_mode = BloqMode::NONE;

  // If true, dsigma must be non-null and the kernel accounts for sigma's
  // dependence on parameters in gradient and Hessian (Parts 1/2/3).
  bool sigma_depends_on_par = false;

  // If true, d2sigma must be non-null. Activates the d2sigma exact-Hessian
  // term in accumulate_aloq / accumulate_bloq (gated also by use_deriv2_exact
  // and sigma_depends_on_par).
  bool d2sigma_present = false;

  // Per-part Hessian toggles (mirror R opt.hessian flags).
  bool aloq_part1 = true;
  bool aloq_part2 = true;
  bool aloq_part3 = true;
  bool bloq_part1 = true;
  bool bloq_part2 = true;
  bool bloq_part3 = true;
};

// Accumulate ALOQ-row contributions for one condition; accumulators are in/out.
// Every row reduces to grad = A * dwr + B * dw0 + C * dlogs, with dwr, dw0, dlogs
// linear in dpred and dsigma, so the seeds stop that chain rule one step earlier.
void accumulate_aloq_residual(
    int n_obs,                 // ALOQ data rows
    int n_par,                 // parameters spanned by gradient and Hessian
    const double* pred,        // [n_obs]
    const double* dpred,       // [n_obs * n_par], row-major
    const double* d2pred,      // [n_obs * n_par * n_par], row-major, or nullptr
    const double* y_data,      // [n_obs]
    const double* sigma,       // [n_obs]
    const double* dsigma,      // [n_obs * n_par], nullptr unless sigma_depends_on_par
    const double* d2sigma,     // [n_obs * n_par * n_par] or nullptr
    const double* lloq,        // [n_obs], w0 for M4BEAL; nullptr under NONE
    const AccumOpts& opts,
    double& value_acc,
    double& chi2_acc,          // sum of wr^2 alone, censored rows excluded
    double* grad_acc,          // [n_par]
    double* hess_acc,          // [n_par * n_par], column-major
    double* seed_pred = nullptr,   // [n_obs] dvalue/dpred, the reverse-mode seed
    double* seed_sigma = nullptr); // [n_obs] dvalue/dsigma


// Accumulate BLOQ-row (val <= lloq) contributions for one condition, no-op under
// NONE or M1. Arguments as in accumulate_aloq_residual; y_data is LOQ-substituted
// as in res(), `pmax(value, lloq)`.
void accumulate_bloq_residual(
    int n_obs,
    int n_par,
    const double* pred,
    const double* dpred,
    const double* d2pred,
    const double* y_data,
    const double* sigma,
    const double* dsigma,
    const double* d2sigma,
    const double* lloq,
    const AccumOpts& opts,
    double& value_acc,
    double* grad_acc,
    double* hess_acc,
    double* seed_pred = nullptr,
    double* seed_sigma = nullptr);

}  // namespace dmod

#endif  // DMOD_RESIDUAL_KERNEL_H
