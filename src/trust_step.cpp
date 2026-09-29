// One trust-region step and one quasi-Newton update, exported so that an
// optimiser driven from R solves its subproblem exactly the way trust() does.
// The multiple-shooting driver is the caller: it condenses its problem onto
// the parameters and then needs nothing but a step on them.
//
// The step is the reflective scheme of trust_kernel.cpp, iteration for
// iteration: parscale, the Coleman-Li scaling, the Moré-Sorensen solve on the
// eigendecomposition, the stepback that keeps the iterate strictly interior.

#include <Rcpp.h>
#include "trust_subproblem.h"
#include "trust_driver.h"
#include "trust_qn.h"
#include <vector>
#include <cmath>
#include <string>

using namespace Rcpp;
using dmod::trust_internal::affine_scaling;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::model_value;
using dmod::trust_internal::stepback;
using dmod::trust_internal::trust_sub;
using dmod::trust_internal::qn_update_B;
using dmod::trust_driver::push_interior;

// [[Rcpp::export]]
List trust_step_impl(NumericVector x, NumericVector g, NumericMatrix H,
                     double r, NumericVector lower, NumericVector upper,
                     NumericVector parscale, double thetamax) {
  const int K = x.size();
  if (g.size() != K || H.nrow() != K || H.ncol() != K ||
      lower.size() != K || upper.size() != K || parscale.size() != K)
    stop("trust_step_impl: dimensions of x, g, H, bounds and parscale differ.");

  std::vector<double> z(K), g_z(K), lbz(K), ubz(K), absv(K), jv(K), sqrtv(K);
  for (int i = 0; i < K; ++i) {
    z[i]   = x[i] * parscale[i];
    g_z[i] = g[i] / parscale[i];
    lbz[i] = lower[i] * parscale[i];
    ubz[i] = upper[i] * parscale[i];
  }
  affine_scaling(K, z.data(), g_z.data(), lbz.data(), ubz.data(),
                 absv.data(), jv.data());
  double opt_measure = 0.0;
  for (int i = 0; i < K; ++i) {
    sqrtv[i] = std::sqrt(absv[i]);
    opt_measure = std::max(opt_measure, std::fabs(absv[i] * g_z[i]));
  }

  std::vector<double> ghat(K), Bhat((std::size_t) K * K);
  for (int i = 0; i < K; ++i) ghat[i] = sqrtv[i] * g_z[i];
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      Bhat[i + (std::size_t) j * K] =
          sqrtv[i] * sqrtv[j] * H(i, j) / (parscale[i] * parscale[j]);
  for (int i = 0; i < K; ++i)
    Bhat[i + (std::size_t) i * K] += std::fabs(g_z[i]) * jv[i];

  std::vector<double> vals(K), vecs((std::size_t) K * K), shat(K);
  eigen_sym_local(Bhat.data(), K, vals.data(), vecs.data());
  bool is_newton = false, is_hard = false, is_easy = false;
  double pred_unused = 0.0;
  trust_sub(K, ghat.data(), vals.data(), vecs.data(), r, shat.data(),
            &pred_unused, &is_newton, &is_hard, &is_easy);

  const double theta_frac =
      std::min(std::max(thetamax, 1.0 - opt_measure), 1.0 - 1e-12);
  std::vector<double> shat_step(K), s_step(K);
  const char* label = "full";
  double m_value = 0.0;
  stepback(K, z.data(), lbz.data(), ubz.data(), sqrtv.data(), ghat.data(),
           Bhat.data(), shat.data(), r, theta_frac, shat_step.data(),
           s_step.data(), &m_value, &label);

  double stepnorm = 0.0;
  for (int i = 0; i < K; ++i) stepnorm += shat_step[i] * shat_step[i];
  stepnorm = std::sqrt(stepnorm);

  std::vector<double> z_try(K);
  for (int i = 0; i < K; ++i) z_try[i] = z[i] + s_step[i];
  push_interior(K, z_try, lbz, ubz);

  // The step as taken, after the push off the bounds, and the plain quadratic
  // model g's + s'Hs/2 along it beside the scaled one trust() scores with.
  NumericVector step(K), xtry(K), vv(K);
  double quad = 0.0;
  for (int i = 0; i < K; ++i) {
    xtry[i] = z_try[i] / parscale[i];
    step[i] = xtry[i] - x[i];
    vv[i] = absv[i];
  }
  for (int i = 0; i < K; ++i) {
    double Hs = 0.0;
    for (int j = 0; j < K; ++j) Hs += H(i, j) * step[j];
    quad += g[i] * step[i] + 0.5 * step[i] * Hs;
  }
  if (x.hasAttribute("names")) {
    step.names() = x.names();
    xtry.names() = x.names();
    vv.names()   = x.names();
  }
  return List::create(
      Named("step")       = step,
      Named("xtry")       = xtry,
      Named("model")      = m_value,
      Named("quad")       = quad,
      Named("stepnorm")   = stepnorm,
      Named("optMeasure") = opt_measure,
      Named("absv")       = vv,
      Named("steptype")   = std::string(label),
      Named("newton")     = is_newton);
}

// [[Rcpp::export]]
List qn_update_impl(NumericMatrix B, NumericVector s, NumericVector y,
                    std::string method) {
  const int K = B.nrow();
  if (B.ncol() != K || s.size() != K || y.size() != K)
    stop("qn_update_impl: dimensions of B, s and y differ.");
  int kind;
  if (method == "bfgs") kind = dmod::trust_internal::QN_BFGS;
  else if (method == "sr1") kind = dmod::trust_internal::QN_SR1;
  else stop("qn_update_impl: method must be \"bfgs\" or \"sr1\".");
  std::vector<double> Bv(B.begin(), B.end()), sv(s.begin(), s.end()),
                      yv(y.begin(), y.end());
  const bool changed = qn_update_B(kind, K, Bv, sv, yv);
  NumericMatrix out(K, K);
  std::copy(Bv.begin(), Bv.end(), out.begin());
  out.attr("dimnames") = B.attr("dimnames");
  return List::create(Named("B") = out, Named("updated") = changed);
}
