// The modes of the components of a joint FOCEI fit (.fitFoceiJoint()) in lock step: per round one
// evaluation of the data terms of the components still moving, then per component a trust step on
// its data and the normal prior e' P e / 2, accepted by its own ratio; a component whose step finds
// its modes takes its last Newton step without a trial and drops out. The rules are those of the
// inner search of the FOCEI kernel.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

List trust_step_impl(NumericVector x, NumericVector g, NumericMatrix H, double r, NumericVector lower,
                     NumericVector upper, NumericVector parscale, double thetamax);

namespace {

// A component at its etas: the value, gradient and Hessian of its data and prior; `ok` false where
// its prediction failed or its value is not finite.
struct Unit {
  std::vector<double> eta, grad, hess;
  double value = 0.0;
  bool ok = false;
};

// The data terms `d` (list(value, gradient, hessian) or NULL) plus the prior with precision P.
Unit withPrior(SEXP d, const std::vector<double>& eta, const NumericMatrix& P) {
  Unit u;
  u.eta = eta;
  if (Rf_isNull(d)) return u;
  List l(d);
  const int K = eta.size();
  NumericVector g = l["gradient"];
  NumericMatrix H = l["hessian"];
  if (g.size() != K || H.nrow() != K || H.ncol() != K) return u;
  u.value = as<double>(l["value"]);
  u.grad.assign(K, 0.0);
  u.hess.assign((std::size_t) K * K, 0.0);
  double quad = 0.0;
  for (int j = 0; j < K; ++j) {
    double pe = 0.0;
    for (int k = 0; k < K; ++k) pe += P(j, k) * eta[k];
    quad += eta[j] * pe;
    u.grad[j] = g[j] + pe;
    for (int k = 0; k < K; ++k) u.hess[j + (std::size_t) K * k] = H(j, k) + P(j, k);
  }
  u.value += quad / 2;
  u.ok = std::isfinite(u.value);
  return u;
}

// The Newton step H^-1 g by an LU factor, as solve() in R; false where H is singular.
bool newtonStep(const Unit& u, std::vector<double>& s) {
  const int K = u.eta.size(), one = 1;
  std::vector<double> A(u.hess);
  std::vector<int> ipiv(K);
  s = u.grad;
  int info = 0;
  F77_CALL(dgesv)(&K, &one, A.data(), &K, ipiv.data(), s.data(), &K, &info);
  if (info != 0) return false;
  for (double v : s) if (!std::isfinite(v)) return false;
  return true;
}

}  // namespace

// `data_cb(idx, etas)`: per component idx (1-based) at its etas, list(value, gradient, hessian) of its
// data terms in the order of its etas, NULL where its prediction fails. `warm` holds the etas the
// search starts from, `P` the prior precisions, `control` rinit, rmax, iterlim, steptol and mtol.
// Returns the etas at the modes per component.
// [[Rcpp::export]]
List focei_units_lockstep(Function data_cb, List warm, List P, List control) {
  const int U = warm.size();
  const double rinit = as<double>(control["rinit"]), rmax = as<double>(control["rmax"]),
               steptol = as<double>(control["steptol"]), mtol = as<double>(control["mtol"]);
  const int iterlim = as<int>(control["iterlim"]);
  std::vector<NumericMatrix> Pm(U);
  for (int i = 0; i < U; ++i) Pm[i] = as<NumericMatrix>(P[i]);

  // one evaluation of the components `idx` at the etas `etas`
  auto evaluate = [&](const std::vector<int>& idx, const std::vector<std::vector<double>>& etas) {
    IntegerVector which(idx.size());
    List at(idx.size());
    for (std::size_t j = 0; j < idx.size(); ++j) {
      which[j] = idx[j] + 1;
      at[j] = NumericVector(etas[j].begin(), etas[j].end());
    }
    List d = data_cb(which, at);
    std::vector<Unit> out(idx.size());
    for (std::size_t j = 0; j < idx.size(); ++j) out[j] = withPrior(d[j], etas[j], Pm[idx[j]]);
    return out;
  };

  std::vector<int> active(U);
  std::vector<std::vector<double>> w(U);
  for (int i = 0; i < U; ++i) {
    active[i] = i;
    w[i] = as<std::vector<double>>(warm[i]);
  }
  std::vector<Unit> cur = evaluate(active, w);
  for (const Unit& u : cur)
    if (!u.ok) stop("FOCEI: the prediction fails at the modes the inner search starts from.");
  std::vector<double> r(U, rinit);

  for (int it = 1; it <= iterlim && !active.empty(); ++it) {
    std::vector<int> moving;
    std::vector<std::vector<double>> steps, trials;
    std::vector<double> pred, stepnorm;
    std::vector<char> newton;
    for (int i : active) {
      const Unit& u = cur[i];
      const int K = u.eta.size();
      NumericMatrix Hs(K, K);
      for (int j = 0; j < K; ++j)
        for (int k = 0; k < K; ++k)
          Hs(j, k) = 0.5 * (u.hess[j + (std::size_t) K * k] + u.hess[k + (std::size_t) K * j]);
      List st = trust_step_impl(NumericVector(K), NumericVector(u.grad.begin(), u.grad.end()), Hs, r[i],
                                NumericVector(K, R_NegInf), NumericVector(K, R_PosInf), NumericVector(K, 1.0), 0.95);
      const double p = -as<double>(st["quad"]), sn = as<double>(st["stepnorm"]);
      const bool nw = as<bool>(st["newton"]);
      // found where the Newton step is below steptol, or a step on the boundary gains less than
      // mtol of the value: rejections shrank the region below the distance to the modes
      if ((nw && sn <= steptol) || (!nw && p <= mtol * std::max(1.0, std::fabs(u.value)))) {
        std::vector<double> s;
        double ns = 0.0;
        if (newtonStep(u, s)) {
          for (double v : s) ns += v * v;
          if (std::sqrt(ns) <= r[i])
            for (int k = 0; k < K; ++k) w[i][k] -= s[k];
        }
        continue;
      }
      NumericVector sv = st["step"];
      std::vector<double> e(w[i]);
      for (int k = 0; k < K; ++k) e[k] += sv[k];
      moving.push_back(i);
      steps.push_back(as<std::vector<double>>(sv));
      trials.push_back(e);
      pred.push_back(p);
      stepnorm.push_back(sn);
      newton.push_back(nw);
    }
    active = moving;
    if (active.empty()) break;
    std::vector<Unit> trial = evaluate(active, trials);
    for (std::size_t j = 0; j < active.size(); ++j) {
      const int i = active[j];
      const Unit& t = trial[j];
      const double rho = t.ok ? (cur[i].value - t.value) / pred[j] : R_NegInf;
      // near the modes, where a Newton step predicts less than 1e-6 of the value, also a smaller
      // gradient accepts it: the solver's error in the value there exceeds what the step gains
      bool byGradient = false;
      if (rho <= 0 && t.ok && newton[j] && pred[j] <= 1e-6 * std::max(1.0, std::fabs(cur[i].value))) {
        double g0 = 0.0, g1 = 0.0;
        for (std::size_t k = 0; k < t.grad.size(); ++k) {
          g0 += cur[i].grad[k] * cur[i].grad[k];
          g1 += t.grad[k] * t.grad[k];
        }
        byGradient = g1 < g0;
      }
      const bool ok = rho > 0 || byGradient;
      if (!ok || (!byGradient && rho < 0.25)) r[i] /= 4;
      else if (rho > 0.75 && stepnorm[j] >= 0.99 * r[i]) r[i] = std::min(2 * r[i], rmax);
      if (!ok) continue;
      w[i] = trials[j];
      cur[i] = t;
    }
  }
  List out(U);
  for (int i = 0; i < U; ++i) out[i] = NumericVector(w[i].begin(), w[i].end());
  return out;
}
