// Rcpp-side scaffolding shared by the two exported trust-region kernels,
// trust_kernel.cpp and trustL1_kernel.cpp: argument parsing for bounds and
// parscale, iteration reporting, the per-iteration trace buffers, and the
// objfun callback wrapper.
//
// The numerical core lives in trust_subproblem.h, which is deliberately free of
// Rcpp so it can be reasoned about (and unit-tested) as plain LAPACK code. This
// header is the other half: everything that only exists because the optimiser is
// driven from R.

#ifndef DMOD_TRUST_DRIVER_H
#define DMOD_TRUST_DRIVER_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace dmod { namespace trust_driver {

using namespace Rcpp;

const double kInf = std::numeric_limits<double>::infinity();

// Named entries slot in by name; an unnamed vector broadcasts its first element.
inline void fill_bound(Nullable<NumericVector> nv, std::vector<double>& out,
                       const CharacterVector& parnames, int K) {
  if (nv.isNull()) return;
  NumericVector v(nv);
  if (v.size() == 0) return;
  if (v.hasAttribute("names")) {
    CharacterVector vn = v.names();
    for (int j = 0; j < v.size(); ++j) {
      std::string nm = as<std::string>(vn[j]);
      for (int i = 0; i < K; ++i)
        if (as<std::string>(parnames[i]) == nm) { out[i] = v[j]; break; }
    }
  } else {
    std::fill(out.begin(), out.end(), v[0]);
  }
}

inline void fill_parscale(Nullable<NumericVector> nv, std::vector<double>& ps,
                          int K, const char* who) {
  if (nv.isNull()) return;
  NumericVector v(nv);
  if (v.size() != K) stop("%s: parscale and parinit not same length", who);
  for (int i = 0; i < K; ++i) {
    if (!(v[i] > 0)) stop("%s: parscale not all positive", who);
    if (!std::isfinite(v[i]) || !std::isfinite(1.0 / v[i]))
      stop("%s: parscale or 1/parscale not all finite", who);
    ps[i] = v[i];
  }
}

// Nudge an iterate off any bound it sits on, so the affine scaling stays
// non-degenerate. Mirrors fides' make_non_degenerate.
inline void push_interior(int K, std::vector<double>& z,
                          const std::vector<double>& lb,
                          const std::vector<double>& ub) {
  const double macheps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < K; ++i) {
    const double eps = 100.0 * macheps * std::max(1.0, std::fabs(z[i]));
    const bool fl = std::isfinite(lb[i]), fu = std::isfinite(ub[i]);
    if (fl && fu && ub[i] - lb[i] <= 2.0 * eps) { z[i] = 0.5 * (lb[i] + ub[i]); continue; }
    if (fu && ub[i] - z[i] < eps) z[i] = ub[i] - eps;
    if (fl && z[i] - lb[i] < eps) z[i] = lb[i] + eps;
  }
}

// printIter console output and the traceFile CSV.
struct Reporter {
  bool          printIter = false;
  int           iter_width = 1;
  std::string   path;
  std::ofstream ofs;
  int           K = 0;
  const CharacterVector* parnames = nullptr;

  void init(bool print, int iterlim, Nullable<CharacterVector> traceFile,
            int K_, const CharacterVector& names) {
    printIter  = print;
    iter_width = static_cast<int>(std::to_string(iterlim).size());
    K          = K_;
    parnames   = &names;
    if (traceFile.isNotNull()) {
      CharacterVector tf(traceFile);
      if (tf.size() > 0) path = as<std::string>(tf[0]);
    }
  }
  void operator()(int it, double v, const NumericVector& x, bool head) {
    if (printIter)
      Rcpp::Rcout << "Iteration: " << std::setw(iter_width) << it
                  << "      Objective value: " << v << "\n";
    if (path.empty()) return;
    if (head) {
      ofs.open(path.c_str());
      ofs << "Iteration,Obj";
      for (int i = 0; i < K; ++i) ofs << "," << as<std::string>((*parnames)[i]);
      ofs << "\n";
    }
    ofs << std::setprecision(15) << it << "," << v;
    for (int i = 0; i < K; ++i) ofs << "," << x[i];
    ofs << "\n";
  }
  ~Reporter() { if (ofs.is_open()) ofs.close(); }
};

// Per-iteration trace buffers. Every vector carries exactly one entry per
// completed iteration, which `attach` relies on.
struct Blather {
  std::vector<double> argpath, argtry;
  std::vector<std::string> steptype, stepback;
  std::vector<int>    accept;
  std::vector<double> r, rho, valpath, valtry, preddiff, stepnorm;

  void attach(List& result, int n, int K, const CharacterVector& parnames,
              bool minimize) const {
    NumericMatrix argpath_M(n, K), argtry_M(n, K);
    for (int it = 0; it < n; ++it)
      for (int c = 0; c < K; ++c) {
        argpath_M(it, c) = argpath[(std::size_t) it * K + c];
        argtry_M (it, c) = argtry [(std::size_t) it * K + c];
      }
    argpath_M.attr("dimnames") = List::create(R_NilValue, parnames);
    argtry_M.attr("dimnames")  = List::create(R_NilValue, parnames);

    CharacterVector steptype_out(steptype.size()), stepback_out(stepback.size());
    for (std::size_t i = 0; i < steptype.size(); ++i) steptype_out[i] = steptype[i];
    for (std::size_t i = 0; i < stepback.size(); ++i) stepback_out[i] = stepback[i];

    LogicalVector accept_out(accept.size());
    for (std::size_t i = 0; i < accept.size(); ++i) accept_out[i] = (accept[i] != 0);

    NumericVector preddiff_out(preddiff.size());
    for (std::size_t i = 0; i < preddiff.size(); ++i)
      preddiff_out[i] = minimize ? preddiff[i] : -preddiff[i];

    result["argpath"]  = argpath_M;
    result["argtry"]   = argtry_M;
    result["steptype"] = steptype_out;
    result["stepback"] = stepback_out;
    result["accept"]   = accept_out;
    result["r"]        = NumericVector(r.begin(), r.end());
    result["rho"]      = NumericVector(rho.begin(), rho.end());
    result["valpath"]  = NumericVector(valpath.begin(), valpath.end());
    result["valtry"]   = NumericVector(valtry.begin(), valtry.end());
    result["preddiff"] = preddiff_out;
    result["stepnorm"] = NumericVector(stepnorm.begin(), stepnorm.end());
  }
};

// Consecutive rejected steps with a flat objective before the run counts as
// stagnant. Five rejections cut the radius by ~4^5; a still-flat objective then
// means nothing is left to resolve, typically the ODE solver's noise floor.
const int kStallLimit = 5;

inline const char* subproblem_label(bool is_newton, bool is_hard, bool is_easy) {
  if (is_newton)          return "Newton";
  if (is_hard && is_easy) return "hard-easy";
  if (is_hard)            return "hard-hard";
  return "easy-easy";
}

// Turn an R-level failure into eval_ok = false, but let a user interrupt
// through -- a bare catch(...) would swallow Ctrl-C and count it as a failed
// evaluation.
inline bool eval_objfun(Function& objfun, const NumericVector& x, List& out) {
  try {
    out = as<List>(objfun(x));
  } catch (Rcpp::internal::InterruptedException&) {
    throw;
  } catch (...) {
    return false;
  }
  return true;
}

}}  // namespace dmod::trust_driver

#endif  // DMOD_TRUST_DRIVER_H
