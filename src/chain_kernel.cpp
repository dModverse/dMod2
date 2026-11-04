// The native chain of normL2(data, g * x * p, errmodel = e) by the plan of .chainPlan(): per
// condition the stages of the transformation through their own compiled functions, the
// observables and the error model through theirs, and the data term; the solve in between is the
// batch entry of the model, which R calls between chain_params() and chain_data(). Indices are
// 0-based. A function that throws or a value that is not finite returns NULL, and the caller
// takes the R path, which reports it.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "pimpl_kernel.h"
#include "residual_kernel.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#if !defined(_WIN32)
#include <pthread.h>
#endif

using namespace Rcpp;

namespace {

// The raw AD entries of a cppFUN() object, `<model>_eval_ad` and `<model>_eval_ad2`.
using EvalAd = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*);
using EvalAd2 = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, double*, double*,
                         double*, int*, int*, int*, int*, int*);

#if !defined(_WIN32)
bool chainForked = false;
struct ChainForkGuard {
  ChainForkGuard() { pthread_atfork(nullptr, nullptr, [] { chainForked = true; }); }
} chainForkGuard;
#endif

// The requested team, at most the processors, also for fewer conditions, so that libgomp keeps its
// threads and cppDE their scratch; one thread inside a parallel region or in a forked child.
int chainThreads(int requested, int n) {
  int nt = n > 1 ? std::max(1, requested) : 1;
#ifndef _OPENMP
  nt = 1;
#else
  nt = std::min(nt, omp_get_num_procs());
  if (omp_in_parallel()) nt = 1;
#endif
#if !defined(_WIN32)
  if (chainForked) nt = 1;
#endif
  return nt;
}

DL_FUNC address(SEXP x) {
  if (TYPEOF(x) != EXTPTRSXP) return nullptr;
  return R_ExternalPtrAddrFn(x);
}

// The values y [n_obs, n_out], first derivatives dy [n_obs, n_out, nd] and with `d2` second
// derivatives d2y [n_obs, n_out, nd, nd], column-major as the entries write them.
struct Ad {
  std::vector<double> y, dy, d2y;
};

bool evalAd(DL_FUNC f, bool d2, std::vector<double>& x, std::vector<double>& p, double* dX, double* dP,
            double* dX2, double* dP2, int n_obs, int n_vars, int n_params, int n_out, int nd, Ad& o) {
  const std::size_t m = std::max(nd, 1), ny = (std::size_t) n_obs * n_out;
  o.y.assign(ny, 0.0);
  o.dy.assign(ny * m, 0.0);
  if (d2) o.d2y.assign(ny * m * m, 0.0);
  if (x.empty()) x.push_back(0.0);
  if (p.empty()) p.push_back(0.0);
  double none = 0.0;
  try {
    if (d2) {
      int hx = dX2 != nullptr, hp = dP2 != nullptr;
      reinterpret_cast<EvalAd2>(f)(x.data(), p.data(), dX, dP, hx ? dX2 : &none, hp ? dP2 : &none, &hx, &hp,
                                   o.y.data(), o.dy.data(), o.d2y.data(), &n_obs, &n_vars, &n_params, &n_out, &nd);
    } else {
      reinterpret_cast<EvalAd>(f)(x.data(), p.data(), dX, dP, o.y.data(), o.dy.data(), &n_obs, &n_vars, &n_params,
                                  &n_out, &nd);
    }
  } catch (...) {
    return false;
  }
  for (double v : o.y) if (!std::isfinite(v)) return false;
  return true;
}

// The columns `idx` of an array [T, n, nd^k] as an array [T, idx, nd^k], k = 0, 1, 2.
std::vector<double> columns(const double* a, int T, int n, const std::vector<int>& idx, std::size_t blocks) {
  const int m = idx.size();
  std::vector<double> out((std::size_t) T * m * blocks);
  for (std::size_t b = 0; b < blocks; ++b)
    for (int j = 0; j < m; ++j)
      std::memcpy(out.data() + (std::size_t) T * (j + (std::size_t) m * b),
                  a + (std::size_t) T * (idx[j] + (std::size_t) n * b), sizeof(double) * T);
  return out;
}

int position(const List& l, const char* name) {
  CharacterVector nm = l.names();
  for (int i = 0; i < nm.size(); ++i)
    if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) return i;
  return -1;
}

// A stage of the transformation: a compiled function (kind 0) of the entries `in` of the vector
// before it, or a steady state (kind 1) of them, solved by `pimpl` with derivatives of order
// `order`; its `nfo` outputs, then the entries `pass` of the vector before it. `guess` are the
// entries a steady state starts from, -1 where it has none.
struct PStage {
  int kind = 0, nfo = 0, order = 0;
  DL_FUNC fn = nullptr, vjp = nullptr;
  dmod::pimpl::State* pimpl = nullptr;
  std::vector<int> in, pass, guess;
};

// The root of a steady-state stage at its inputs: x [ndep], D [ndep, n_in], D2 [ndep, n_in, n_in].
struct Root {
  const double *x = nullptr, *D = nullptr, *D2 = nullptr;
  int ndep = 0;
};

// A vector with its derivatives along nd directions: v [n], J [n, nd], H [n, nd, nd].
struct Vec {
  int n = 0;
  std::vector<double> v, J, H;
};

// A condition's stages, its outer entries `p_oidx` and per entry its direction `p_odir` and its
// tangent row `p_ovrow`, -1 where it has none.
struct PPlan {
  std::vector<PStage> st;
  std::vector<int> oidx, odir, ovrow;
};

bool readPlan(const List& pc, PPlan& P) {
  P.oidx = as<std::vector<int>>(pc["p_oidx"]);
  P.odir = as<std::vector<int>>(pc["p_odir"]);
  P.ovrow = as<std::vector<int>>(pc["p_ovrow"]);
  List st = pc["stages"];
  P.st.assign(st.size(), PStage());
  for (int s = 0; s < st.size(); ++s) {
    List l = st[s];
    PStage& S = P.st[s];
    S.kind = as<int>(l["kind"]);
    S.nfo = as<int>(l["nfo"]);
    S.in = as<std::vector<int>>(l["in"]);
    S.pass = as<std::vector<int>>(l["pass"]);
    S.guess = as<std::vector<int>>(l["guess"]);
    if (S.kind == 0) {
      S.fn = address(l["fn"]);
      if (S.fn == nullptr) return false;
    } else {
      S.pimpl = dmod::pimpl::stateOf(l["state"]);
      S.order = as<int>(l["order"]);
      if (S.pimpl == nullptr) return false;
    }
    if (!Rf_isNull(l["vjp"])) S.vjp = address(l["vjp"]);
  }
  return true;
}

// The roots of the steady-state stages of a condition, from its entry of `roots` (a list per stage,
// NULL or list(x, D, D2)); the index of the first such stage without one, or -1.
int readRoots(SEXP roots, const PPlan& P, std::vector<Root>& R) {
  R.assign(P.st.size(), Root());
  for (std::size_t s = 0; s < P.st.size(); ++s) {
    if (P.st[s].kind != 1) continue;
    SEXP r = (!Rf_isNull(roots) && (R_xlen_t) s < Rf_xlength(roots)) ? VECTOR_ELT(roots, s) : R_NilValue;
    if (Rf_isNull(r)) return s;
    List l(r);
    NumericVector x = l["x"];
    R[s].x = REAL(x);
    R[s].ndep = x.size();
    SEXP D = l["D"], D2 = l["D2"];
    if (!Rf_isNull(D)) R[s].D = REAL(D);
    if (!Rf_isNull(D2)) R[s].D2 = REAL(D2);
  }
  return -1;
}

// The outer entries of a condition with their derivatives: unit directions, or with `V` the rows
// of the tangents.
Vec outerVec(const PPlan& P, const NumericVector& v, int nd, bool deriv2, const NumericMatrix* V) {
  Vec u;
  u.n = P.oidx.size();
  u.v.resize(u.n);
  u.J.assign((std::size_t) u.n * nd, 0.0);
  if (deriv2) u.H.assign((std::size_t) u.n * nd * nd, 0.0);
  for (int j = 0; j < u.n; ++j) {
    u.v[j] = v[P.oidx[j]];
    if (V != nullptr) {
      for (int d = 0; d < nd && P.ovrow[j] >= 0; ++d) u.J[j + (std::size_t) u.n * d] = (*V)(P.ovrow[j], d);
    } else if (nd > 0 && P.odir[j] >= 0) {
      u.J[j + (std::size_t) u.n * P.odir[j]] = 1.0;
    }
  }
  return u;
}

// The stages [0, upto) of a condition on `u`, which becomes their result; with `keep` the input
// of every stage. False where a function fails or a root is missing.
bool runStages(const PPlan& P, const std::vector<Root>& R, Vec& u, int upto, int nd, bool deriv2,
               std::vector<Vec>* keep) {
  for (int s = 0; s < upto; ++s) {
    const PStage& S = P.st[s];
    const int nin = S.in.size(), np = S.pass.size(), nfo = S.nfo;
    const std::size_t m = std::max(nd, 1), un = u.n;
    Vec o;
    o.n = nfo + np;
    o.v.assign(o.n, 0.0);
    o.J.assign((std::size_t) o.n * nd, 0.0);
    if (deriv2) o.H.assign((std::size_t) o.n * nd * nd, 0.0);
    const std::size_t on = o.n, nd2 = (std::size_t) nd * nd;
    if (S.kind == 0) {
      std::vector<double> pin(nin), dP((std::size_t) std::max(nin, 1) * m, 0.0), dP2;
      for (int j = 0; j < nin; ++j) {
        pin[j] = u.v[S.in[j]];
        for (int d = 0; d < nd; ++d) dP[j + (std::size_t) nin * d] = u.J[S.in[j] + un * d];
      }
      // the input of the first stage has no second derivatives
      const bool h2 = deriv2 && nd > 0 && s > 0;
      if (h2) {
        dP2.assign((std::size_t) std::max(nin, 1) * nd2, 0.0);
        for (std::size_t b = 0; b < nd2; ++b)
          for (int j = 0; j < nin; ++j) dP2[j + (std::size_t) nin * b] = u.H[S.in[j] + un * b];
      }
      std::vector<double> none;
      Ad r;
      if (!evalAd(S.fn, deriv2, none, pin, nullptr, nd > 0 ? dP.data() : nullptr, nullptr, h2 ? dP2.data() : nullptr,
                  1, 0, nin, nfo, nd, r))
        return false;
      for (int i = 0; i < nfo; ++i) {
        o.v[i] = r.y[i];
        for (int d = 0; d < nd; ++d) o.J[i + on * d] = r.dy[i + (std::size_t) nfo * d];
        if (deriv2)
          for (std::size_t b = 0; b < nd2; ++b) o.H[i + on * b] = r.d2y[i + (std::size_t) nfo * b];
      }
    } else {
      // the chain rule through the implicit function: D J, and D H + J' D2 J per root entry
      const Root& Q = R[s];
      if (Q.x == nullptr || (nd > 0 && Q.D == nullptr) || (deriv2 && nd > 0 && Q.D2 == nullptr)) return false;
      const int nq = Q.ndep;
      for (int i = 0; i < nq; ++i) o.v[i] = Q.x[i];
      // an input at -1 is a resting state absent from the vector: 0, without derivatives
      for (int d = 0; d < nd; ++d)
        for (int j = 0; j < nin; ++j) {
          if (S.in[j] < 0) continue;
          const double Jj = u.J[S.in[j] + un * d];
          if (Jj != 0.0)
            for (int i = 0; i < nq; ++i) o.J[i + on * d] += Q.D[i + (std::size_t) nq * j] * Jj;
        }
      if (deriv2)
        for (int b = 0; b < nd; ++b)
          for (int a = 0; a < nd; ++a) {
            const std::size_t ab = a + (std::size_t) nd * b;
            for (int j = 0; j < nin; ++j) {
              if (S.in[j] < 0) continue;
              const double Hj = u.H[S.in[j] + un * ab], Ja = u.J[S.in[j] + un * a];
              for (int i = 0; i < nq; ++i) o.H[i + on * ab] += Q.D[i + (std::size_t) nq * j] * Hj;
              if (Ja == 0.0) continue;
              for (int l = 0; l < nin; ++l) {
                if (S.in[l] < 0) continue;
                const double Jb = u.J[S.in[l] + un * b];
                if (Jb == 0.0) continue;
                for (int i = 0; i < nq; ++i)
                  o.H[i + on * ab] += Q.D2[i + (std::size_t) nq * (j + (std::size_t) nin * l)] * Ja * Jb;
              }
            }
          }
    }
    for (int q = 0; q < np; ++q) {
      const std::size_t i = nfo + q, k = S.pass[q];
      o.v[i] = u.v[k];
      for (int d = 0; d < nd; ++d) o.J[i + on * d] = u.J[k + un * d];
      if (deriv2)
        for (std::size_t b = 0; b < nd2; ++b) o.H[i + on * b] = u.H[k + un * b];
    }
    if (keep != nullptr) keep->push_back(std::move(u));
    u = std::move(o);
  }
  return true;
}

}  // namespace


// The transformations of all conditions and the argument lists of the batch solve. `values` holds
// the outer values per condition, or one vector for all; `tangents`, a matrix over the free names
// of a sweep plan, the directions instead of the free names themselves. The steady-state stages
// are solved on the way, all conditions at a stage at once. Returns list(pv, pj, ph, args, roots,
// warnings): per condition the inner values, their first [n_out, nd] and second derivatives
// [n_out, nd, nd], the solve arguments and per stage its root (list(x, D, D2), NULL for a
// compiled stage). NULL where a function fails or a root has no derivatives; an error where no
// root is found.
// [[Rcpp::export]]
SEXP chain_params(List plan, List values, bool deriv, bool deriv2, int threads, SEXP tangents) {
  List conds = plan["cond"];
  const int K = conds.size();
  const bool shared = values.size() == 1;
  const bool sweep = !Rf_isNull(tangents);
  NumericMatrix V = sweep ? NumericMatrix(tangents) : NumericMatrix(0, 0);
  std::vector<PPlan> pp(K);
  std::vector<std::vector<Root>> rt(K);
  std::vector<std::vector<dmod::pimpl::Root>> solved(K);
  std::vector<Vec> base(K), res(K);
  std::vector<int> nd(K), need(K, -1);
  for (int c = 0; c < K; ++c) {
    List pc = conds[c];
    if (!readPlan(pc, pp[c])) return R_NilValue;
    nd[c] = deriv ? (sweep ? V.ncol() : as<int>(pc["nd"])) : 0;
    rt[c].assign(pp[c].st.size(), Root());
    solved[c].assign(pp[c].st.size(), dmod::pimpl::Root());
    NumericVector v = values[shared ? 0 : c];
    base[c] = outerVec(pp[c], v, nd[c], deriv2, sweep ? &V : nullptr);
  }

  std::vector<char> ok(K, 1);
  std::vector<std::string> warnings;
  const int nt = chainThreads(threads, K);
  for (;;) {
    for (int c = 0; c < K; ++c) {
      need[c] = -1;
      for (std::size_t s = 0; s < pp[c].st.size() && need[c] < 0; ++s)
        if (pp[c].st[s].kind == 1 && rt[c][s].x == nullptr) need[c] = s;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
    for (int c = 0; c < K; ++c) {
      res[c] = base[c];
      const int upto = need[c] >= 0 ? need[c] : (int) pp[c].st.size();
      ok[c] = runStages(pp[c], rt[c], res[c], upto, nd[c], deriv2, nullptr);
    }
    for (int c = 0; c < K; ++c) if (!ok[c]) return R_NilValue;
    std::vector<dmod::pimpl::Request> req;
    std::vector<int> at;
    for (int c = 0; c < K; ++c) {
      if (need[c] < 0) continue;
      const PStage& S = pp[c].st[need[c]];
      dmod::pimpl::Request q;
      q.st = S.pimpl;
      q.order = S.order;
      for (int j : S.in) q.pv.push_back(j >= 0 ? res[c].v[j] : 0.0);
      for (int j : S.guess) q.guess.push_back(j >= 0 ? res[c].v[j] : 1.0);
      req.push_back(std::move(q));
      at.push_back(c);
    }
    if (req.empty()) break;
    std::vector<dmod::pimpl::Root> got;
    std::string error;
    if (!dmod::pimpl::roots(req, got, nt, error, warnings)) stop(error);
    for (std::size_t j = 0; j < at.size(); ++j) {
      const int c = at[j], s = need[c];
      if (req[j].order >= 1 && (!got[j].hasD || (req[j].order >= 2 && !got[j].hasD2))) return R_NilValue;
      solved[c][s] = std::move(got[j]);
      const dmod::pimpl::Root& q = solved[c][s];
      rt[c][s].x = q.x.data();
      rt[c][s].ndep = q.x.size();
      rt[c][s].D = q.hasD ? q.D.data() : nullptr;
      rt[c][s].D2 = q.hasD2 ? q.D2.data() : nullptr;
    }
  }

  List roots(K);
  for (int c = 0; c < K; ++c) {
    List per(pp[c].st.size());
    for (std::size_t s = 0; s < pp[c].st.size(); ++s) {
      if (pp[c].st[s].kind != 1) continue;
      const dmod::pimpl::Root& q = solved[c][s];
      per[s] = List::create(Named("x") = NumericVector(q.x.begin(), q.x.end()),
                            Named("D") = q.hasD ? (SEXP) NumericVector(q.D.begin(), q.D.end()) : R_NilValue,
                            Named("D2") = q.hasD2 ? (SEXP) NumericVector(q.D2.begin(), q.D2.end()) : R_NilValue);
    }
    roots[c] = per;
  }

  List outV(K), outJ(K), outH(K), args(K);
  for (int c = 0; c < K; ++c) {
    List pc = conds[c];
    const std::size_t no = res[c].n, n = nd[c];
    const Vec& u = res[c];
    outV[c] = NumericVector(u.v.begin(), u.v.end());
    outJ[c] = NumericVector(u.J.begin(), u.J.end());
    if (deriv2) outH[c] = NumericVector(u.H.begin(), u.H.end());
    // the prepared arguments of the solve with the numbers of this evaluation
    SEXP prepared = pc["x_args"];
    List a = Rf_shallow_duplicate(prepared);
    IntegerVector xp = pc["x_parm"];
    const int atParms = position(a, "parms"), atTangent = position(a, "tangent"), atHessian = position(a, "hessian");
    if (atParms < 0 || (deriv && atTangent < 0) || (deriv2 && atHessian < 0)) return R_NilValue;
    const int nreq = xp.size();
    NumericVector parms(nreq);
    for (int q = 0; q < nreq; ++q) parms[q] = u.v[xp[q]];
    a[atParms] = parms;
    if (deriv) {
      NumericMatrix tg(nreq, (int) n);
      for (std::size_t d = 0; d < n; ++d)
        for (int q = 0; q < nreq; ++q) tg[q + nreq * d] = u.J[xp[q] + no * d];
      RObject old = a[atTangent];
      tg.attr("dimnames") = old.attr("dimnames");
      a[atTangent] = tg;
    }
    if (deriv2) {
      RObject old = a[atHessian];
      NumericVector hs((std::size_t) nreq * n * n);
      for (std::size_t b = 0; b < n * n; ++b)
        for (int q = 0; q < nreq; ++q) hs[q + nreq * b] = u.H[xp[q] + no * b];
      hs.attr("dim") = old.attr("dim");
      hs.attr("dimnames") = old.attr("dimnames");
      a[atHessian] = hs;
    }
    args[c] = a;
  }
  return List::create(Named("pv") = outV, Named("pj") = outJ, Named("ph") = outH, Named("args") = args,
                      Named("roots") = roots, Named("warnings") = wrap(warnings));
}


namespace {

// What the parallel region reads of a condition, gathered where the R API is allowed.
struct In {
  DL_FUNC yf = nullptr, ef = nullptr;
  const double *tsol = nullptr, *xv = nullptr, *xs = nullptr, *xh = nullptr, *pv = nullptr, *pj = nullptr,
               *ph = nullptr;
  int ntime = 0, nstate = 0, nd = 0, npout = 0, nyo = 0, neo = 0;
  std::vector<int> yvar, ypar, evar, epar, t, o, e, sna;
  std::vector<double> sfix, ydat, lloq;
  std::vector<char> bloq;
};

// The inputs of every condition from the plan, chain_params() and the batch solve; with `data`
// the data rows. False where a solve failed or the plan does not fit the results.
bool gather(const List& conds, const List& pre, const List& raw, bool deriv, bool deriv2, bool data,
            std::vector<In>& in) {
  const int K = conds.size();
  List pvL = pre["pv"], pjL = pre["pj"], phL = pre["ph"];
  in.assign(K, In());
  for (int c = 0; c < K; ++c) {
    List pc = conds[c];
    In& I = in[c];
    List r = raw[c];
    const int dg = position(r, "diagnostics");
    if (dg >= 0) {
      List d = r[dg];
      const int rc = position(d, "return_code");
      if (rc >= 0 && as<int>(d[rc]) != 0) return false;
    }
    const int atVar = position(r, "variable"), atTan = position(r, "tangent"), atHes = position(r, "hessian"),
              atTime = position(r, "time");
    if (atVar < 0 || atTime < 0) return false;
    NumericMatrix xv = r[atVar];
    NumericVector tsol = r[atTime];
    I.xv = REAL(xv);
    I.tsol = REAL(tsol);
    I.ntime = xv.nrow();
    I.nstate = xv.ncol();
    if (tsol.size() != I.ntime) return false;
    I.nd = deriv ? as<int>(pc["nd"]) : 0;
    if (I.nd > 0) {
      if (atTan < 0 || (deriv2 && atHes < 0)) return false;
      NumericVector xs = r[atTan];
      I.xs = REAL(xs);
      if (deriv2) {
        NumericVector xh = r[atHes];
        I.xh = REAL(xh);
      }
    }
    NumericVector pv = pvL[c], pj = pjL[c];
    I.pv = REAL(pv);
    I.pj = pj.size() ? REAL(pj) : nullptr;
    if (deriv2 && I.nd > 0) {
      NumericVector ph = phL[c];
      I.ph = REAL(ph);
    }
    I.npout = pv.size();
    I.yf = address(pc["y_fn"]);
    if (I.yf == nullptr) return false;
    I.yvar = as<std::vector<int>>(pc["y_var"]);
    I.ypar = as<std::vector<int>>(pc["y_par"]);
    I.nyo = as<int>(pc["n_yout"]);
    if (!Rf_isNull(pc["e_fn"])) {
      I.ef = address(pc["e_fn"]);
      if (I.ef == nullptr) return false;
      I.evar = as<std::vector<int>>(pc["e_var"]);
      I.epar = as<std::vector<int>>(pc["e_par"]);
      I.neo = as<int>(pc["n_eout"]);
    }
    if (!data) continue;
    // the rows of the data times in the times the solve returns, which may start before them
    NumericVector tdat = pc["t_data"];
    I.t.resize(tdat.size());
    for (int i = 0; i < tdat.size(); ++i) {
      const double* hit = std::lower_bound(tsol.begin(), tsol.end(), tdat[i]);
      if (hit == tsol.end() || *hit != tdat[i]) return false;
      I.t[i] = hit - tsol.begin();
    }
    I.o = as<std::vector<int>>(pc["o_idx"]);
    I.e = as<std::vector<int>>(pc["e_idx"]);
    I.sna = as<std::vector<int>>(pc["sigma_na"]);
    I.sfix = as<std::vector<double>>(pc["sigma_fixed"]);
    I.ydat = as<std::vector<double>>(pc["y_data"]);
    I.lloq = as<std::vector<double>>(pc["lloq"]);
    IntegerVector b = pc["bloq"];
    I.bloq.assign(b.begin(), b.end());
  }
  return true;
}

// The observables of a condition on the grid of its solve and, with `withErr`, the error model on
// them; inputs from the states and the transformation.
bool observe(const In& I, bool deriv2, bool withErr, Ad& Y, Ad& E) {
  const int T = I.ntime, nd = I.nd;
  const std::size_t nd2 = (std::size_t) nd * nd;
  auto params = [&](const std::vector<int>& idx, std::vector<double>& p, std::vector<double>& dP,
                    std::vector<double>& dP2) {
    const int m = idx.size();
    p.resize(m);
    for (int j = 0; j < m; ++j) p[j] = I.pv[idx[j]];
    if (nd == 0) return;
    dP.assign((std::size_t) std::max(m, 1) * nd, 0.0);
    for (int d = 0; d < nd; ++d)
      for (int j = 0; j < m; ++j) dP[j + (std::size_t) m * d] = I.pj[idx[j] + (std::size_t) I.npout * d];
    if (!deriv2) return;
    dP2.assign((std::size_t) std::max(m, 1) * nd2, 0.0);
    for (std::size_t b = 0; b < nd2; ++b)
      for (int j = 0; j < m; ++j) dP2[j + m * b] = I.ph[idx[j] + I.npout * b];
  };
  const bool d2 = deriv2 && nd > 0;
  std::vector<double> p, dP, dP2;
  params(I.ypar, p, dP, dP2);
  std::vector<double> x = columns(I.xv, T, I.nstate, I.yvar, 1), dX, dX2;
  if (nd > 0) dX = columns(I.xs, T, I.nstate, I.yvar, nd);
  if (d2) dX2 = columns(I.xh, T, I.nstate, I.yvar, nd2);
  if (!evalAd(I.yf, deriv2, x, p, nd > 0 ? dX.data() : nullptr, nd > 0 ? dP.data() : nullptr,
              d2 ? dX2.data() : nullptr, d2 ? dP2.data() : nullptr, T, (int) I.yvar.size(), (int) I.ypar.size(),
              I.nyo, nd, Y))
    return false;
  if (!withErr) return true;
  if (I.ef == nullptr) return false;
  std::vector<double> ep, edP, edP2;
  params(I.epar, ep, edP, edP2);
  std::vector<double> ex = columns(Y.y.data(), T, I.nyo, I.evar, 1), edX, edX2;
  if (nd > 0) edX = columns(Y.dy.data(), T, I.nyo, I.evar, nd);
  if (d2) edX2 = columns(Y.d2y.data(), T, I.nyo, I.evar, nd2);
  return evalAd(I.ef, deriv2, ex, ep, nd > 0 ? edX.data() : nullptr, nd > 0 ? edP.data() : nullptr,
                d2 ? edX2.data() : nullptr, d2 ? edP2.data() : nullptr, T, (int) I.evar.size(),
                (int) I.epar.size(), I.neo, nd, E);
}

dmod::BloqMode bloqMode(const std::string& s) {
  return s == "M1" ? dmod::BloqMode::M1 : s == "M3" ? dmod::BloqMode::M3 : s == "M4NM" ? dmod::BloqMode::M4NM
       : s == "M4BEAL" ? dmod::BloqMode::M4BEAL : dmod::BloqMode::NONE;
}

// A matrix [T, 1 + n] with the times first and the values `v` [T, n] behind them, named
// c("time", names), with the derivatives `d` [T, n, nd] and `d2` [T, n, nd, nd] as attributes.
SEXP frame(const double* t, const Ad& a, int T, int n, int nd, bool deriv2, const CharacterVector& names,
           const CharacterVector& dirs) {
  NumericMatrix m(T, n + 1);
  std::memcpy(&m(0, 0), t, sizeof(double) * T);
  if (n > 0) std::memcpy(&m(0, 1), a.y.data(), sizeof(double) * T * n);
  CharacterVector cn(n + 1);
  cn[0] = "time";
  for (int j = 0; j < n; ++j) cn[j + 1] = names[j];
  m.attr("dimnames") = List::create(R_NilValue, cn);
  NumericVector d1(a.dy.begin(), a.dy.begin() + (std::size_t) T * n * nd);
  d1.attr("dim") = IntegerVector::create(T, n, nd);
  d1.attr("dimnames") = List::create(R_NilValue, names, dirs);
  m.attr("deriv") = d1;
  if (deriv2) {
    NumericVector d(a.d2y.begin(), a.d2y.begin() + (std::size_t) T * n * nd * nd);
    d.attr("dim") = IntegerVector::create(T, n, nd, nd);
    d.attr("dimnames") = List::create(R_NilValue, names, dirs, dirs);
    m.attr("deriv2") = d;
  }
  return m;
}

}  // namespace


// Observables, error model and data term of all conditions from the transformations of
// chain_params() and the results of the batch solve. Returns list(value, chi2, gradient, hessian,
// local): over the `n_global` parameters of the plan, or per condition in its directions; with
// `want_seed` also seed_pred and seed_sigma, the derivatives of the value in the prediction and
// sigma of every data row, in the order of the conditions and of their data. NULL where a solve or
// a function fails or a prediction at a data point is not finite.
// [[Rcpp::export]]
SEXP chain_data(List plan, List pre, List raw, bool deriv, bool deriv2, bool build_hessian, bool local,
                int threads, std::string bloq_mode, double t_df, bool want_seed) {
  List conds = plan["cond"];
  const int K = conds.size();
  const int nGlobal = as<int>(plan["n_global"]);
  std::vector<In> in;
  if (!gather(conds, pre, raw, deriv, deriv2, true, in)) return R_NilValue;
  std::vector<int> rowOff(K + 1, 0);
  for (int c = 0; c < K; ++c) rowOff[c + 1] = rowOff[c] + (int) in[c].t.size();
  NumericVector seedP(want_seed ? rowOff[K] : 0), seedS(want_seed ? rowOff[K] : 0);
  double *SPd = REAL(seedP), *SSd = REAL(seedS);

  dmod::AccumOpts base;
  base.build_hessian = build_hessian;
  base.use_deriv2_exact = deriv2;
  base.bloq_mode = bloqMode(bloq_mode);
  if (t_df > 0) base.t_df = t_df;

  std::vector<char> ok(K, 1);
  std::vector<double> value(K, 0.0), chi2(K, 0.0);
  std::vector<std::vector<double>> grad(K), hess(K);
  const int nt = chainThreads(threads, K);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int c = 0; c < K; ++c) {
    const In& I = in[c];
    const int T = I.ntime, nd = I.nd;
    const std::size_t nd2 = (std::size_t) nd * nd;
    bool anyNa = false;
    for (int v : I.sna) anyNa = anyNa || v;
    Ad Y, E;
    if (!observe(I, deriv2, anyNa, Y, E)) {
      ok[c] = 0;
      continue;
    }

    // the data rows, those above the limit first
    const int n = I.t.size();
    std::vector<int> perm(n);
    int na = 0;
    for (int i = 0; i < n; ++i) if (!I.bloq[i]) perm[na++] = i;
    int nb = na;
    for (int i = 0; i < n; ++i) if (I.bloq[i]) perm[nb++] = i;
    const bool d2 = deriv2 && nd > 0;
    std::vector<double> pred(n), sig(n), ydat(n), lloq(n), dpred((std::size_t) n * nd, 0.0), dsig, d2pred, d2sig;
    if (anyNa && nd > 0) dsig.assign((std::size_t) n * nd, 0.0);
    if (d2) d2pred.assign((std::size_t) n * nd2, 0.0);
    if (d2 && anyNa) d2sig.assign((std::size_t) n * nd2, 0.0);
    bool hasDs = false;
    for (int j = 0; j < n; ++j) {
      const int i = perm[j], ti = I.t[i], oi = I.o[i];
      pred[j] = Y.y[ti + (std::size_t) T * oi];
      ydat[j] = I.ydat[i];
      lloq[j] = I.lloq[i];
      for (int d = 0; d < nd; ++d) dpred[(std::size_t) j * nd + d] = Y.dy[ti + (std::size_t) T * (oi + (std::size_t) I.nyo * d)];
      if (d2)
        for (int a = 0; a < nd; ++a)
          for (int b = 0; b < nd; ++b)
            d2pred[j * nd2 + (std::size_t) a * nd + b] =
                Y.d2y[ti + (std::size_t) T * (oi + (std::size_t) I.nyo * (a + (std::size_t) nd * b))];
      if (!I.sna[i]) {
        sig[j] = I.sfix[i];
        continue;
      }
      const int ei = I.e[i];
      sig[j] = E.y[ti + (std::size_t) T * ei];
      for (int d = 0; d < nd; ++d) {
        dsig[(std::size_t) j * nd + d] = E.dy[ti + (std::size_t) T * (ei + (std::size_t) I.neo * d)];
        hasDs = true;
      }
      if (d2)
        for (int a = 0; a < nd; ++a)
          for (int b = 0; b < nd; ++b)
            d2sig[j * nd2 + (std::size_t) a * nd + b] =
                E.d2y[ti + (std::size_t) T * (ei + (std::size_t) I.neo * (a + (std::size_t) nd * b))];
    }
    for (int j = 0; j < n; ++j) if (!std::isfinite(pred[j]) || !std::isfinite(sig[j])) ok[c] = 0;
    if (!ok[c]) continue;

    grad[c].assign(nd, 0.0);
    if (build_hessian) hess[c].assign(nd2, 0.0);
    dmod::AccumOpts opts = base;
    opts.sigma_depends_on_par = hasDs;
    opts.d2sigma_present = d2 && hasDs;
    double v = 0.0, ch = 0.0;
    std::vector<double> sp(want_seed ? n : 0, 0.0), ss(want_seed ? n : 0, 0.0);
    if (na > 0)
      dmod::accumulate_aloq_residual(na, nd, pred.data(), dpred.data(), d2 ? d2pred.data() : nullptr, ydat.data(),
                                     sig.data(), hasDs ? dsig.data() : nullptr,
                                     d2 && hasDs ? d2sig.data() : nullptr, lloq.data(), opts, v, ch,
                                     grad[c].data(), build_hessian ? hess[c].data() : nullptr,
                                     want_seed ? sp.data() : nullptr, want_seed ? ss.data() : nullptr);
    if (n > na) {
      const std::size_t off = na, off2 = off * nd, off3 = off * nd2;
      dmod::accumulate_bloq_residual(n - na, nd, pred.data() + off, dpred.data() + off2,
                                     d2 ? d2pred.data() + off3 : nullptr, ydat.data() + off, sig.data() + off,
                                     hasDs ? dsig.data() + off2 : nullptr,
                                     d2 && hasDs ? d2sig.data() + off3 : nullptr, lloq.data() + off, opts, v,
                                     grad[c].data(), build_hessian ? hess[c].data() : nullptr,
                                     want_seed ? sp.data() + off : nullptr, want_seed ? ss.data() + off : nullptr);
    }
    // the seeds back in the order of the data
    for (int j = 0; j < n && want_seed; ++j) {
      SPd[rowOff[c] + perm[j]] = sp[j];
      SSd[rowOff[c] + perm[j]] = ss[j];
    }
    value[c] = v;
    chi2[c] = ch;
  }
  for (int c = 0; c < K; ++c) if (!ok[c]) return R_NilValue;

  // in the order of the conditions, so that the sum does not depend on the threads
  double vTot = 0.0, cTot = 0.0;
  for (int c = 0; c < K; ++c) { vTot += value[c]; cTot += chi2[c]; }
  List out = List::create(Named("value") = vTot, Named("chi2") = cTot);
  if (want_seed) {
    out["seed_pred"] = seedP;
    out["seed_sigma"] = seedS;
  }
  if (local) {
    List loc(K);
    for (int c = 0; c < K; ++c) {
      List pc = conds[c];
      CharacterVector nm = pc["dir_names"];
      const int nd = in[c].nd;
      NumericVector g(grad[c].begin(), grad[c].end());
      if (nd > 0) g.names() = nm;
      RObject h = R_NilValue;
      if (build_hessian) {
        NumericMatrix H(nd, nd);
        if (nd > 0) {
          std::memcpy(&H(0, 0), hess[c].data(), sizeof(double) * nd * nd);
          H.attr("dimnames") = List::create(nm, nm);
        }
        h = H;
      }
      loc[c] = List::create(Named("value") = value[c], Named("chi2") = chi2[c], Named("gradient") = g,
                            Named("hessian") = h);
    }
    out["local"] = loc;
    return out;
  }
  std::vector<double> g(nGlobal, 0.0), H;
  if (build_hessian) H.assign((std::size_t) nGlobal * nGlobal, 0.0);
  for (int c = 0; c < K; ++c) {
    List pc = conds[c];
    IntegerVector dg = pc["dir_global"];
    const int nd = in[c].nd;
    for (int a = 0; a < nd; ++a) {
      const int ga = dg[a];
      if (ga < 0) continue;
      g[ga] += grad[c][a];
      if (!build_hessian) continue;
      for (int b = 0; b < nd; ++b) {
        const int gb = dg[b];
        if (gb >= 0) H[ga + (std::size_t) nGlobal * gb] += hess[c][a + (std::size_t) nd * b];
      }
    }
  }
  out["gradient"] = NumericVector(g.begin(), g.end());
  if (build_hessian) {
    NumericMatrix Hm(nGlobal, nGlobal);
    if (nGlobal > 0) std::memcpy(&Hm(0, 0), H.data(), sizeof(double) * H.size());
    out["hessian"] = Hm;
  }
  return out;
}


// The prediction and error model of all conditions from the transformations of chain_params()
// and the results of the batch solve: per condition list(pred, err), matrices with the times
// first and the derivatives as attributes, `err` NULL without an error model. NULL where a solve
// or a function fails.
// [[Rcpp::export]]
SEXP chain_frames(List plan, List pre, List raw, bool deriv2, int threads) {
  List conds = plan["cond"];
  const int K = conds.size();
  std::vector<In> in;
  if (!gather(conds, pre, raw, true, deriv2, false, in)) return R_NilValue;
  std::vector<char> ok(K, 1);
  std::vector<Ad> Y(K), E(K);
  const int nt = chainThreads(threads, K);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int c = 0; c < K; ++c) ok[c] = observe(in[c], deriv2, in[c].ef != nullptr, Y[c], E[c]);
  for (int c = 0; c < K; ++c) if (!ok[c]) return R_NilValue;

  List out(K);
  for (int c = 0; c < K; ++c) {
    List pc = conds[c];
    const In& I = in[c];
    CharacterVector dirs = pc["dir_names"], yn = pc["y_names"];
    SEXP pred = PROTECT(frame(I.tsol, Y[c], I.ntime, I.nyo, I.nd, deriv2 && I.nd > 0, yn, dirs));
    SEXP err = R_NilValue;
    if (I.ef != nullptr) {
      CharacterVector en = pc["e_names"];
      err = frame(I.tsol, E[c], I.ntime, I.neo, I.nd, deriv2 && I.nd > 0, en, dirs);
    }
    PROTECT(err);
    out[c] = List::create(Named("pred") = pred, Named("err") = err);
    UNPROTECT(2);
  }
  return out;
}


// ---- Forward over reverse ------------------------------------------------------------------
// The sweeps of a plan with tangents (.chainPlan(tangents = )): chain_rows() gives the data rows
// along the tangents, chain_back() pulls seeds on them through the error model and the
// observables onto the states and fills the arguments of the backward solve, chain_pull() takes
// the answer of that solve through the transformation onto the free names. Seeds come in blocks:
// a seed and its derivatives along the K tangents.

namespace {

// The raw forward over reverse entry of a cppFUN() object, `<model>_vjp_ad`.
using VjpAd = void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                       double*, int*, int*, int*, int*, int*, int*);

// Cotangents of the inputs wx [n_obs, n_vars, S], wp [n_params, S] and their derivatives along
// the tangents dwx [n_obs, n_vars, S, K], dwp [n_params, S, K].
struct Pull {
  std::vector<double> wx, wp, dwx, dwp;
};

bool vjpAd(DL_FUNC f, std::vector<double>& x, std::vector<double>& p, double* vx, double* vp, double* w,
           double* dw, int n_obs, int n_vars, int n_params, int n_out, int S, int K, Pull& o) {
  const std::size_t s1 = std::max(S, 1), k1 = std::max(K, 1);
  std::vector<double> y((std::size_t) std::max(n_obs * n_out, 1));
  o.wx.assign((std::size_t) n_obs * n_vars * s1 + 1, 0.0);
  o.wp.assign((std::size_t) n_params * s1 + 1, 0.0);
  o.dwx.assign((std::size_t) n_obs * n_vars * s1 * k1 + 1, 0.0);
  o.dwp.assign((std::size_t) n_params * s1 * k1 + 1, 0.0);
  if (x.empty()) x.push_back(0.0);
  if (p.empty()) p.push_back(0.0);
  try {
    reinterpret_cast<VjpAd>(f)(x.data(), p.data(), w, vx, vp, dw, y.data(), o.wx.data(), o.wp.data(), o.dwx.data(),
                               o.dwp.data(), &n_obs, &n_vars, &n_params, &n_out, &S, &K);
  } catch (...) {
    return false;
  }
  for (const auto* v : {&o.wx, &o.wp, &o.dwx, &o.dwp})
    for (double a : *v) if (!std::isfinite(a)) return false;
  return true;
}

// The raw first-order entry `<model>_vjp`: the cotangents of the inputs wx, wp alone.
using Vjp1 = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*);

bool vjp1(DL_FUNC f, std::vector<double>& x, std::vector<double>& p, double* w, int n_obs, int n_vars, int n_params,
          int n_out, int S, Pull& o) {
  const std::size_t s1 = std::max(S, 1);
  std::vector<double> y((std::size_t) std::max(n_obs * n_out, 1));
  o.wx.assign((std::size_t) n_obs * n_vars * s1 + 1, 0.0);
  o.wp.assign((std::size_t) n_params * s1 + 1, 0.0);
  o.dwx.clear();
  o.dwp.clear();
  if (x.empty()) x.push_back(0.0);
  if (p.empty()) p.push_back(0.0);
  try {
    reinterpret_cast<Vjp1>(f)(x.data(), p.data(), w, y.data(), o.wx.data(), o.wp.data(), &n_obs, &n_vars, &n_params,
                              &n_out, &S);
  } catch (...) {
    return false;
  }
  for (const auto* v : {&o.wx, &o.wp})
    for (double a : *v) if (!std::isfinite(a)) return false;
  return true;
}

// The rows `idx` of the derivatives [n, K] of an array of rows `n`, as [idx, K].
std::vector<double> rowsOf(const double* a, int n, const std::vector<int>& idx, int K) {
  const int m = idx.size();
  std::vector<double> out((std::size_t) std::max(m, 1) * std::max(K, 1), 0.0);
  for (int k = 0; k < K; ++k)
    for (int j = 0; j < m; ++j) out[j + (std::size_t) m * k] = a[idx[j] + (std::size_t) n * k];
  return out;
}

std::vector<double> valuesOf(const double* a, const std::vector<int>& idx) {
  std::vector<double> out(idx.size());
  for (std::size_t j = 0; j < idx.size(); ++j) out[j] = a[idx[j]];
  return out;
}

bool anyModelSigma(const In& I) {
  for (int v : I.sna) if (v) return true;
  return false;
}

}  // namespace


// The data rows of all conditions along the tangents, in the order of the conditions and of their
// data: prediction p, sigma s, their derivatives pd, sd [rows, K], the data y, `cens` below the
// limit and `model` for a sigma of the error model, `n` the rows per condition. NULL where a solve
// or a function fails or a value is not finite.
// [[Rcpp::export]]
SEXP chain_rows(List plan, List pre, List raw, int threads) {
  List conds = plan["cond"];
  const int C = conds.size();
  std::vector<In> in;
  if (!gather(conds, pre, raw, true, false, true, in)) return R_NilValue;
  std::vector<int> off(C + 1, 0);
  for (int c = 0; c < C; ++c) off[c + 1] = off[c] + (int) in[c].t.size();
  const int N = off[C], K = C ? in[0].nd : 0;
  NumericVector p(N), s(N), y(N);
  NumericMatrix pd(N, K), sd(N, K);
  LogicalVector cens(N), model(N);
  IntegerVector n(C);
  double *P = REAL(p), *Sg = REAL(s), *Yd = REAL(y), *PD = REAL(pd), *SD = REAL(sd);
  int *Ce = LOGICAL(cens), *Mo = LOGICAL(model);
  std::vector<char> ok(C, 1);
  const int nt = chainThreads(threads, C);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int c = 0; c < C; ++c) {
    const In& I = in[c];
    const bool withErr = anyModelSigma(I);
    Ad Yv, Ev;
    if (!observe(I, false, withErr, Yv, Ev)) {
      ok[c] = 0;
      continue;
    }
    const int T = I.ntime;
    for (std::size_t i = 0; i < I.t.size(); ++i) {
      const int r = off[c] + i, ti = I.t[i], oi = I.o[i];
      P[r] = Yv.y[ti + (std::size_t) T * oi];
      for (int k = 0; k < K; ++k) PD[r + (std::size_t) N * k] = Yv.dy[ti + (std::size_t) T * (oi + (std::size_t) I.nyo * k)];
      Yd[r] = I.ydat[i];
      Ce[r] = I.bloq[i];
      Mo[r] = I.sna[i] != 0;
      if (!I.sna[i]) {
        Sg[r] = I.sfix[i];
      } else {
        const int ei = I.e[i];
        Sg[r] = Ev.y[ti + (std::size_t) T * ei];
        for (int k = 0; k < K; ++k)
          SD[r + (std::size_t) N * k] = Ev.dy[ti + (std::size_t) T * (ei + (std::size_t) I.neo * k)];
      }
      if (!std::isfinite(P[r]) || !std::isfinite(Sg[r])) ok[c] = 0;
    }
  }
  for (int c = 0; c < C; ++c) {
    if (!ok[c]) return R_NilValue;
    n[c] = off[c + 1] - off[c];
  }
  return List::create(Named("p") = p, Named("pd") = pd, Named("s") = s, Named("sd") = sd, Named("y") = y,
                      Named("cens") = cens, Named("model") = model, Named("n") = n);
}


// The seeds `seedP` on the predictions and `seedS` on the sigmas of the rows of chain_rows(),
// [rows, S (K + 1)], through the error model and the observables. Returns list(args, cp, dcp): per
// condition the arguments of the backward solve with its cotangent [T, states, S] and curvature
// [T, states, S, K], and the cotangent of the transformation's outputs so far [outputs, S] with
// its derivatives [outputs, S, K]. With `first` the first order of an adjoint plan: no tangents,
// the raw `_vjp` entries, and the store of the forward solve handed to the backward one. NULL
// where a function fails.
// [[Rcpp::export]]
SEXP chain_back(List plan, List pre, List raw, NumericMatrix seedP, NumericMatrix seedS, int S, int threads,
                bool first) {
  List conds = plan["cond"];
  const int C = conds.size();
  std::vector<In> in;
  if (!gather(conds, pre, raw, !first, false, true, in)) return R_NilValue;
  std::vector<DL_FUNC> yv(C, nullptr), ev(C, nullptr);
  std::vector<int> off(C + 1, 0);
  const char *yName = first ? "y_vjp1" : "y_vjp", *eName = first ? "e_vjp1" : "e_vjp";
  for (int c = 0; c < C; ++c) {
    List pc = conds[c];
    yv[c] = address(pc[yName]);
    if (yv[c] == nullptr) return R_NilValue;
    if (in[c].ef != nullptr && !Rf_isNull(pc[eName])) ev[c] = address(pc[eName]);
    if (anyModelSigma(in[c]) && ev[c] == nullptr) return R_NilValue;
    off[c + 1] = off[c] + (int) in[c].t.size();
  }
  const int N = seedP.nrow(), K = C ? in[0].nd : 0, W = S * (K + 1);
  if (N != off[C] || seedP.ncol() != W || seedS.nrow() != N || seedS.ncol() != W) return R_NilValue;
  const double *SP = REAL(seedP), *SS = REAL(seedS);
  const std::size_t s1 = S, k1 = std::max(K, 1);
  std::vector<std::vector<double>> cot(C), curv(C), cp(C), dcp(C);
  std::vector<char> ok(C, 1);
  const int nt = chainThreads(threads, C);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int c = 0; c < C; ++c) {
    const In& I = in[c];
    const bool withErr = anyModelSigma(I);
    Ad Yv, Ev;
    if (!observe(I, false, withErr, Yv, Ev)) {
      ok[c] = 0;
      continue;
    }
    const std::size_t T = I.ntime, nyo = I.nyo, neo = I.neo, npo = I.npout, nst = I.nstate;
    std::vector<double> wY(T * nyo * s1, 0.0), dwY(T * nyo * s1 * k1, 0.0), wE, dwE;
    if (withErr) {
      wE.assign(T * neo * s1, 0.0);
      dwE.assign(T * neo * s1 * k1, 0.0);
    }
    for (std::size_t i = 0; i < I.t.size(); ++i) {
      const std::size_t r = off[c] + i, ti = I.t[i], oi = I.o[i];
      for (std::size_t s = 0; s < s1; ++s) {
        const std::size_t b = s * (K + 1);
        wY[ti + T * (oi + nyo * s)] += SP[r + N * b];
        for (int k = 0; k < K; ++k) dwY[ti + T * (oi + nyo * (s + s1 * k))] += SP[r + N * (b + 1 + k)];
        if (!I.sna[i]) continue;
        const std::size_t ei = I.e[i];
        wE[ti + T * (ei + neo * s)] += SS[r + N * b];
        for (int k = 0; k < K; ++k) dwE[ti + T * (ei + neo * (s + s1 * k))] += SS[r + N * (b + 1 + k)];
      }
    }
    cp[c].assign(npo * s1, 0.0);
    dcp[c].assign(npo * s1 * k1, 0.0);
    auto toOutputs = [&](const Pull& u, const std::vector<int>& par) {
      const std::size_t m = par.size();
      for (std::size_t s = 0; s < s1; ++s)
        for (std::size_t j = 0; j < m; ++j) {
          cp[c][par[j] + npo * s] += u.wp[j + m * s];
          for (int k = 0; k < K; ++k) dcp[c][par[j] + npo * (s + s1 * k)] += u.dwp[j + m * (s + s1 * k)];
        }
    };
    if (withErr) {
      std::vector<double> x = columns(Yv.y.data(), T, nyo, I.evar, 1), vx = columns(Yv.dy.data(), T, nyo, I.evar, K);
      std::vector<double> p = valuesOf(I.pv, I.epar), vp = rowsOf(I.pj, npo, I.epar, K);
      Pull u;
      if (!(first ? vjp1(ev[c], x, p, wE.data(), T, I.evar.size(), I.epar.size(), neo, S, u)
                  : vjpAd(ev[c], x, p, vx.data(), vp.data(), wE.data(), dwE.data(), T, I.evar.size(), I.epar.size(),
                          neo, S, K, u))) {
        ok[c] = 0;
        continue;
      }
      const std::size_t m = I.evar.size();
      for (std::size_t s = 0; s < s1; ++s)
        for (std::size_t j = 0; j < m; ++j)
          for (std::size_t t = 0; t < T; ++t) {
            wY[t + T * (I.evar[j] + nyo * s)] += u.wx[t + T * (j + m * s)];
            for (int k = 0; k < K; ++k)
              dwY[t + T * (I.evar[j] + nyo * (s + s1 * k))] += u.dwx[t + T * (j + m * (s + s1 * k))];
          }
      toOutputs(u, I.epar);
    }
    std::vector<double> x = columns(I.xv, T, nst, I.yvar, 1), vx = columns(I.xs, T, nst, I.yvar, K);
    std::vector<double> p = valuesOf(I.pv, I.ypar), vp = rowsOf(I.pj, npo, I.ypar, K);
    Pull u;
    if (!(first ? vjp1(yv[c], x, p, wY.data(), T, I.yvar.size(), I.ypar.size(), nyo, S, u)
                : vjpAd(yv[c], x, p, vx.data(), vp.data(), wY.data(), dwY.data(), T, I.yvar.size(), I.ypar.size(),
                        nyo, S, K, u))) {
      ok[c] = 0;
      continue;
    }
    toOutputs(u, I.ypar);
    const std::size_t m = I.yvar.size();
    cot[c].assign(T * nst * s1, 0.0);
    curv[c].assign(T * nst * s1 * k1, 0.0);
    for (std::size_t s = 0; s < s1; ++s)
      for (std::size_t j = 0; j < m; ++j)
        for (std::size_t t = 0; t < T; ++t) {
          cot[c][t + T * (I.yvar[j] + nst * s)] += u.wx[t + T * (j + m * s)];
          for (int k = 0; k < K; ++k)
            curv[c][t + T * (I.yvar[j] + nst * (s + s1 * k))] += u.dwx[t + T * (j + m * (s + s1 * k))];
        }
  }
  for (int c = 0; c < C; ++c) if (!ok[c]) return R_NilValue;

  List args(C), cpL(C), dcpL(C);
  for (int c = 0; c < C; ++c) {
    List pc = conds[c];
    const In& I = in[c];
    SEXP prepared = pc["xr_args"];
    List a = Rf_shallow_duplicate(prepared);
    const int atParms = position(a, "parms"), atTangent = position(a, "tangent"), atCot = position(a, "cotangent"),
              atCurv = position(a, "curvature");
    if (atParms < 0 || atCot < 0 || (!first && (atTangent < 0 || atCurv < 0))) return R_NilValue;
    IntegerVector xp = pc["x_parm"];
    const int nreq = xp.size();
    const std::size_t npo = I.npout;
    NumericVector parms(nreq);
    for (int q = 0; q < nreq; ++q) parms[q] = I.pv[xp[q]];
    a[atParms] = parms;
    RObject oldC = a[atCot];
    NumericVector w(cot[c].begin(), cot[c].end());
    w.attr("dim") = IntegerVector::create(I.ntime, I.nstate, S);
    w.attr("dimnames") = oldC.attr("dimnames");
    a[atCot] = w;
    if (first) {
      // the checkpoints of the forward solve, so that the backward one integrates nothing
      List r = raw[c];
      const int atStore = position(r, "store");
      if (atStore < 0) return R_NilValue;
      a.push_back(r[atStore], "store");
      args[c] = a;
      cpL[c] = NumericVector(cp[c].begin(), cp[c].end());
      dcpL[c] = NumericVector(0);
      continue;
    }
    NumericMatrix tg(nreq, K);
    for (int q = 0; q < nreq; ++q)
      for (int k = 0; k < K; ++k) tg[q + (std::size_t) nreq * k] = I.pj[xp[q] + npo * k];
    RObject oldT = a[atTangent];
    tg.attr("dimnames") = oldT.attr("dimnames");
    a[atTangent] = tg;
    NumericVector dw(curv[c].begin(), curv[c].end());
    dw.attr("dim") = IntegerVector::create(I.ntime, I.nstate, S, K);
    a[atCurv] = dw;
    args[c] = a;
    cpL[c] = NumericVector(cp[c].begin(), cp[c].end());
    dcpL[c] = NumericVector(dcp[c].begin(), dcp[c].end());
  }
  return List::create(Named("args") = args, Named("cp") = cpL, Named("dcp") = dcpL);
}


namespace {

// Cotangents w [n, S] of a stage's outputs and their derivatives dw [n, S, K] along the tangents,
// pulled through the stage onto its input `u` (the vector before it with its tangents u.J [n, K]).
// First order (`first`) through the raw `_vjp` entry, else through `_vjp_ad`; a steady state by
// D' w and D' dw + (D2 u.J)' w. False where a function fails or a root is missing.
bool pullStage(const PStage& S, const Root& Q, const Vec& u, int nout, int Sn, int K, bool first,
               std::vector<double>& w, std::vector<double>& dw) {
  const int nin = S.in.size(), nfo = S.nfo, n = u.n;
  const std::size_t s1 = Sn;
  std::vector<double> wp((std::size_t) n * s1, 0.0), dwp((std::size_t) n * s1 * K, 0.0);
  for (std::size_t q = 0; q < S.pass.size(); ++q)
    for (std::size_t s = 0; s < s1; ++s) {
      const std::size_t i = nfo + q, k = S.pass[q];
      wp[k + n * s] += w[i + nout * s];
      for (int t = 0; t < K; ++t) dwp[k + n * (s + s1 * t)] += dw[i + nout * (s + s1 * t)];
    }
  if (S.kind == 0) {
    if (S.vjp == nullptr) return false;
    std::vector<double> pin(nin), vp((std::size_t) std::max(nin, 1) * std::max(K, 1), 0.0),
        wf((std::size_t) nfo * s1), dwf((std::size_t) std::max<std::size_t>(nfo * s1 * K, 1), 0.0);
    for (int j = 0; j < nin; ++j) {
      pin[j] = u.v[S.in[j]];
      for (int t = 0; t < K; ++t) vp[j + (std::size_t) nin * t] = u.J[S.in[j] + (std::size_t) n * t];
    }
    for (std::size_t s = 0; s < s1; ++s)
      for (int i = 0; i < nfo; ++i) {
        wf[i + nfo * s] = w[i + nout * s];
        for (int t = 0; t < K; ++t) dwf[i + nfo * (s + s1 * t)] = dw[i + nout * (s + s1 * t)];
      }
    std::vector<double> none;
    Pull r;
    const bool ok = first ? vjp1(S.vjp, none, pin, wf.data(), 1, 0, nin, nfo, Sn, r)
                          : vjpAd(S.vjp, none, pin, nullptr, vp.data(), wf.data(), dwf.data(), 1, 0, nin, nfo, Sn, K, r);
    if (!ok) return false;
    for (int j = 0; j < nin; ++j)
      for (std::size_t s = 0; s < s1; ++s) {
        wp[S.in[j] + n * s] += r.wp[j + nin * s];
        for (int t = 0; t < K; ++t) dwp[S.in[j] + n * (s + s1 * t)] += r.dwp[j + nin * (s + s1 * t)];
      }
  } else {
    const int nq = Q.ndep;
    if (Q.D == nullptr || (K > 0 && Q.D2 == nullptr)) return false;
    for (int j = 0; j < nin; ++j) {
      if (S.in[j] < 0) continue;
      for (std::size_t s = 0; s < s1; ++s) {
        double a = 0.0;
        for (int i = 0; i < nq; ++i) a += Q.D[i + (std::size_t) nq * j] * w[i + nout * s];
        wp[S.in[j] + n * s] += a;
        for (int t = 0; t < K; ++t) {
          double b = 0.0;
          for (int i = 0; i < nq; ++i) b += Q.D[i + (std::size_t) nq * j] * dw[i + nout * (s + s1 * t)];
          for (int l = 0; l < nin; ++l) {
            if (S.in[l] < 0) continue;
            const double Jl = u.J[S.in[l] + (std::size_t) n * t];
            if (Jl == 0.0) continue;
            for (int i = 0; i < nq; ++i)
              b += Q.D2[i + (std::size_t) nq * (j + (std::size_t) nin * l)] * w[i + nout * s] * Jl;
          }
          dwp[S.in[j] + n * (s + s1 * t)] += b;
        }
      }
    }
  }
  w.swap(wp);
  dw.swap(dwp);
  return true;
}

}  // namespace

// The answers `rraw` of the backward solve on the arguments of chain_back() `back`, added to its
// cotangents of the transformation's outputs and pulled through the stages of the transformation
// at `values` along `tangents`, with the `roots` of chain_params(). Per condition list(rows, ct):
// the free names it reaches, 0-based rows of `tangents`, and their cotangents [rows, S (K + 1)] in
// blocks. With `first` the first order of an adjoint plan, through the raw `_vjp` entries. NULL
// where a solve or a function fails.
// [[Rcpp::export]]
SEXP chain_pull(List plan, List values, NumericMatrix tangents, SEXP roots, List back, List rraw, int S, int threads,
                bool first) {
  List conds = plan["cond"];
  const int C = conds.size();
  const bool shared = values.size() == 1;
  const int K = tangents.ncol();
  const std::size_t s1 = S;
  List cpL = back["cp"], dcpL = back["dcp"];
  std::vector<PPlan> pp(C);
  std::vector<std::vector<Root>> rt(C);
  std::vector<Vec> u0(C);
  std::vector<std::vector<double>> cp(C), dcp(C);
  for (int c = 0; c < C; ++c) {
    List pc = conds[c];
    if (!readPlan(pc, pp[c])) return R_NilValue;
    if (readRoots(Rf_isNull(roots) ? R_NilValue : VECTOR_ELT(roots, c), pp[c], rt[c]) >= 0) return R_NilValue;
    NumericVector v = values[shared ? 0 : c];
    u0[c] = outerVec(pp[c], v, K, false, &tangents);
    IntegerVector xp = pc["x_parm"];
    const int nreq = xp.size();
    cp[c] = as<std::vector<double>>(cpL[c]);
    dcp[c] = as<std::vector<double>>(dcpL[c]);
    const std::size_t npo = cp[c].size() / s1;
    // the cotangent of the solve's inputs [inputs, S] and its derivatives [inputs, K, S]
    List r = rraw[c];
    const int dg = position(r, "diagnostics");
    if (dg >= 0) {
      List d = r[dg];
      const int rc = position(d, "return_code");
      if (rc >= 0 && as<int>(d[rc]) != 0) return R_NilValue;
    }
    const int atCot = position(r, "cotangent"), atCurv = position(r, "curvature");
    if (atCot < 0 || (!first && atCurv < 0)) return R_NilValue;
    NumericVector w = r[atCot], dw = first ? NumericVector(0) : NumericVector(r[atCurv]);
    if ((std::size_t) w.size() != (std::size_t) nreq * s1 || (std::size_t) dw.size() != (std::size_t) nreq * K * s1)
      return R_NilValue;
    for (std::size_t s = 0; s < s1; ++s)
      for (int q = 0; q < nreq; ++q) {
        cp[c][xp[q] + npo * s] += w[q + (std::size_t) nreq * s];
        for (int k = 0; k < K; ++k) dcp[c][xp[q] + npo * (s + s1 * k)] += dw[q + (std::size_t) nreq * (k + K * s)];
      }
  }
  std::vector<char> ok(C, 1);
  const int nt = chainThreads(threads, C);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int c = 0; c < C; ++c) {
    const PPlan& P = pp[c];
    std::vector<Vec> ins;
    Vec u = u0[c];
    if (!runStages(P, rt[c], u, P.st.size(), K, false, &ins)) {
      ok[c] = 0;
      continue;
    }
    int nout = u.n;
    for (int s = (int) P.st.size() - 1; s >= 0 && ok[c]; --s) {
      ok[c] = pullStage(P.st[s], rt[c][s], ins[s], nout, S, K, first, cp[c], dcp[c]);
      nout = ins[s].n;
    }
  }
  List out(C);
  for (int c = 0; c < C; ++c) {
    if (!ok[c]) return R_NilValue;
    const std::vector<int>& vrow = pp[c].ovrow;
    const int nin = vrow.size(), Wd = S * (K + 1);
    const std::vector<double> &wp = cp[c], &dwp = dcp[c];
    for (double a : wp) if (!std::isfinite(a)) return R_NilValue;
    for (double a : dwp) if (!std::isfinite(a)) return R_NilValue;
    std::vector<int> keep;
    for (int j = 0; j < nin; ++j) {
      if (vrow[j] < 0) continue;
      bool any = false;
      for (std::size_t s = 0; s < s1 && !any; ++s) {
        any = wp[j + nin * s] != 0.0;
        for (int k = 0; k < K && !any; ++k) any = dwp[j + nin * (s + s1 * k)] != 0.0;
      }
      if (any) keep.push_back(j);
    }
    const int m = keep.size();
    NumericMatrix ct(m, Wd);
    IntegerVector rows(m);
    for (int q = 0; q < m; ++q) {
      const int j = keep[q];
      rows[q] = vrow[j];
      for (std::size_t s = 0; s < s1; ++s) {
        ct(q, s * (K + 1)) = wp[j + nin * s];
        for (int k = 0; k < K; ++k) ct(q, s * (K + 1) + 1 + k) = dwp[j + nin * (s + s1 * k)];
      }
    }
    out[c] = List::create(Named("rows") = rows, Named("ct") = ct);
  }
  return out;
}


// The local data terms `local` of chain_data(), one per request, summed per unit: request r adds to
// unit owner[r] (0-based), its gradient and Hessian by the names of the unit's `etaNames`. Returns
// per unit list(value, gradient, hessian) in the order of its names.
// [[Rcpp::export]]
List chain_unit_sums(List local, IntegerVector owner, List etaNames) {
  const int U = etaNames.size(), R = local.size();
  std::vector<double> value(U, 0.0);
  std::vector<std::vector<double>> grad(U), hess(U);
  std::vector<std::unordered_map<std::string, int>> pos(U);
  for (int u = 0; u < U; ++u) {
    CharacterVector nm = etaNames[u];
    const int K = nm.size();
    for (int k = 0; k < K; ++k) pos[u][std::string(nm[k])] = k;
    grad[u].assign(K, 0.0);
    hess[u].assign((std::size_t) K * K, 0.0);
  }
  for (int r = 0; r < R; ++r) {
    const int u = owner[r];
    List l = local[r];
    value[u] += as<double>(l["value"]);
    NumericVector g = l["gradient"];
    const int n = g.size();
    if (!n) continue;
    CharacterVector gn = g.names();
    std::vector<int> at(n);
    for (int j = 0; j < n; ++j) {
      auto hit = pos[u].find(std::string(gn[j]));
      at[j] = hit == pos[u].end() ? -1 : hit->second;
    }
    NumericMatrix h = l["hessian"];
    const std::size_t K = grad[u].size();
    for (int j = 0; j < n; ++j) {
      if (at[j] < 0) continue;
      grad[u][at[j]] += g[j];
      for (int i = 0; i < n; ++i)
        if (at[i] >= 0) hess[u][at[i] + K * at[j]] += h(i, j);
    }
  }
  List out(U);
  for (int u = 0; u < U; ++u) {
    CharacterVector nm = etaNames[u];
    const int K = nm.size();
    NumericVector g(grad[u].begin(), grad[u].end());
    g.names() = nm;
    NumericMatrix H(K, K);
    std::copy(hess[u].begin(), hess[u].end(), H.begin());
    H.attr("dimnames") = List::create(nm, nm);
    out[u] = List::create(Named("value") = value[u], Named("gradient") = g, Named("hessian") = H);
  }
  return out;
}
