// Alignment models: the derivatives of the readings of every condition in its parameters from the
// derivatives of the compiled model in the symbols, and the vjp contractions of the backward sweep.
//
// A reading r of a condition depends on the parameters of the condition through the predictors of
// the symbols, eta_a(r) = sum_p D_a[r, p] theta_p, with D_a the dense coefficients of symbol a on
// the columns of its parameters. With t[r, o, a] = d f_o / d eta_a and h[r, o, a, b] the second
// derivatives of output o (prediction, sigma),
//
//   d f_o(r) / d theta_p            = sum_a t[r, o, a] D_a[r, p],
//   d2 f_o(r) / d theta_p d theta_q = sum_{a, b} h[r, o, a, b] D_a[r, p] D_b[r, q],
//
// and the vjp of a seed w[r, k] is J' w, at second order sum_r w[r] Hess(f_o(r)) v plus J' wdot
// for a tangent v of the parameters. The model arrays cover all rows of all conditions; a condition
// reads its rows from `start` on. R objects are allocated in the serial part, the OpenMP region
// writes into their memory only.

#include <Rcpp.h>
#include <vector>
#include <algorithm>

#if !defined(_WIN32)
#include <pthread.h>
namespace {
bool alignmentForked = false;
struct AlignmentForkGuard {
  AlignmentForkGuard() { pthread_atfork(nullptr, nullptr, [] { alignmentForked = true; }); }
} alignmentForkGuard;
}
static inline bool alignment_forked() { return alignmentForked; }
#else
static inline bool alignment_forked() { return false; }
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace Rcpp;

namespace {

// One symbol of one condition: the local parameters of its columns and its dense coefficients.
struct SymbolBlock {
  std::vector<int> cols;      // local parameter index of every column, 0-based
  const double* D = nullptr;  // [n x cols.size()], column-major
};

struct CondBlock {
  int start = 0, n = 0, P = 0;
  std::vector<SymbolBlock> sym;
  std::vector<int> slot;      // per local parameter its output column, -1 when not kept
  int F = 0;                  // number of kept parameters
};

int nThreads(int threads) {
  int n = std::max(1, threads);
#ifndef _OPENMP
  n = 1;
#endif
  if (alignment_forked()) n = 1;
  return n;
}

// The blocks of the conditions from the R description: per condition `start`, `n`, `P`, `cols`
// (per symbol 1-based local indices), `dense` (per symbol a numeric matrix) and `keep` (1-based
// local indices of the parameters in the output, in output order; none without it).
std::vector<CondBlock> readBlocks(const List& conds, int S) {
  std::vector<CondBlock> out(conds.size());
  for (int c = 0; c < conds.size(); ++c) {
    List cd = conds[c];
    CondBlock& B = out[c];
    B.start = as<int>(cd["start"]);
    B.n = as<int>(cd["n"]);
    B.P = as<int>(cd["P"]);
    List cols = cd["cols"], dense = cd["dense"];
    if (cols.size() != S || dense.size() != S) stop("alignment kernel: one block per symbol.");
    B.sym.resize(S);
    for (int a = 0; a < S; ++a) {
      IntegerVector ca = cols[a];
      NumericMatrix Da = dense[a];
      if (Da.nrow() != B.n || Da.ncol() != ca.size()) stop("alignment kernel: coefficients of the wrong shape.");
      B.sym[a].cols.resize(ca.size());
      for (int j = 0; j < ca.size(); ++j) B.sym[a].cols[j] = ca[j] - 1;
      B.sym[a].D = Da.size() ? &Da[0] : nullptr;
    }
    B.slot.assign(B.P, -1);
    B.F = 0;
    if (cd.containsElementNamed("keep")) {
      IntegerVector keep = cd["keep"];
      for (int j = 0; j < keep.size(); ++j) B.slot[keep[j] - 1] = j;
      B.F = keep.size();
    }
  }
  return out;
}

} // namespace

// The predictors of all rows, [N x S]: for the rows of a condition sum_j D_a[r, j] theta[cols_a[j]]
// with `theta` its parameter values, one vector per condition over its local parameters.
// [[Rcpp::export]]
NumericMatrix alignmentPredictors(List conds, List theta, int S, int threads = 1) {
  std::vector<CondBlock> B = readBlocks(conds, S);
  const int C = B.size();
  int N = 0;
  for (int c = 0; c < C; ++c) N = std::max(N, B[c].start + B[c].n);
  NumericMatrix out(N, S);
  double* o = N * S ? &out[0] : nullptr;
  std::vector<const double*> th(C);
  for (int c = 0; c < C; ++c) {
    NumericVector tc = theta[c];
    if (tc.size() != B[c].P) stop("alignment kernel: parameter values of the wrong length.");
    th[c] = tc.size() ? &tc[0] : nullptr;
  }
  const int nt = nThreads(threads);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
#endif
  for (int c = 0; c < C; ++c) {
    const CondBlock& b = B[c];
    for (int a = 0; a < S; ++a) {
      const SymbolBlock& sa = b.sym[a];
      for (std::size_t j = 0; j < sa.cols.size(); ++j) {
        const double v = th[c][sa.cols[j]];
        if (v == 0.0) continue;
        for (int r = 0; r < b.n; ++r) o[b.start + r + (std::size_t) N * a] += sa.D[r + (std::size_t) b.n * j] * v;
      }
    }
  }
  return out;
}

// Per condition the derivative arrays of prediction and sigma over the kept parameters,
// [n, 1, F] and with `hessian` [n, 1, F, F]; `tangent` is [N, 2, S], `hessian` [N, 2, S, S].
// [[Rcpp::export]]
List alignmentForwardBlocks(NumericVector tangent, Nullable<NumericVector> hessian, List conds,
                            List keepNames, int threads = 1) {
  IntegerVector dt = tangent.attr("dim");
  const int N = dt[0], S = dt[2];
  std::vector<CondBlock> B = readBlocks(conds, S);
  const int C = B.size();
  const bool second = hessian.isNotNull();
  NumericVector H2;
  if (second) H2 = NumericVector(hessian);
  const double* t = &tangent[0];
  const double* h = second ? &H2[0] : nullptr;

  // allocation and dimnames in the serial part
  List out(C);
  std::vector<double*> Jp(C), Js(C), Hp(C, nullptr), Hs(C, nullptr);
  for (int c = 0; c < C; ++c) {
    const int n = B[c].n, F = B[c].F;
    List dn = List::create(R_NilValue, CharacterVector::create("y"), keepNames[c]);
    NumericVector jp(n * F), js(n * F);
    jp.attr("dim") = IntegerVector::create(n, 1, F); jp.attr("dimnames") = dn;
    js.attr("dim") = IntegerVector::create(n, 1, F); js.attr("dimnames") = dn;
    Jp[c] = n * F ? &jp[0] : nullptr; Js[c] = n * F ? &js[0] : nullptr;
    List oc = List::create(Named("prediction") = jp, Named("sigma") = js);
    if (second) {
      List dn2 = List::create(R_NilValue, CharacterVector::create("y"), keepNames[c], keepNames[c]);
      NumericVector hp(n * F * F), hs(n * F * F);
      hp.attr("dim") = IntegerVector::create(n, 1, F, F); hp.attr("dimnames") = dn2;
      hs.attr("dim") = IntegerVector::create(n, 1, F, F); hs.attr("dimnames") = dn2;
      Hp[c] = n * F ? &hp[0] : nullptr; Hs[c] = n * F ? &hs[0] : nullptr;
      oc["prediction2"] = hp; oc["sigma2"] = hs;
    }
    out[c] = oc;
  }

  const int nt = nThreads(threads);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
#endif
  for (int c = 0; c < C; ++c) {
    const CondBlock& b = B[c];
    const int n = b.n, F = b.F;
    double* J[2] = {Jp[c], Js[c]};
    double* H[2] = {Hp[c], Hs[c]};
    for (int r = 0; r < n; ++r) {
      const int R = b.start + r;
      for (int a = 0; a < S; ++a) {
        const SymbolBlock& sa = b.sym[a];
        for (std::size_t j = 0; j < sa.cols.size(); ++j) {
          const int p = b.slot[sa.cols[j]];
          if (p < 0) continue;
          const double d = sa.D[r + (std::size_t) n * j];
          if (d == 0.0) continue;
          for (int o = 0; o < 2; ++o)
            J[o][r + (std::size_t) n * p] += t[R + (std::size_t) N * (o + 2 * a)] * d;
          if (!second) continue;
          for (int bb = 0; bb < S; ++bb) {
            const SymbolBlock& sb = b.sym[bb];
            for (std::size_t k = 0; k < sb.cols.size(); ++k) {
              const int q = b.slot[sb.cols[k]];
              if (q < 0) continue;
              const double dd = d * sb.D[r + (std::size_t) n * k];
              if (dd == 0.0) continue;
              for (int o = 0; o < 2; ++o)
                H[o][r + (std::size_t) n * (p + (std::size_t) F * q)] +=
                  h[R + (std::size_t) N * (o + 2 * (a + (std::size_t) S * bb))] * dd;
            }
          }
        }
      }
    }
  }
  return out;
}

// Per condition the vjp of output `o` (1 prediction, 2 sigma) over its local parameters,
// [P x K]: per seed block J' w, and with tangents `V` (per condition [P x D] or NULL) the
// derivatives sum_r w[r, first] Hess(f(r)) v_d added to the columns `rest`. `W` per condition is
// [n x K]; `first` (length nb) and `rest` ([D x nb]) are 1-based.
// [[Rcpp::export]]
List alignmentVjpBlocks(NumericVector tangent, Nullable<NumericVector> hessian, int o, List conds,
                        List W, List V, IntegerVector first, IntegerMatrix rest, int threads = 1) {
  IntegerVector dt = tangent.attr("dim");
  const int N = dt[0], S = dt[2];
  std::vector<CondBlock> B = readBlocks(conds, S);
  const int C = B.size();
  const int nb = first.size(), D = rest.nrow();
  const bool second = hessian.isNotNull() && D > 0;
  NumericVector H2;
  if (second) H2 = NumericVector(hessian);
  const double* t = &tangent[0];
  const double* h = second ? &H2[0] : nullptr;
  const int oi = o - 1;

  List out(C);
  std::vector<double*> ct(C);
  std::vector<const double*> w(C), v(C, nullptr);
  std::vector<int> K(C);
  for (int c = 0; c < C; ++c) {
    NumericMatrix Wc = W[c];
    if (Wc.nrow() != B[c].n) stop("alignment kernel: a seed of the wrong length.");
    K[c] = Wc.ncol();
    w[c] = Wc.size() ? &Wc[0] : nullptr;
    if (second) {
      NumericMatrix Vc = V[c];
      if (Vc.nrow() != B[c].P || Vc.ncol() < D) stop("alignment kernel: tangents of the wrong shape.");
      v[c] = Vc.size() ? &Vc[0] : nullptr;
    }
    NumericMatrix oc(B[c].P, K[c]);
    ct[c] = oc.size() ? &oc[0] : nullptr;
    out[c] = oc;
  }
  std::vector<int> fi(nb), re((std::size_t) D * nb);
  for (int s = 0; s < nb; ++s) fi[s] = first[s] - 1;
  for (int s = 0; s < nb; ++s) for (int d = 0; d < D; ++d) re[d + (std::size_t) D * s] = rest(d, s) - 1;

  const int nt = nThreads(threads);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
#endif
  for (int c = 0; c < C; ++c) {
    const CondBlock& b = B[c];
    const int n = b.n, P = b.P, Kc = K[c];
    double* out_c = ct[c];
    const double* wc = w[c];
    // Z_b[r, d] = sum_q D_b[r, q] V[col_q, d], the tangent of every predictor along the D directions
    std::vector<double> Z;
    if (second) {
      Z.assign((std::size_t) n * D * S, 0.0);
      for (int bb = 0; bb < S; ++bb) {
        const SymbolBlock& sb = b.sym[bb];
        for (std::size_t k = 0; k < sb.cols.size(); ++k) {
          const int q = sb.cols[k];
          for (int r = 0; r < n; ++r) {
            const double dq = sb.D[r + (std::size_t) n * k];
            if (dq == 0.0) continue;
            for (int d = 0; d < D; ++d)
              Z[r + (std::size_t) n * (d + (std::size_t) D * bb)] += dq * v[c][q + (std::size_t) P * d];
          }
        }
      }
    }
    for (int a = 0; a < S; ++a) {
      const SymbolBlock& sa = b.sym[a];
      for (std::size_t j = 0; j < sa.cols.size(); ++j) {
        const int p = sa.cols[j];
        for (int r = 0; r < n; ++r) {
          const double d = sa.D[r + (std::size_t) n * j];
          if (d == 0.0) continue;
          const int R = b.start + r;
          const double ta = t[R + (std::size_t) N * (oi + 2 * a)] * d;
          for (int k = 0; k < Kc; ++k) out_c[p + (std::size_t) P * k] += ta * wc[r + (std::size_t) n * k];
          if (!second) continue;
          for (int bb = 0; bb < S; ++bb) {
            const double hab = h[R + (std::size_t) N * (oi + 2 * (a + (std::size_t) S * bb))] * d;
            if (hab == 0.0) continue;
            for (int s = 0; s < nb; ++s) {
              const double ws = wc[r + (std::size_t) n * fi[s]] * hab;
              if (ws == 0.0) continue;
              for (int dd = 0; dd < D; ++dd)
                out_c[p + (std::size_t) P * re[dd + (std::size_t) D * s]] +=
                  ws * Z[r + (std::size_t) n * (dd + (std::size_t) D * bb)];
            }
          }
        }
      }
    }
  }
  return out;
}
