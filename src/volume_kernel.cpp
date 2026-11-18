// The volume terms of the FOCEI marginal of one condition (.volumeTerms()), summed over the rows. Per row
// r with loss derivatives L_ij = d^(i+j) l / dp^i ds^j, the first derivatives a, b of p and s in the
// etas, Ha = H^-1 a, Hb = H^-1 b and the quadratic forms qaa = a'Ha, qab = a'Hb, qbb = b'Hb:
//   F = (Wp; Ws), the derivatives of (a; b) in theta,   C = (Vp; Vs), those in the etas,
//   T = (pA; sA), the derivatives of (p; s) in theta,   c = (a; b),
//   G = [Lpp Lps; Lps Lss] (x) H^-1,   M = Gv L3 with Gv = [Ha Hb 0; 0 Ha Hb], L3 = [L30 L21; L21 L12; L12 L03],
//   Z = [cpp cps; cps css], h = (c1, c2), the contractions of the third and fourth loss derivatives with q.
// The gradient is F'g + T'h over theta and C'g + c'h over the etas, g = 2 G (a; b); the Hessian parts
// are the sums of X' w Y over the rows, and dH by theta and by the etas is vec(x V' + V x') +
// vec(y W' + W y') over F = (x, y) plus the loss terms over T.

#include <Rcpp.h>

#include <cstddef>
#include <vector>

using namespace Rcpp;

namespace {

// out (ni x nj) += X' w Y for X (dk x ni), w (dk x dl), Y (dl x nj), all column-major
void addXtwY(double* out, const double* X, const double* w, const double* Y, int dk, int dl, int ni, int nj,
             double scale, std::vector<double>& wY) {
  wY.assign((std::size_t) dk * nj, 0.0);
  for (int j = 0; j < nj; ++j)
    for (int l = 0; l < dl; ++l) {
      const double y = Y[l + (std::size_t) dl * j];
      if (y == 0.0) continue;
      for (int k = 0; k < dk; ++k) wY[k + (std::size_t) dk * j] += w[k + (std::size_t) dk * l] * y;
    }
  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i) {
      double s = 0.0;
      for (int k = 0; k < dk; ++k) s += X[k + (std::size_t) dk * i] * wY[k + (std::size_t) dk * j];
      out[i + (std::size_t) ni * j] += scale * s;
    }
}

// rows (d x m) of row r of an array [n, d, m] stacked under rows of a second one: out (2d x m)
void stackRow(double* out, const double* A, const double* B, int r, int n, int d, int m) {
  for (int j = 0; j < m; ++j)
    for (int k = 0; k < d; ++k) {
      out[k + (std::size_t) 2 * d * j] = A[r + (std::size_t) n * (k + (std::size_t) d * j)];
      out[d + k + (std::size_t) 2 * d * j] = B[r + (std::size_t) n * (k + (std::size_t) d * j)];
    }
}

} // namespace

// [[Rcpp::export]]
List volume_terms(NumericVector loss, NumericVector a, NumericVector b, NumericVector Wp, NumericVector Ws,
                  NumericVector Vp, NumericVector Vs, NumericVector pA, NumericVector sA, SEXP pAA, SEXP sAA,
                  NumericMatrix Hinv, int n, int q, bool hessian) {
  const int K = Hinv.nrow(), K2 = 2 * K, KK = K * K;
  if (loss.size() != (R_xlen_t) n * 25 || a.size() != (R_xlen_t) n * K || b.size() != (R_xlen_t) n * K ||
      Wp.size() != (R_xlen_t) n * K * q || Ws.size() != (R_xlen_t) n * K * q ||
      Vp.size() != (R_xlen_t) n * KK || Vs.size() != (R_xlen_t) n * KK ||
      pA.size() != (R_xlen_t) n * q || sA.size() != (R_xlen_t) n * q)
    stop("volume_terms: the row arrays do not match n, K and q.");
  const double* pAAp = nullptr;
  const double* sAAp = nullptr;
  if (hessian) {
    NumericVector pAAv(pAA), sAAv(sAA);
    if (pAAv.size() != (R_xlen_t) n * q * q || sAAv.size() != (R_xlen_t) n * q * q)
      stop("volume_terms: the second derivatives in theta do not match n and q.");
    pAAp = REAL(pAA); sAAp = REAL(sAA);
  }
  const double* H = Hinv.begin();
  auto at = [n](const double* x, int r, int k) { return x[r + (std::size_t) n * k]; };

  NumericVector g(q), gE(K);
  NumericMatrix local(hessian ? q : 0, hessian ? q : 0), N(hessian ? q : 0, hessian ? K : 0),
      Q(hessian ? K : 0, hessian ? K : 0), D(hessian ? KK : 0, hessian ? q : 0),
      Gam(hessian ? KK : 0, hessian ? K : 0);
  std::vector<double> X2((std::size_t) q * q, 0.0), Q2((std::size_t) KK, 0.0);
  std::vector<double> ar(K), br(K), Ha(K), Hb(K), gw(K2), V(K), Wv(K);
  std::vector<double> S((std::size_t) K2 * q), C((std::size_t) K2 * K), T((std::size_t) 2 * q),
      cm((std::size_t) 2 * K), G((std::size_t) K2 * K2), M((std::size_t) K2 * 2), tM((std::size_t) 2 * K2),
      Z(4), work;

  for (int r = 0; r < n; ++r) {
    auto L = [&](int i, int j) { return loss[r + (std::size_t) n * (i + 5 * j)]; };
    const double Lpp = L(2, 0), Lps = L(1, 1), Lss = L(0, 2);
    for (int k = 0; k < K; ++k) { ar[k] = at(a.begin(), r, k); br[k] = at(b.begin(), r, k); }
    double qaa = 0.0, qab = 0.0, qbb = 0.0;
    for (int k = 0; k < K; ++k) {
      double sa = 0.0, sb = 0.0;
      for (int l = 0; l < K; ++l) { sa += ar[l] * H[l + (std::size_t) K * k]; sb += br[l] * H[l + (std::size_t) K * k]; }
      Ha[k] = sa; Hb[k] = sb;
      qaa += sa * ar[k]; qab += sa * br[k]; qbb += sb * br[k];
    }
    const double c1 = qaa * L(3, 0) + 2 * qab * L(2, 1) + qbb * L(1, 2);
    const double c2 = qaa * L(2, 1) + 2 * qab * L(1, 2) + qbb * L(0, 3);
    for (int k = 0; k < K; ++k) {
      gw[k] = 2 * (Lpp * Ha[k] + Lps * Hb[k]);
      gw[K + k] = 2 * (Lps * Ha[k] + Lss * Hb[k]);
    }
    stackRow(S.data(), Wp.begin(), Ws.begin(), r, n, K, q);
    stackRow(C.data(), Vp.begin(), Vs.begin(), r, n, K, K);
    for (int j = 0; j < q; ++j) { T[2 * j] = at(pA.begin(), r, j); T[2 * j + 1] = at(sA.begin(), r, j); }
    for (int m = 0; m < K; ++m) { cm[2 * m] = ar[m]; cm[2 * m + 1] = br[m]; }

    for (int j = 0; j < q; ++j) {
      double s = T[2 * j] * c1 + T[2 * j + 1] * c2;
      for (int k = 0; k < K2; ++k) s += S[k + (std::size_t) K2 * j] * gw[k];
      g[j] += s;
    }
    for (int m = 0; m < K; ++m) {
      double s = ar[m] * c1 + br[m] * c2;
      for (int k = 0; k < K2; ++k) s += C[k + (std::size_t) K2 * m] * gw[k];
      gE[m] += s;
    }
    if (!hessian) continue;

    const double Lm[4] = {Lpp, Lps, Lps, Lss};
    for (int l = 0; l < K2; ++l)
      for (int k = 0; k < K2; ++k)
        G[k + (std::size_t) K2 * l] = Lm[(k / K) + 2 * (l / K)] * H[(k % K) + (std::size_t) K * (l % K)];
    const double L3[3][2] = {{L(3, 0), L(2, 1)}, {L(2, 1), L(1, 2)}, {L(1, 2), L(0, 3)}};
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < K; ++k) {
        M[k + (std::size_t) K2 * j] = Ha[k] * L3[0][j] + Hb[k] * L3[1][j];
        M[K + k + (std::size_t) K2 * j] = Ha[k] * L3[1][j] + Hb[k] * L3[2][j];
      }
    for (int k = 0; k < K2; ++k) for (int j = 0; j < 2; ++j) tM[j + 2 * (std::size_t) k] = M[k + (std::size_t) K2 * j];
    const double cpp = qaa * L(4, 0) + 2 * qab * L(3, 1) + qbb * L(2, 2);
    const double cps = qaa * L(3, 1) + 2 * qab * L(2, 2) + qbb * L(1, 3);
    const double css = qaa * L(2, 2) + 2 * qab * L(1, 3) + qbb * L(0, 4);
    Z[0] = cpp; Z[1] = cps; Z[2] = cps; Z[3] = css;

    // local = 2 S'GS + 2 (X2 + X2') + T'ZT + c1 pAA + c2 sAA, X2 = S'MT
    double* lo = local.begin();
    addXtwY(lo, S.data(), G.data(), S.data(), K2, K2, q, q, 2.0, work);
    addXtwY(X2.data(), S.data(), M.data(), T.data(), K2, 2, q, q, 1.0, work);
    addXtwY(lo, T.data(), Z.data(), T.data(), 2, 2, q, q, 1.0, work);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < q; ++i) {
        const std::size_t o = r + (std::size_t) n * (i + (std::size_t) q * j);
        lo[i + (std::size_t) q * j] += c1 * pAAp[o] + c2 * sAAp[o];
      }
    // N = 2 S'GC + 2 (S'Mc + T'tM C) + T'Zc + (c1 Wp + c2 Ws)'
    double* No = N.begin();
    addXtwY(No, S.data(), G.data(), C.data(), K2, K2, q, K, 2.0, work);
    addXtwY(No, S.data(), M.data(), cm.data(), K2, 2, q, K, 2.0, work);
    addXtwY(No, T.data(), tM.data(), C.data(), 2, K2, q, K, 2.0, work);
    addXtwY(No, T.data(), Z.data(), cm.data(), 2, 2, q, K, 1.0, work);
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < q; ++j) {
        const std::size_t o = r + (std::size_t) n * (k + (std::size_t) K * j);
        No[j + (std::size_t) q * k] += c1 * Wp[o] + c2 * Ws[o];
      }
    // Q = 2 C'GC + 2 (Q2 + Q2') + c'Zc + c1 Vp + c2 Vs, Q2 = C'Mc
    double* Qo = Q.begin();
    addXtwY(Qo, C.data(), G.data(), C.data(), K2, K2, K, K, 2.0, work);
    addXtwY(Q2.data(), C.data(), M.data(), cm.data(), K2, 2, K, K, 1.0, work);
    addXtwY(Qo, cm.data(), Z.data(), cm.data(), 2, 2, K, K, 1.0, work);
    for (int l = 0; l < K; ++l)
      for (int k = 0; k < K; ++k) {
        const std::size_t o = r + (std::size_t) n * (k + (std::size_t) K * l);
        Qo[k + (std::size_t) K * l] += c1 * Vp[o] + c2 * Vs[o];
      }

    // dH: D over theta (x = Wp, y = Ws, loss terms over pA, sA), Gam over the etas (Vp, Vs; a, b)
    for (int k = 0; k < K; ++k) { V[k] = Lpp * ar[k] + Lps * br[k]; Wv[k] = Lps * ar[k] + Lss * br[k]; }
    auto dH = [&](double* out, const double* X, const double* Y, const double* P, const double* Sg, int m,
                  bool rowsOfXY) {
      for (int j = 0; j < m; ++j) {
        const double pj = rowsOfXY ? P[r + (std::size_t) n * j] : P[j];
        const double sj = rowsOfXY ? Sg[r + (std::size_t) n * j] : Sg[j];
        const double dpp = L(3, 0) * pj + L(2, 1) * sj;
        const double dps = L(2, 1) * pj + L(1, 2) * sj;
        const double dss = L(1, 2) * pj + L(0, 3) * sj;
        for (int l = 0; l < K; ++l) {
          const double xl = X[r + (std::size_t) n * (l + (std::size_t) K * j)];
          const double yl = Y[r + (std::size_t) n * (l + (std::size_t) K * j)];
          for (int k = 0; k < K; ++k) {
            const double xk = X[r + (std::size_t) n * (k + (std::size_t) K * j)];
            const double yk = Y[r + (std::size_t) n * (k + (std::size_t) K * j)];
            out[k + (std::size_t) K * l + (std::size_t) KK * j] +=
                V[l] * xk + V[k] * xl + Wv[l] * yk + Wv[k] * yl +
                dpp * ar[k] * ar[l] + dps * (ar[k] * br[l] + br[k] * ar[l]) + dss * br[k] * br[l];
          }
        }
      }
    };
    dH(D.begin(), Wp.begin(), Ws.begin(), pA.begin(), sA.begin(), q, true);
    dH(Gam.begin(), Vp.begin(), Vs.begin(), ar.data(), br.data(), K, false);
  }

  if (!hessian) return List::create(_["g"] = g, _["gE"] = gE);
  for (int j = 0; j < q; ++j)
    for (int i = 0; i < q; ++i)
      local(i, j) += 2.0 * (X2[i + (std::size_t) q * j] + X2[j + (std::size_t) q * i]);
  for (int l = 0; l < K; ++l)
    for (int k = 0; k < K; ++k)
      Q(k, l) += 2.0 * (Q2[k + (std::size_t) K * l] + Q2[l + (std::size_t) K * k]);
  return List::create(_["g"] = g, _["gE"] = gE, _["local"] = local, _["N"] = N, _["Q"] = Q, _["D"] = D,
                      _["Gam"] = Gam);
}
