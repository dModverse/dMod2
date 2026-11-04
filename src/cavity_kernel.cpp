// The cavities of expectation consistent inference (.sparseExpectation()): per group the Gaussian
// N(mu, V) with V = (P + diag(Lam))^-1 and mu = V (P m + gam), its marginal precisions
// t = 1 / diag(V) - Lam and centres n = (mu / diag(V) - gam) / t. V by an LU factor with the
// condition test of solve() in R.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

// [[Rcpp::export]]
List ec_cavity(List P, NumericMatrix Lam, NumericMatrix gam, NumericMatrix Pm, CharacterVector stems,
               List dimnames) {
  const int G = Lam.nrow(), K = Lam.ncol();
  NumericMatrix tt(G, K), nn(G, K), mu(G, K);
  List cov(G);
  std::vector<double> A((std::size_t) K * K), work(4 * (std::size_t) K);
  std::vector<int> ipiv(K), iwork(K);
  const double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < G; ++i) {
    NumericMatrix Pi = P[i];
    if (Pi.nrow() != K || Pi.ncol() != K) stop("ec_cavity: a precision matrix is not K x K.");
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < K; ++j) A[j + (std::size_t) K * k] = Pi(j, k) + (j == k ? Lam(i, j) : 0.0);
    double anorm = 0.0;
    for (int k = 0; k < K; ++k) {
      double s = 0.0;
      for (int j = 0; j < K; ++j) s += std::fabs(A[j + (std::size_t) K * k]);
      anorm = std::max(anorm, s);
    }
    int info = 0;
    F77_CALL(dgetrf)(&K, &K, A.data(), &K, ipiv.data(), &info);
    if (info > 0) stop("Lapack routine dgesv: system is exactly singular: U[%d,%d] = 0", info, info);
    double rcond = 0.0;
    F77_CALL(dgecon)("1", &K, A.data(), &K, &anorm, &rcond, work.data(), iwork.data(), &info FCONE);
    if (rcond < eps) stop("system is computationally singular: reciprocal condition number = %g", rcond);
    NumericMatrix V(K, K);
    for (int j = 0; j < K; ++j) V(j, j) = 1.0;
    F77_CALL(dgetrs)("N", &K, &K, A.data(), &K, ipiv.data(), V.begin(), &K, &info FCONE);
    for (int j = 0; j < K; ++j) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += V(j, k) * (Pm(i, k) + gam(i, k));
      mu(i, j) = s;
    }
    for (int j = 0; j < K; ++j) {
      tt(i, j) = 1.0 / V(j, j) - Lam(i, j);
      nn(i, j) = (mu(i, j) / V(j, j) - gam(i, j)) / tt(i, j);
    }
    V.attr("dimnames") = List::create(stems, stems);
    cov[i] = V;
  }
  tt.attr("dimnames") = dimnames;
  nn.attr("dimnames") = dimnames;
  mu.attr("dimnames") = dimnames;
  return List::create(Named("t") = tt, Named("n") = nn, Named("mu") = mu, Named("cov") = cov);
}
