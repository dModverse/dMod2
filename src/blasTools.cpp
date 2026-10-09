// Batched dgemm on batch-first [B, M, N] arrays, the layout of cppDE
// sensitivities. The batch axis varies fastest in column-major storage, so
// per-batch slices are gathered into scratch except in bmm_lb.

#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <vector>

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif


// C = alpha * A * B + beta * C, column-major.
inline void dgemm_nn(int M, int N, int K,
                     double alpha,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double beta,
                     double* C, int ldc) {
  F77_CALL(dgemm)("N", "N", &M, &N, &K, &alpha,
           A, &lda, B, &ldb, &beta, C, &ldc FCONE FCONE);
}


// bmm_lb: [B,M,K] x [K,N] -> [B,M,N]. [B,M,K] and [B*M,K] share column-major
// memory, so one dgemm on the reshaped left operand covers all batches.
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  dgemm_nn(Bn * M, N, K,
           1.0,
           A.begin(), Bn * M,
           B.begin(), K,
           0.0,
           C.begin(), Bn * M);

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}


// bmm_rb: [M,K] x [B,K,N] -> [B,M,N]. Slice B[b,,] has row stride Bn, so it is
// gathered into contiguous scratch, multiplied by A and scattered into C[b,,].
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  std::vector<double> Bbuf((size_t) K * N);
  std::vector<double> Cbuf((size_t) M * N);

  const double* Bdata = B.begin();
  double* Cdata = C.begin();

  for (int b = 0; b < Bn; ++b) {
    // Gather B[b,,]: element (k,n) at Bdata[b + k*Bn + n*Bn*K]
    for (int n = 0; n < N; ++n) {
      for (int k = 0; k < K; ++k) {
        Bbuf[k + n * K] = Bdata[b + k * Bn + n * (size_t) Bn * K];
      }
    }

    dgemm_nn(M, N, K,
             1.0,
             A.begin(), M,
             Bbuf.data(), K,
             0.0,
             Cbuf.data(), M);

    // Scatter into C[b,,]: element (m,n) at Cdata[b + m*Bn + n*Bn*M]
    for (int n = 0; n < N; ++n) {
      for (int m = 0; m < M; ++m) {
        Cdata[b + m * Bn + n * (size_t) Bn * M] = Cbuf[m + n * M];
      }
    }
  }

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}


// bmm_bb: [B,M,K] x [B,K,N] -> [B,M,N], both operands gathered per slice.
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  std::vector<double> Abuf((size_t) M * K);
  std::vector<double> Bbuf((size_t) K * N);
  std::vector<double> Cbuf((size_t) M * N);

  const double* Adata = A.begin();
  const double* Bdata = B.begin();
  double* Cdata = C.begin();

  for (int b = 0; b < Bn; ++b) {
    // Gather A[b,,]
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        Abuf[m + k * M] = Adata[b + m * Bn + k * (size_t) Bn * M];
      }
    }
    // Gather B[b,,]
    for (int n = 0; n < N; ++n) {
      for (int k = 0; k < K; ++k) {
        Bbuf[k + n * K] = Bdata[b + k * Bn + n * (size_t) Bn * K];
      }
    }

    dgemm_nn(M, N, K,
             1.0,
             Abuf.data(), M,
             Bbuf.data(), K,
             0.0,
             Cbuf.data(), M);

    // Scatter into C[b,,]
    for (int n = 0; n < N; ++n) {
      for (int m = 0; m < M; ++m) {
        Cdata[b + m * Bn + n * (size_t) Bn * M] = Cbuf[m + n * M];
      }
    }
  }

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}
