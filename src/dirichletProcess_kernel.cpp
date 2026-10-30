// Cluster random effects of one symbol: the exact sum over all set partitions of the groups and
// the anchor (the population value, deviation 0).
//
// Group g enters through its quadratic factor exp(-a_g / 2 (eta_g - m_g)^2). A partition puts
// the groups and the anchor into blocks; the block of the anchor has deviation 0, every other
// block one common deviation v ~ N(0, tau^2). The prior of a partition is the Ewens distribution
// of the Chinese restaurant process with concentration c and the anchor seated first,
//
//   p(partition) = c^(K - 1) prod_b (n_b - 1)! / prod_{i = 1..G} (c + i),
//
// K the number of blocks, n_b their sizes with the anchor counted. With A, B, C the sums of
// a, a m and a m^2 over the groups of a block, its integral is exp(-C / 2) for the anchor
// block and (1 + tau^2 A)^(-1/2) exp(-C / 2 + B^2 / (2 (A + 1 / tau^2))) for a slab block.
//
// The sum over partitions is a subset recursion: with F(S) the sum over the partitions of the
// set S of the products of their block weights w,
//
//   F(S) = sum_{T in S, T holds the lowest element of S} w(T) F(S \ T),
//
// O(3^n) for n = G + 1 elements. The probability that T is a block is w(T) F(U \ T) / F(U),
// U the full set; posterior moments, M-step statistics and the derivatives of log I by
// Fisher's identity follow from these block probabilities. The most probable partition is
// the same recursion with max for sum.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

// running log-sum-exp
struct LogSum {
  double m = NEG_INF, s = 0.0;
  void add(double v) {
    if (v == NEG_INF) return;
    if (v > m) { s = s * std::exp(m - v) + 1.0; m = v; } else s += std::exp(v - m);
  }
  double value() const { return m == NEG_INF ? NEG_INF : m + std::log(s); }
};

// labels of a partition given per set its chosen block: 0 for the anchor block, 1, 2, ... for
// the slab blocks in the order of their first group
void backtrack(const std::vector<unsigned>& choice, unsigned full, int G, unsigned anchor,
               IntegerVector& labels, std::vector<unsigned>& blocks) {
  unsigned S = full;
  blocks.clear();
  while (S) { unsigned T = choice[S]; blocks.push_back(T); S &= ~T; }
  std::vector<int> label(G, -1);
  int next = 0;
  for (int g = 0; g < G; ++g) {
    if (label[g] >= 0) continue;
    for (unsigned T : blocks) {
      if (!(T & (1u << g))) continue;
      int l = (T & anchor) ? 0 : ++next;
      for (int h = 0; h < G; ++h) if (T & (1u << h)) label[h] = l;
    }
  }
  for (int g = 0; g < G; ++g) labels[g] = label[g];
}

}  // namespace

// [[Rcpp::export]]
List dirichlet_process_partitions(NumericVector a, NumericVector m, double concentration, double tau) {
  const int G = a.size();
  if (G > 14) stop("em: dirichletProcess() sums over all partitions of at most 14 groups; more groups are "
                   "fitted by method = \"saem\".");
  const int n = G + 1;
  const unsigned N = 1u << n, full = N - 1u, anchor = 1u << G;
  const double tau2 = tau * tau, lconc = std::log(concentration);

  // block sums and weights
  std::vector<double> A(N, 0.0), B(N, 0.0), C(N, 0.0), lw(N, NEG_INF);
  std::vector<int> size(N, 0);
  for (unsigned T = 1; T < N; ++T) {
    unsigned low = T & (~T + 1u), rest = T ^ low;
    int i = __builtin_ctz(low);
    A[T] = A[rest]; B[T] = B[rest]; C[T] = C[rest]; size[T] = size[rest] + 1;
    if (i < G) { A[T] += a[i]; B[T] += a[i] * m[i]; C[T] += a[i] * m[i] * m[i]; }
    const double prior = std::lgamma(static_cast<double>(size[T]));
    if (T & anchor) {
      lw[T] = -C[T] / 2.0 + prior;
    } else {
      const double P = A[T] + 1.0 / tau2;
      lw[T] = -0.5 * std::log1p(tau2 * A[T]) - C[T] / 2.0 + B[T] * B[T] / (2.0 * P) + prior + lconc;
    }
  }

  // F(S) by the subset recursion, the most probable partition of S by the same recursion with max
  std::vector<double> logF(N, NEG_INF), best(N, NEG_INF);
  std::vector<unsigned> choice(N, 0u);
  logF[0] = 0.0; best[0] = 0.0;
  for (unsigned S = 1; S < N; ++S) {
    unsigned low = S & (~S + 1u), rest = S ^ low;
    LogSum acc;
    for (unsigned sub = rest;; sub = (sub - 1u) & rest) {
      unsigned T = sub | low, R = S ^ T;
      acc.add(lw[T] + logF[R]);
      double vb = lw[T] + best[R];
      if (vb > best[S]) { best[S] = vb; choice[S] = T; }
      if (sub == 0u) break;
    }
    logF[S] = acc.value();
  }

  double norm = 0.0, dnorm = 0.0;
  for (int i = 1; i <= G; ++i) { norm += std::log(concentration + i); dnorm += 1.0 / (concentration + i); }
  const double logI = logF[full] - norm;

  // block probabilities and the moments
  NumericVector mean(G), cen2(G), spike(G);
  double nfree = 0.0, slab2 = 0.0;
  for (unsigned T = 1; T < N; ++T) {
    double p = std::exp(lw[T] + logF[full ^ T] - logF[full]);
    if (p == 0.0) continue;
    double mu = 0.0, var = 0.0;
    if (!(T & anchor)) {
      const double P = A[T] + 1.0 / tau2;
      mu = B[T] / P; var = 1.0 / P;
      nfree += p; slab2 += p * (var + mu * mu);
    }
    for (int g = 0; g < G; ++g) {
      if (!(T & (1u << g))) continue;
      mean[g] += p * mu;
      cen2[g] += p * (var + (mu - m[g]) * (mu - m[g]));
      if (T & anchor) spike[g] += p;
    }
  }
  NumericVector dm(G), da(G);
  for (int g = 0; g < G; ++g) { dm[g] = a[g] * (mean[g] - m[g]); da[g] = -cen2[g] / 2.0; }

  // the most probable partition and the posterior means of its block values
  IntegerVector labels(G);
  NumericVector values(G);
  std::vector<unsigned> blocks;
  backtrack(choice, full, G, anchor, labels, blocks);
  for (unsigned T : blocks) {
    if (T & anchor) continue;
    double v = B[T] / (A[T] + 1.0 / tau2);
    for (int g = 0; g < G; ++g) if (T & (1u << g)) values[g] = v;
  }

  return List::create(_["logI"] = logI, _["dm"] = dm, _["da"] = da, _["mean"] = mean, _["cen2"] = cen2,
                      _["spike"] = spike, _["dconcentration"] = nfree / concentration - dnorm,
                      _["dtau"] = -nfree / tau + slab2 / (tau2 * tau), _["free"] = nfree, _["slab2"] = slab2,
                      _["labels"] = labels, _["values"] = values);
}
