// Sparse-grid quadrature for adaptive NLME marginal integration.
//
// Smolyak combination of nested Genz-Keister rules for the weight exp(-z^2) (1D total mass
// sqrt(pi)).
//
// Conventions:
//   - Output `nodes` is z-space, batch-first [B, K] (column-major
//     NumericMatrix; matches dMod's batch-first invariant for sensitivities).
//   - Genz-Keister: 1D level l -> the smallest nested rule exact to degree
//     2l - 1; A(L, K) = sum_{|i| <= L} D^{i_1} (x) ... (x) D^{i_K} with the
//     difference rules D^l = U^l - U^{l-1}, so the grids of successive levels
//     are nested and A(L, K) is exact for polynomials of total degree 2(L-K)+1.
//   - The R wrapper folds the z -> eta change of variable
//     (eta = eta_hat + sqrt(2) * L^{-T} z, L = chol(H_i)) into log_weights.

#include <Rcpp.h>

#include <vector>
#include <map>
#include <cmath>

using namespace Rcpp;


// Append to `out` every length-K positive-integer multi-index with sum q.
static void enumerate_compositions(int K, int q,
                                   std::vector<int>& cur,
                                   int idx,
                                   int remaining,
                                   std::vector<std::vector<int>>& out) {
  if (idx == K - 1) {
    if (remaining >= 1) {
      cur[idx] = remaining;
      out.push_back(cur);
    }
    return;
  }
  int hi = remaining - (K - 1 - idx);
  for (int v = 1; v <= hi; ++v) {
    cur[idx] = v;
    enumerate_compositions(K, q, cur, idx + 1, remaining - v, out);
  }
}


// Genz-Keister rules for the standard normal weight (Genz and Keister 1996), nested: 1, 3, 9,
// 19 and 35 nodes, exact to degree 1, 5, 15, 29 and 51. Nodes x >= 0 of each symmetric rule and
// their weights; computed in 300-bit arithmetic as Kronrod-Patterson extensions.
static const std::vector<std::vector<double>> gk_x = {
  {0},
  {0, 1.7320508075688772},
  {0, 0.74109534999454085, 1.7320508075688772, 2.8612795760570582, 4.1849560176727323},
  {0, 0.74109534999454085, 1.230423634027306, 1.7320508075688772, 2.5960831150492023,
   2.8612795760570582, 3.2053337944991944, 4.1849560176727323, 5.1870160399136562,
   6.3633944943363696},
  {0, 0.24899229757996061, 0.74109534999454085, 1.230423634027306, 1.7320508075688772,
   2.2336260616769414, 2.5960831150492023, 2.8612795760570582, 3.2053337944991944,
   3.6353185190372783, 4.1849560176727323, 4.7364330859522967, 5.1870160399136562,
   5.6981777684881099, 6.3633944943363696, 7.1221067008046166, 7.9807717985905606,
   9.0169397898903032}};
static const std::vector<std::vector<double>> gk_w = {
  {1},
  {0.66666666666666663, 0.16666666666666666},
  {0.25396825396825395, 0.27007432957793787, 0.094850948509485097, 0.0079963254708935328,
   9.4269457556517484e-05},
  {0.30346719985420589, 0.20832499164960888, 0.061151730125247675, 0.064096054686807583,
   0.018085234254798452, -0.0063372247933737354, 0.0028848804365067511, 6.012336945984782e-05,
   6.094808731468983e-07, 8.6296846022298859e-10},
  {0.00051489450806878432, 0.19176011588804442, 0.14807083115521599, 0.092364726716986312,
   0.045273685465150516, 0.015673473751851151, 0.0031554462691875639, 0.0023113452403522102,
   0.00081895392750226486, 0.00027524214116785158, 3.5729348198975102e-05,
   2.7342206801187829e-06, 2.4676421345798077e-07, 2.1394194479561105e-08,
   4.6011760348656186e-10, 3.097222357606316e-12, 5.4500412650636897e-15,
   1.0541326582333341e-18}};
static const int gk_degree[] = {1, 5, 15, 29, 51};

// Index of the smallest Genz-Keister rule exact to degree 2l - 1, -1 beyond the table.
static int gk_rule(int l) {
  for (int r = 0; r < 5; ++r) if (gk_degree[r] >= 2 * l - 1) return r;
  return -1;
}

// Difference rule D^l = U^l - U^{l-1} of level l on the nodes z = x / sqrt(2) of U^l, weights
// sqrt(pi) w; empty when both levels share their rule.
static void gk_difference(int l, std::vector<double>& z, std::vector<double>& w) {
  z.clear(); w.clear();
  int r = gk_rule(l), rp = l > 1 ? gk_rule(l - 1) : -1;
  if (r == rp) return;
  const double s2 = std::sqrt(2.0), sp = std::sqrt(M_PI);
  const std::vector<double>& x = gk_x[r];
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    for (std::size_t j = 0; j < x.size(); ++j) {
      if (j == 0 && sgn < 0) continue;
      double wp = 0.0;
      if (rp >= 0) {
        const std::vector<double>& xp = gk_x[rp];
        for (std::size_t k = 0; k < xp.size(); ++k) if (xp[k] == x[j]) wp = gk_w[rp][k];
      }
      double d = gk_w[r][j] - wp;
      if (d == 0.0) continue;
      z.push_back(sgn * x[j] / s2);
      w.push_back(sp * d);
    }
  }
}

// The Smolyak grid of `level` in K dimensions (K <= level <= K + 25) for the weight exp(-z'z)
// from nested Genz-Keister rules, exact to total degree 2 (level - K) + 1: nodes and signed weights.
// [[Rcpp::export]]
List sparseGridGK(int K, int level) {
  if (K < 1) Rcpp::stop("sparseGridGK: K must be >= 1.");
  if (level < K)
    Rcpp::stop("sparseGridGK: level (%d) is below dimension K (%d).", level, K);
  if (gk_rule(level - K + 1) < 0)
    Rcpp::stop("sparseGridGK: level (%d) above K + 25 in K = %d.", level, K);
  int top = level - K + 1;
  std::vector<std::vector<double>> dz(top + 1), dw(top + 1);
  for (int l = 1; l <= top; ++l) gk_difference(l, dz[l], dw[l]);

  auto key_of = [](const std::vector<double>& z) {
    std::vector<long long> key(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) key[i] = static_cast<long long>(std::llround(z[i] * 1e12));
    return key;
  };
  std::map<std::vector<long long>, std::pair<std::vector<double>, double>> grid;
  std::vector<int> cur(K, 0);
  std::vector<std::vector<int>> midx;
  for (int q = K; q <= level; ++q) {
    midx.clear();
    enumerate_compositions(K, q, cur, 0, q, midx);
    for (const auto& i : midx) {
      bool empty = false;
      for (int d = 0; d < K; ++d) if (dz[i[d]].empty()) { empty = true; break; }
      if (empty) continue;
      std::vector<int> idx(K, 0);
      while (true) {
        std::vector<double> z(K);
        double w = 1.0;
        for (int d = 0; d < K; ++d) { z[d] = dz[i[d]][idx[d]]; w *= dw[i[d]][idx[d]]; }
        auto key = key_of(z);
        auto it = grid.find(key);
        if (it == grid.end()) grid[key] = std::make_pair(z, w);
        else it->second.second += w;
        int d = K - 1;
        while (d >= 0) {
          if (++idx[d] < static_cast<int>(dz[i[d]].size())) break;
          idx[d] = 0;
          --d;
        }
        if (d < 0) break;
      }
    }
  }
  std::vector<const std::pair<std::vector<double>, double>*> kept;
  for (const auto& kv : grid) if (kv.second.second != 0.0) kept.push_back(&kv.second);
  int B = static_cast<int>(kept.size());
  NumericMatrix nodes(B, K);
  NumericVector weights(B);
  for (int b = 0; b < B; ++b) {
    for (int d = 0; d < K; ++d) nodes(b, d) = kept[b]->first[d];
    weights[b] = kept[b]->second;
  }
  return List::create(_["nodes"] = nodes, _["weights"] = weights, _["K"] = K, _["level"] = level);
}
