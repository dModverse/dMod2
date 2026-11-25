// Dimension-adaptive sparse-grid quadrature (Gerstner and Griebel 2003) on the nested Genz-Keister rules
// of quadrature.cpp: the old indices of a downward closed set make the grid, the active ones the error.
// Terms sum w exp(l(z) - l(0) + z'z), the ratio of the integrand exp(l) to the Gaussian at the mode.

#include <Rcpp.h>

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace Rcpp;

namespace {

// Genz-Keister rules for the standard normal weight, as in quadrature.cpp: nodes x >= 0 and weights.
const std::vector<std::vector<double>> GK_X = {
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
const std::vector<std::vector<double>> GK_W = {
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
const int N_RULES = 5;
// The smallest Smolyak level of quadrature.cpp at which each rule enters one coordinate.
const int FIRST_LEVEL[N_RULES] = {1, 2, 4, 9, 16};

// one rule index r in 0..4 per coordinate (1, 3, 9, 19, 35 nodes), the term the tensor product of D_r
typedef std::vector<int> Index;
typedef std::vector<long long> Key;

// Difference rule D_r on z = x / sqrt(2) with weights sqrt(pi) (w_r - w_{r-1}).
void diffRule(int r, std::vector<double>& z, std::vector<double>& w) {
  z.clear(); w.clear();
  const double s2 = std::sqrt(2.0), sp = std::sqrt(M_PI);
  const std::vector<double>& x = GK_X[r];
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    for (std::size_t j = 0; j < x.size(); ++j) {
      if (j == 0 && sgn < 0) continue;
      double wp = 0.0;
      if (r > 0) for (std::size_t k = 0; k < GK_X[r - 1].size(); ++k) if (GK_X[r - 1][k] == x[j]) wp = GK_W[r - 1][k];
      double d = GK_W[r][j] - wp;
      if (d == 0.0) continue;
      z.push_back(sgn * x[j] / s2);
      w.push_back(sp * d);
    }
  }
}

struct Term {
  double delta = 0.0;
  bool computed = false;
  bool old = false;
};

class AdaptiveGrid {
 public:
  int K;
  std::vector<std::vector<double>> dz, dw;
  std::map<Index, Term> terms;
  std::map<Key, std::pair<std::vector<double>, double>> values;  // node -> (z, log value)
  std::vector<Index> waiting;                                      // proposed, not yet computed
  double l0 = NA_REAL;
  int maxLevel;                                                    // no index of a higher level

  AdaptiveGrid(int K_, int maxLevel_) : K(K_), dz(N_RULES), dw(N_RULES), maxLevel(maxLevel_) {
    for (int r = 0; r < N_RULES; ++r) diffRule(r, dz[r], dw[r]);
  }

  static Key keyOf(const std::vector<double>& z) {
    Key k(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) k[i] = std::llround(z[i] * 1e12);
    return k;
  }

  // calls f(z, w) for every node of the tensor term of `i`
  template <class F> void eachNode(const Index& i, F f) const {
    std::vector<int> pos(K, 0);
    std::vector<double> z(K);
    while (true) {
      double w = 1.0;
      for (int d = 0; d < K; ++d) { z[d] = dz[i[d]][pos[d]]; w *= dw[i[d]][pos[d]]; }
      f(z, w);
      int d = K - 1;
      while (d >= 0) {
        if (++pos[d] < static_cast<int>(dz[i[d]].size())) break;
        pos[d] = 0;
        --d;
      }
      if (d < 0) break;
    }
  }

  void add(const Index& i) {
    if (terms.count(i)) return;
    terms[i] = Term();
    waiting.push_back(i);
  }

  static int levelOf(const Index& i) {
    int lv = 0;
    for (int r : i) lv += FIRST_LEVEL[r];
    return lv;
  }

  // forward neighbors of `i` not yet in the set whose backward neighbors are all old become active;
  // above maxLevel they estimate the error but never become old
  bool addNeighbors(const Index& i) {
    bool added = false;
    for (int k = 0; k < K; ++k) {
      if (i[k] + 1 >= N_RULES) continue;
      Index j = i; ++j[k];
      if (terms.count(j)) continue;
      bool allOld = true;
      for (int d = 0; d < K && allOld; ++d) {
        if (j[d] == 0) continue;
        Index b = j; --b[d];
        auto it = terms.find(b);
        allOld = it != terms.end() && it->second.old;
      }
      if (allOld) { add(j); added = true; }
    }
    return added;
  }

  // the indices of `closed` old, their forward neighbors whose backward neighbors are all old active
  void start(const std::set<Index>& closed) {
    for (const Index& i : closed) { add(i); terms[i].old = true; }
    for (const Index& i : closed) addNeighbors(i);
  }

  // nodes of the waiting terms without a value
  NumericMatrix pending() const {
    std::set<Key> need;
    std::vector<std::vector<double>> out;
    for (const Index& i : waiting) eachNode(i, [&](const std::vector<double>& z, double) {
      Key k = keyOf(z);
      if (values.count(k) || need.count(k)) return;
      need.insert(k);
      out.push_back(z);
    });
    NumericMatrix m(out.size(), K);
    for (std::size_t b = 0; b < out.size(); ++b) for (int d = 0; d < K; ++d) m(b, d) = out[b][d];
    return m;
  }

  // log values of the nodes of pending(), in that order; the waiting terms are then computed
  void absorb(const NumericMatrix& z, const NumericVector& lv) {
    if (z.ncol() != K || lv.size() != z.nrow())
      Rcpp::stop("adaptive grid: absorb() needs a K-column node matrix and one log value per row.");
    for (int b = 0; b < z.nrow(); ++b) {
      std::vector<double> zb(K);
      double zz = 0.0;
      for (int d = 0; d < K; ++d) { zb[d] = z(b, d); zz += zb[d] * zb[d]; }
      values[keyOf(zb)] = std::make_pair(zb, lv[b]);
      if (zz == 0.0) l0 = lv[b];
    }
    if (ISNAN(l0)) Rcpp::stop("adaptive grid: the node at the mode has no value.");
    for (const Index& i : waiting) {
      double s = 0.0;
      eachNode(i, [&](const std::vector<double>& zn, double w) {
        auto it = values.find(keyOf(zn));
        if (it == values.end()) Rcpp::stop("adaptive grid: a pending node has no value.");
        double zz = 0.0;
        for (double v : zn) zz += v * v;
        s += w * std::exp(it->second.second - l0 + zz);
      });
      terms[i].delta = s;
      terms[i].computed = true;
    }
    waiting.clear();
  }

  // the number of evaluated nodes after the forward neighbors of `i` wait
  int prospective(const Index& i) const {
    std::set<Key> extra;
    for (int k = 0; k < K; ++k) {
      if (i[k] + 1 >= N_RULES) continue;
      Index j = i; ++j[k];
      if (terms.count(j)) continue;
      eachNode(j, [&](const std::vector<double>& z, double) {
        Key key = keyOf(z);
        if (!values.count(key)) extra.insert(key);
      });
    }
    return static_cast<int>(values.size() + extra.size());
  }

  // the active terms become old, the integral their last refinement: at a budget the change they
  // make stays the error estimate
  void accept() {
    for (auto& kv : terms) if (kv.second.computed) kv.second.old = true;
  }

  // the active term of largest |delta| becomes old and the forward neighbors whose backward neighbors
  // are all old wait: 1; 0 where no active term up to maxLevel is left, 2 where the evaluated nodes
  // would exceed `maxNodes`
  int refine(double maxNodes) {
    std::vector<std::pair<double, Index>> act;
    for (const auto& kv : terms)
      if (kv.second.computed && !kv.second.old && levelOf(kv.first) <= maxLevel)
        act.push_back({std::fabs(kv.second.delta), kv.first});
    std::sort(act.begin(), act.end(), [](const std::pair<double, Index>& a, const std::pair<double, Index>& b) {
      return a.first > b.first;
    });
    // the old set grows by the largest active term even where it has no forward neighbor left
    if (act.empty()) return 0;
    if (prospective(act.front().second) > maxNodes) return 2;
    for (const auto& a : act) {
      terms[a.second].old = true;
      if (addNeighbors(a.second)) return 1;
    }
    return 1;
  }

  // the integral of the old terms, the change by the active ones, their number, and the shell: the
  // terms of the old maximal indices, the change from the grid below where no active term is left
  List status() const {
    double total = 0.0, error = 0.0, shell = 0.0;
    int level = 0, active = 0;
    for (const auto& kv : terms) {
      if (!kv.second.computed) continue;
      if (!kv.second.old) { error += kv.second.delta; ++active; continue; }
      total += kv.second.delta;
      level = std::max(level, levelOf(kv.first));
      bool maximal = true;
      for (int d = 0; d < K && maximal; ++d) {
        if (kv.first[d] + 1 >= N_RULES) continue;
        Index j = kv.first; ++j[d];
        auto it = terms.find(j);
        if (it != terms.end() && it->second.old) maximal = false;
      }
      if (maximal) shell += kv.second.delta;
    }
    return List::create(_["total"] = total, _["error"] = std::fabs(error), _["change"] = error, _["shell"] = shell,
                        _["active"] = active, _["l0"] = l0,
                        _["nodes"] = static_cast<int>(values.size()), _["level"] = level);
  }

  // the combined grid of the old terms: nodes, signed weights for exp(-z'z) and log values
  List grid() const {
    std::map<Key, std::pair<std::vector<double>, double>> acc;
    for (const auto& kv : terms) {
      if (!kv.second.computed || !kv.second.old) continue;
      eachNode(kv.first, [&](const std::vector<double>& z, double w) {
        Key k = keyOf(z);
        auto it = acc.find(k);
        if (it == acc.end()) acc[k] = std::make_pair(z, w);
        else it->second.second += w;
      });
    }
    std::vector<const std::pair<const Key, std::pair<std::vector<double>, double>>*> kept;
    for (const auto& kv : acc) if (kv.second.second != 0.0) kept.push_back(&kv);
    int B = kept.size();
    NumericMatrix z(B, K);
    NumericVector w(B), lv(B);
    CharacterVector key(B);
    for (int b = 0; b < B; ++b) {
      const auto& e = *kept[b];
      std::string s;
      for (int d = 0; d < K; ++d) {
        z(b, d) = e.second.first[d];
        s += std::to_string(e.first[d]) + (d + 1 < K ? " " : "");
      }
      w[b] = e.second.second;
      lv[b] = values.at(e.first).second;
      key[b] = s;
    }
    return List::create(_["z"] = z, _["weights"] = w, _["logValues"] = lv, _["key"] = key);
  }

  // the index set as rows, with attributes "delta" (the terms) and "old"
  IntegerMatrix indices() const {
    IntegerMatrix m(terms.size(), K);
    NumericVector delta(terms.size());
    LogicalVector old(terms.size());
    int b = 0;
    for (const auto& kv : terms) {
      for (int d = 0; d < K; ++d) m(b, d) = kv.first[d];
      delta[b] = kv.second.delta;
      old[b] = kv.second.old;
      ++b;
    }
    m.attr("delta") = delta;
    m.attr("old") = old;
    return m;
  }
};

}  // namespace


// A grid of dimension K over the downward closure of the Smolyak set of `level` (every index with
// sum of FIRST_LEVEL at most `level`, the grid of sparseGridGK(K, level)) and the rule-index rows of
// `start`, such as the index set of an earlier grid; no index above `maxLevel` is added later.
// [[Rcpp::export]]
SEXP adaptiveGridNew(int K, int level, IntegerMatrix start, int maxLevel) {
  if (K < 1) Rcpp::stop("adaptiveGridNew: K must be >= 1.");
  if (start.nrow() > 0 && start.ncol() != K) Rcpp::stop("adaptiveGridNew: `start` needs K columns.");
  XPtr<AdaptiveGrid> g(new AdaptiveGrid(K, std::max(maxLevel, level)), true);
  std::vector<Index> seeds;
  for (int b = 0; b < start.nrow(); ++b) {
    Index i(K);
    for (int d = 0; d < K; ++d) i[d] = std::min(std::max(start(b, d), 0), N_RULES - 1);
    seeds.push_back(i);
  }
  // the Smolyak set, enumerated within the remaining level so that no other index is visited
  Index cur(K, 0);
  std::function<void(int, int)> smolyak = [&](int d, int budget) {
    if (d == K) { seeds.push_back(cur); return; }
    for (int r = 0; r < N_RULES && FIRST_LEVEL[r] + (K - d - 1) <= budget; ++r) {
      cur[d] = r;
      smolyak(d + 1, budget - FIRST_LEVEL[r]);
    }
    cur[d] = 0;
  };
  smolyak(0, std::max(level, K));
  // downward closure
  std::set<Index> closed;
  std::vector<Index> stack(seeds.begin(), seeds.end());
  while (!stack.empty()) {
    Index top = stack.back(); stack.pop_back();
    if (!closed.insert(top).second) continue;
    for (int d = 0; d < K; ++d) if (top[d] > 0) { Index j = top; --j[d]; stack.push_back(j); }
  }
  g->start(closed);
  return g;
}

// [[Rcpp::export]]
NumericMatrix adaptiveGridPending(SEXP grid) { return XPtr<AdaptiveGrid>(grid)->pending(); }

// [[Rcpp::export]]
void adaptiveGridAbsorb(SEXP grid, NumericMatrix z, NumericVector logValues) {
  XPtr<AdaptiveGrid>(grid)->absorb(z, logValues);
}

// [[Rcpp::export]]
int adaptiveGridRefine(SEXP grid, double maxNodes) { return XPtr<AdaptiveGrid>(grid)->refine(maxNodes); }

// [[Rcpp::export]]
void adaptiveGridAccept(SEXP grid) { XPtr<AdaptiveGrid>(grid)->accept(); }

// [[Rcpp::export]]
List adaptiveGridStatus(SEXP grid) { return XPtr<AdaptiveGrid>(grid)->status(); }

// [[Rcpp::export]]
List adaptiveGridNodes(SEXP grid) { return XPtr<AdaptiveGrid>(grid)->grid(); }

// [[Rcpp::export]]
IntegerMatrix adaptiveGridIndices(SEXP grid) { return XPtr<AdaptiveGrid>(grid)->indices(); }
