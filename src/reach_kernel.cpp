// Reachability in the dependency graph of a model, for the random effects that reach the data of
// their group (.etaReach()): edges from -> to between named nodes, and per source whether a
// breadth-first search from it meets a sink.

#include <Rcpp.h>

#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector graph_reach(CharacterVector from, CharacterVector to, CharacterVector sources,
                          CharacterVector sinks) {
  if (from.size() != to.size()) stop("graph_reach: `from` and `to` differ in length");
  std::unordered_map<std::string, int> id;
  auto number = [&](SEXP s) {
    auto it = id.emplace(CHAR(s), (int) id.size());
    return it.first->second;
  };
  std::vector<int> f(from.size()), t(to.size());
  for (R_xlen_t e = 0; e < from.size(); ++e) { f[e] = number(STRING_ELT(from, e)); t[e] = number(STRING_ELT(to, e)); }
  std::vector<int> src(sources.size()), snk(sinks.size());
  for (R_xlen_t i = 0; i < sources.size(); ++i) src[i] = number(STRING_ELT(sources, i));
  for (R_xlen_t i = 0; i < sinks.size(); ++i) snk[i] = number(STRING_ELT(sinks, i));
  const int nodes = id.size();
  std::vector<std::vector<int>> out(nodes);
  for (std::size_t e = 0; e < f.size(); ++e) out[f[e]].push_back(t[e]);
  std::vector<char> sink(nodes, 0);
  for (int s : snk) sink[s] = 1;
  LogicalVector hit(src.size());
  std::vector<char> seen(nodes);
  std::vector<int> queue;
  for (std::size_t i = 0; i < src.size(); ++i) {
    std::fill(seen.begin(), seen.end(), 0);
    queue.assign(1, src[i]);
    seen[src[i]] = 1;
    bool found = false;
    for (std::size_t q = 0; q < queue.size() && !found; ++q) {
      const int v = queue[q];
      if (sink[v]) { found = true; break; }
      for (int w : out[v])
        if (!seen[w]) { seen[w] = 1; queue.push_back(w); }
    }
    hit[i] = found;
  }
  return hit;
}
