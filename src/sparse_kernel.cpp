// The integral over the sparse random effects of em(), in natural parameters:
//
//   I = int exp(-eta' P_i eta / 2 + b_i' eta) prod_k p_k(eta_k) d eta,
//
// per group i the Gaussian part of the data term with precision P_i and linear term b_i, and the
// priors p_k of the spike-and-slab, Laplace and Dirichlet-process coordinates. Every factor is
// integrated against an exponential tilt exp(h eta - t eta^2 / 2); a flat tilt t = 0 is no special
// case. The coupling of the factors through P is integrated by expectation consistent (Opper and
// Winther 2005) as a convergent double loop (Heskes, Opper, Wiegerinck, Winther and Zoeter 2005):
// with the statistics phi = (eta, -eta^2 / 2), the cavities lambda, the sites gamma = nu - lambda
// of the Gaussian part and nu the natural parameters of the marginals,
//
//   log I = A_q(lambda) + A_r(gamma) - A_s(nu),
//
// the outer loop sets nu to the current marginals, the inner one minimises the convex
// Psi(lambda) = A_q(lambda) + A_r(nu - lambda) by Newton with a backtracking line search. At the
// fixed point d log I / d b_i = E_r eta, d log I / d P_i = -E_r[eta eta'] / 2 and the derivatives in
// the hyperparameters are those of the factors at their cavities. Where P couples no two factors,
// the spike-and-slab coordinates of a group counting as one, the integral is the product of the
// factors and exact. Groups that no Dirichlet-process stem joins are integrated one by one.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <cstdint>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
#define FCONE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();
const double LOG_2PI = std::log(2.0 * M_PI);

double logSumExp(double a, double b) {
  if (a == NEG_INF) return b;
  if (b == NEG_INF) return a;
  const double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::fabs(a - b)));
}

// log erfcx(x) = x^2 + log erfc(x), stable for every real x; beyond 10 by the asymptotic series
double logErfcx(double x) {
  if (x > 10.0) {
    const double t = 1.0 / (2.0 * x * x);
    double term = 1.0, s = 0.0;
    for (int k = 1; k <= 14; ++k) { term = -term * (2 * k - 1) * t; s += term; }
    return std::log1p(s) - std::log(x * std::sqrt(M_PI));
  }
  return M_LN2 + x * x + R::pnorm(-x * M_SQRT2, 0.0, 1.0, 1, 1);
}

// The half integral int_0^inf u^k exp(-u^2 / 2 - alpha u) du for k = 0: its log, and the moments
// E u^k, k = 1..4, of the normalised density; by the recurrence I_k + alpha I_{k-1} = (k-1) I_{k-2}
// forward for small alpha and by its continued fraction for the ratios I_k / I_{k-1} beyond.
struct Half { double logI0, e1, e2, e3, e4; };

Half half(double alpha) {
  Half h;
  h.logI0 = 0.5 * std::log(M_PI / 2.0) + logErfcx(alpha / M_SQRT2);
  if (alpha <= 2.0) {
    const double inv = std::exp(-h.logI0);
    h.e1 = inv - alpha;
    h.e2 = 1.0 - alpha * h.e1;
    h.e3 = 2.0 * h.e1 - alpha * h.e2;
    h.e4 = 3.0 * h.e2 - alpha * h.e3;
    return h;
  }
  double r = 0.0, r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0;
  for (int k = 300; k >= 2; --k) {
    r = (k - 1) / (alpha + r);
    if (k - 1 == 4) r4 = r;
    if (k - 1 == 3) r3 = r;
    if (k - 1 == 2) r2 = r;
    if (k - 1 == 1) r1 = r;
  }
  h.e1 = r1; h.e2 = r1 * r2; h.e3 = h.e2 * r3; h.e4 = h.e3 * r4;
  return h;
}

// A coordinate factor against its tilt: log Z, the raw moments m1..m4 and what the reports and
// the M-steps read.
struct Coord {
  bool ok = false;
  double logZ = 0, m1 = 0, m2 = 0, m3 = 0, m4 = 0, var = 0;
  double incl = 0, abs = 0, value = NA_REAL, dprob = 0, dtau = 0, dlambda = 0;
  int selected = 0;
};

// spike-and-slab: (1 - p) delta_0 + p N(0, tau^2); needs t > -1 / tau^2
Coord slab1(double h, double t, double p, double tau) {
  Coord c;
  const double prec = t + 1.0 / (tau * tau);
  if (!(prec > 0)) return c;
  const double logZs = -0.5 * std::log(tau * tau * prec) + h * h / (2.0 * prec);
  c.logZ = logSumExp(std::log1p(-p), std::log(p) + logZs);
  const double w = std::exp(std::log(p) + logZs - c.logZ), ms = h / prec, vs = 1.0 / prec;
  c.m1 = w * ms;
  c.m2 = w * (vs + ms * ms);
  c.m3 = w * (ms * ms * ms + 3.0 * ms * vs);
  c.m4 = w * (ms * ms * ms * ms + 6.0 * ms * ms * vs + 3.0 * vs * vs);
  c.var = w * vs + w * (1.0 - w) * ms * ms;
  c.incl = w;
  c.dprob = w / p - (1.0 - w) / (1.0 - p);
  c.dtau = -w / tau + c.m2 / (tau * tau * tau);
  c.selected = w > 0.5;
  c.value = c.selected ? ms : 0.0;
  c.ok = std::isfinite(c.logZ);
  return c;
}

// Laplace: lambda / 2 exp(-lambda |eta|); needs t > 0
Coord laplace1(double h, double t, double lambda) {
  Coord c;
  if (!(t > 0)) return c;
  const double s = std::sqrt(t);
  const Half p = half((lambda - h) / s), n = half((lambda + h) / s);
  const double lse = logSumExp(p.logI0, n.logI0);
  c.logZ = std::log(lambda / 2.0) - 0.5 * std::log(t) + lse;
  const double wp = std::exp(p.logI0 - lse), wn = std::exp(n.logI0 - lse);
  c.m1 = (wp * p.e1 - wn * n.e1) / s;
  c.m2 = (wp * p.e2 + wn * n.e2) / t;
  c.m3 = (wp * p.e3 - wn * n.e3) / (t * s);
  c.m4 = (wp * p.e4 + wn * n.e4) / (t * t);
  c.abs = (wp * p.e1 + wn * n.e1) / s;
  const double mp = p.e1 / s, mn = -n.e1 / s;
  c.var = wp * ((p.e2 - p.e1 * p.e1) / t + (mp - c.m1) * (mp - c.m1)) +
          wn * ((n.e2 - n.e1 * n.e1) / t + (mn - c.m1) * (mn - c.m1));
  c.dlambda = 1.0 / lambda - c.abs;
  c.ok = std::isfinite(c.logZ);
  return c;
}

// running log-sum-exp
struct LogSum {
  double m = NEG_INF, s = 0.0;
  void add(double v) {
    if (v == NEG_INF) return;
    if (v > m) { s = s * std::exp(m - v) + 1.0; m = v; } else s += std::exp(v - m);
  }
  double value() const { return m == NEG_INF ? NEG_INF : m + std::log(s); }
};

// A Dirichlet-process symbol over the groups against their tilts (h_g, t_g): the sum over all
// set partitions of the groups and the anchor (deviation 0) with the Ewens prior, the anchor
// seated first; a block integral is 1 for the anchor and (1 + tau^2 A)^(-1/2) exp(B^2 / (2 (A +
// 1 / tau^2))) for a slab block, A and B the sums of t and h. O(3^n), n = G + 1.
struct Stem {
  bool ok = false;
  double logZ = 0, dconc = 0, dtau = 0, free = 0, slab2 = 0;
  std::vector<double> m1, m2, m3, m4, spike, value;
  std::vector<int> label;
  // with `cross`, G x G: E eta_g eta_g', E eta_g eta_g'^2, E eta_g^2 eta_g'^2
  std::vector<double> e11, e12, e22;
};

// The groups up to which the moments between groups enter the Newton matrix, O(G 3^(G + 1)).
const int CROSS_MAX = 10;

Stem dpStem(const std::vector<double>& h, const std::vector<double>& t, double conc, double tau,
            bool cross = false) {
  Stem st;
  const int G = h.size();
  if (G > 14) stop("em: dirichletProcess() sums over all partitions of at most 14 groups; more groups are "
                   "fitted by method = \"saem\".");
  const double tau2 = tau * tau;
  double neg = 1.0 / tau2;
  for (int g = 0; g < G; ++g) neg += std::min(t[g], 0.0);
  if (!(neg > 0)) return st;
  const int n = G + 1;
  const unsigned N = 1u << n, full = N - 1u, anchor = 1u << G;
  const double lconc = std::log(conc);
  std::vector<double> A(N, 0.0), B(N, 0.0), lw(N, NEG_INF);
  std::vector<int> size(N, 0);
  for (unsigned T = 1; T < N; ++T) {
    const unsigned low = T & (~T + 1u), rest = T ^ low;
    const int i = __builtin_ctz(low);
    A[T] = A[rest]; B[T] = B[rest]; size[T] = size[rest] + 1;
    if (i < G) { A[T] += t[i]; B[T] += h[i]; }
    const double prior = std::lgamma(static_cast<double>(size[T]));
    if (T & anchor) lw[T] = prior;
    else {
      const double P = A[T] + 1.0 / tau2;
      lw[T] = -0.5 * std::log(tau2 * P) + B[T] * B[T] / (2.0 * P) + prior + lconc;
    }
  }
  std::vector<double> logF(N, NEG_INF), best(N, NEG_INF);
  std::vector<unsigned> choice(N, 0u);
  logF[0] = 0.0; best[0] = 0.0;
  for (unsigned S = 1; S < N; ++S) {
    const unsigned low = S & (~S + 1u), rest = S ^ low;
    LogSum acc;
    for (unsigned sub = rest;; sub = (sub - 1u) & rest) {
      const unsigned T = sub | low, R = S ^ T;
      acc.add(lw[T] + logF[R]);
      const double vb = lw[T] + best[R];
      if (vb > best[S]) { best[S] = vb; choice[S] = T; }
      if (sub == 0u) break;
    }
    logF[S] = acc.value();
  }
  double norm = 0.0, dnorm = 0.0;
  for (int i = 1; i <= G; ++i) { norm += std::log(conc + i); dnorm += 1.0 / (conc + i); }
  st.logZ = logF[full] - norm;
  st.m1.assign(G, 0.0); st.m2.assign(G, 0.0); st.m3.assign(G, 0.0); st.m4.assign(G, 0.0);
  st.spike.assign(G, 0.0); st.value.assign(G, 0.0); st.label.assign(G, 0);
  for (unsigned T = 1; T < N; ++T) {
    const double p = std::exp(lw[T] + logF[full ^ T] - logF[full]);
    if (p == 0.0) continue;
    double mu = 0.0, v = 0.0;
    if (!(T & anchor)) {
      const double P = A[T] + 1.0 / tau2;
      mu = B[T] / P; v = 1.0 / P;
      st.free += p; st.slab2 += p * (v + mu * mu);
    }
    for (int g = 0; g < G; ++g) {
      if (!(T & (1u << g))) continue;
      st.m1[g] += p * mu;
      st.m2[g] += p * (v + mu * mu);
      st.m3[g] += p * (mu * mu * mu + 3.0 * mu * v);
      st.m4[g] += p * (mu * mu * mu * mu + 6.0 * mu * mu * v + 3.0 * v * v);
      if (T & anchor) st.spike[g] += p;
    }
  }
  st.dconc = st.free / conc - dnorm;
  st.dtau = -st.free / tau + st.slab2 / (tau2 * tau);
  if (cross && G <= CROSS_MAX) {
    // the moments of two groups: from a common slab block the moments of its value, from two
    // blocks the product of theirs, Sum_T p_T a_f(T) r(U \ T) with r(R) the expectation of a_h of
    // the block of g' over the partitions of R, by the same recursion as F
    std::vector<double> a1(N, 0.0), a2(N, 0.0);
    for (unsigned T = 1; T < N; ++T) {
      if (T & anchor) continue;
      const double P = A[T] + 1.0 / tau2, mu = B[T] / P;
      a1[T] = mu; a2[T] = 1.0 / P + mu * mu;
    }
    st.e11.assign(G * G, 0.0); st.e12.assign(G * G, 0.0); st.e22.assign(G * G, 0.0);
    std::vector<double> pT(N, 0.0);
    for (unsigned T = 1; T < N; ++T) pT[T] = std::exp(lw[T] + logF[full ^ T] - logF[full]);
    for (unsigned T = 1; T < N; ++T) {
      if ((T & anchor) || pT[T] == 0.0) continue;
      const double P = A[T] + 1.0 / tau2, mu = B[T] / P, v = 1.0 / P;
      const double v11 = v + mu * mu, v12 = mu * mu * mu + 3.0 * mu * v,
                   v22 = mu * mu * mu * mu + 6.0 * mu * mu * v + 3.0 * v * v;
      for (int g = 0; g < G; ++g) {
        if (!(T & (1u << g))) continue;
        for (int k = 0; k < G; ++k) {
          if (!(T & (1u << k))) continue;
          st.e11[g + G * k] += pT[T] * v11; st.e12[g + G * k] += pT[T] * v12; st.e22[g + G * k] += pT[T] * v22;
        }
      }
    }
    std::vector<double> r1(N), r2(N);
    for (int k = 0; k < G; ++k) {
      const unsigned bk = 1u << k;
      r1[0] = r2[0] = 0.0;
      for (unsigned S = 1; S < N; ++S) {
        double s1 = 0.0, s2 = 0.0;
        if (S & bk) {
          const unsigned low = S & (~S + 1u), rest = S ^ low;
          for (unsigned sub = rest;; sub = (sub - 1u) & rest) {
            const unsigned T = sub | low, R = S ^ T;
            const double q = std::exp(lw[T] + logF[R] - logF[S]);
            if (T & bk) { s1 += q * a1[T]; s2 += q * a2[T]; } else { s1 += q * r1[R]; s2 += q * r2[R]; }
            if (sub == 0u) break;
          }
        }
        r1[S] = s1; r2[S] = s2;
      }
      for (unsigned T = 1; T < N; ++T) {
        if ((T & anchor) || (T & bk) || pT[T] == 0.0) continue;
        const unsigned R = full ^ T;
        for (int g = 0; g < G; ++g) {
          if (!(T & (1u << g))) continue;
          st.e11[g + G * k] += pT[T] * a1[T] * r1[R];
          st.e12[g + G * k] += pT[T] * a1[T] * r2[R];
          st.e22[g + G * k] += pT[T] * a2[T] * r2[R];
        }
      }
    }
  }
  // the most probable partition: labels 0 for the anchor block, 1, 2, ... for the slab blocks in
  // the order of their first group, with the posterior means of their values
  std::vector<unsigned> blocks;
  for (unsigned S = full; S; S &= ~choice[S]) blocks.push_back(choice[S]);
  std::vector<int> lab(G, -1);
  int next = 0;
  for (int g = 0; g < G; ++g) {
    if (lab[g] >= 0) continue;
    for (unsigned T : blocks) {
      if (!(T & (1u << g))) continue;
      const int l = (T & anchor) ? 0 : ++next;
      const double v = (T & anchor) ? 0.0 : B[T] / (A[T] + 1.0 / tau2);
      for (int k = 0; k < G; ++k) if (T & (1u << k)) { lab[k] = l; st.value[k] = v; }
    }
  }
  st.label = lab;
  st.ok = std::isfinite(st.logZ);
  return st;
}

// dense helpers on column-major n x n matrices, by LAPACK
// the lower Cholesky factor of A in place; false where A is not positive definite
bool cholesky(std::vector<double>& A, int n) {
  int info = 0;
  F77_CALL(dpotrf)("L", &n, A.data(), &n, &info FCONE);
  return info == 0;
}

// X = A^-1 X for the factor L of A = L L' and the nrhs columns of X
void cholSolve(const std::vector<double>& L, int n, double* X, int nrhs = 1) {
  int info = 0;
  F77_CALL(dpotrs)("L", &n, &nrhs, L.data(), &n, X, &n, &info FCONE);
}

// x = A^-1 x by an LU factor with partial pivoting; false where a pivot falls below 1e-13 of the
// largest entry of A
bool luSolve(std::vector<double> A, int n, double* x) {
  double scale = 0.0;
  for (double a : A) scale = std::max(scale, std::fabs(a));
  if (!(scale > 0) || !std::isfinite(scale)) return false;
  std::vector<int> piv(n);
  int info = 0, one = 1;
  F77_CALL(dgetrf)(&n, &n, A.data(), &n, piv.data(), &info);
  if (info != 0) return false;
  for (int j = 0; j < n; ++j) if (!(std::fabs(A[j + j * n]) > 1e-13 * scale)) return false;
  F77_CALL(dgetrs)("N", &n, &one, A.data(), &n, piv.data(), x, &n, &info FCONE);
  return info == 0;
}

// The Gaussian part of one group under sites: Q = P + diag(lt), its inverse C and the mean
// C (b + lh), log |Q|, by a Cholesky factor of the Jacobi-scaled Q; false outside its domain.
struct Gauss { std::vector<double> C, mu; double logdet = 0; };

bool gauss(const std::vector<double>& P, const std::vector<double>& b, const std::vector<double>& lt,
           const std::vector<double>& lh, int n, Gauss& out) {
  std::vector<double> d(n), Q(n * n);
  for (int j = 0; j < n; ++j) {
    const double q = P[j + j * n] + lt[j];
    if (!(q > 0)) return false;
    d[j] = 1.0 / std::sqrt(q);
  }
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) Q[i + j * n] = d[i] * (P[i + j * n] + (i == j ? lt[j] : 0.0)) * d[j];
  if (!cholesky(Q, n)) return false;
  out.logdet = 0.0;
  for (int j = 0; j < n; ++j) out.logdet += 2.0 * std::log(Q[j + j * n]) - 2.0 * std::log(d[j]);
  int info = 0;
  F77_CALL(dpotri)("L", &n, Q.data(), &n, &info FCONE);
  if (info != 0) return false;
  out.C.assign(n * n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = j; i < n; ++i) out.C[i + j * n] = out.C[j + i * n] = d[i] * Q[i + j * n] * d[j];
  // the mean, refined once against the residual of Q mu = b + lh
  out.mu.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = 0; k < n; ++k) s += out.C[i + k * n] * (b[k] + lh[k]);
    out.mu[i] = s;
  }
  std::vector<double> res(n);
  for (int i = 0; i < n; ++i) {
    double s = b[i] + lh[i];
    for (int k = 0; k < n; ++k) s -= (P[i + k * n] + (i == k ? lt[k] : 0.0)) * out.mu[k];
    res[i] = s;
  }
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = 0; k < n; ++k) s += out.C[i + k * n] * res[k];
    out.mu[i] += s;
  }
  return true;
}

// The spike-and-slab coordinates of a group jointly on their precision Pf (K x K) with linear terms
// h: the sum over the subsets F of slab coordinates of closed-form Gaussian integrals, eta = 0 off
// F, with A_F = Pf_FF + diag(1 / tau^2). A depth-first search over the coordinates extends
// W = L^-1, L the Cholesky factor of A_F, by one row per included coordinate, O(2^K K^2) in all.
// As A_F^-1 = sum_r w_r w_r' over the rows of W and row r depends on the coordinates up to r only,
// an included coordinate adds w_r w_r' times the weight of the subsets below it; these rows and the
// subset means go in blocks into the second moment by dsyrk. A subtree whose subsets weigh less than
// exp(-PRUNE) of the largest found is left out, at most 2^K exp(-PRUNE) of Z: with S the subset and
// the coordinates still open, h_F' A_F^-1 h_F is at most that of S, and every later diagonal of L at
// least the conditional precision 1 / (A_S^-1)_jj; A_S^-1 follows the exclusions by rank-one
// downdates. The child of the larger bound goes first, so that large subsets are found early.
struct Joint {
  bool ok = false, over = false;  // over: the search left its budget of nodes
  double logZ = 0;
  std::vector<double> mean, cov, incl, slab2, dprob, dtau, value;
  std::vector<int> selected;
};

class SlabSum {
public:
  SlabSum(const std::vector<double>& Pf, const std::vector<double>& h, const std::vector<double>& p,
          const std::vector<double>& tau, int K, double budget)
      : budget(budget), K(K), Pf(Pf), h(h), W(K * K, 0.0), muStack((K + 1) * K, 0.0), F(K), z(K), buf(K * NB, 0.0),
        C(K * K, 0.0), mean(K, 0.0), incl(K, 0.0), best(K, 0.0), bestIn(K, 0), in(K, 0), l(K), t(K), v(K),
        Vs((K + 1) * K * K, 0.0), ms((K + 1) * K, 0.0), qs(K + 1, 0.0) {
    lin.resize(K); lout.resize(K); itau2.resize(K);
    for (int k = 0; k < K; ++k) {
      lin[k] = std::log(p[k]) - std::log(tau[k]);
      lout[k] = std::log1p(-p[k]);
      itau2[k] = 1.0 / (tau[k] * tau[k]);
    }
    // A^-1 of all coordinates and its mean and quadratic form, the bound at the root
    std::vector<double> A(K * K);
    for (int c = 0; c < K; ++c)
      for (int r = 0; r < K; ++r) A[r + c * K] = Pf[r + c * K] + (r == c ? itau2[r] : 0.0);
    int info = 0;
    F77_CALL(dpotrf)("L", &K, A.data(), &K, &info FCONE);
    if (info == 0) F77_CALL(dpotri)("L", &K, A.data(), &K, &info FCONE);
    bounded = info == 0;
    if (!bounded) return;
    for (int c = 0; c < K; ++c)
      for (int r = c; r < K; ++r) Vs[r + c * K] = Vs[c + r * K] = A[r + c * K];
    for (int r = 0; r < K; ++r) {
      double sm = 0.0;
      for (int c = 0; c < K; ++c) sm += Vs[r + c * K] * h[c];
      ms[r] = sm;
      qs[0] += h[r] * sm;
    }
  }

  // the number of nodes of the search, by Knuth's estimator over `probes` random paths with the
  // pruning bound of the subset a greedy path reaches first; the final bound only prunes more
  double size(int probes = 32) {
    if (!bounded) return std::ldexp(1.0, K + 1);
    probing = true;
    greedy = true;
    pathWeight = 1.0; estimate = 0.0;
    node(0, 0, 0.0, 0.0, 0.0, 0);
    greedy = false;
    double total = 0.0;
    for (int r = 0; r < probes; ++r) {
      pathWeight = 1.0; estimate = 1.0;
      node(0, 0, 0.0, 0.0, 0.0, 0);
      total += estimate;
    }
    probing = false;
    M = NEG_INF; bestLog = NEG_INF; failed = false; nodes = 0;
    return total / probes;
  }

  Joint run(const std::vector<double>& p, const std::vector<double>& tau) {
    Joint J;
    if (size() > budget) { J.over = true; return J; }
    node(0, 0, 0.0, 0.0, 0.0, 0);
    J.over = over;
    if (failed) return J;
    flush();
    J.logZ = M + std::log(S);
    J.mean.assign(K, 0.0); J.cov.assign(K * K, 0.0); J.incl.assign(K, 0.0); J.slab2.assign(K, 0.0);
    for (int k = 0; k < K; ++k) { J.mean[k] = mean[k] / S; J.incl[k] = incl[k] / S; }
    for (int c = 0; c < K; ++c)
      for (int r = c; r < K; ++r) {
        const double e = C[r + c * K] / S - J.mean[r] * J.mean[c];
        J.cov[r + c * K] = e; J.cov[c + r * K] = e;
      }
    J.dprob.assign(K, 0.0); J.dtau.assign(K, 0.0); J.value = best; J.selected = bestIn;
    for (int k = 0; k < K; ++k) {
      J.slab2[k] = J.cov[k + k * K] + J.mean[k] * J.mean[k];
      J.dprob[k] = J.incl[k] / p[k] - (1.0 - J.incl[k]) / (1.0 - p[k]);
      J.dtau[k] = -J.incl[k] / tau[k] + J.slab2[k] / (tau[k] * tau[k] * tau[k]);
    }
    J.ok = std::isfinite(J.logZ);
    return J;
  }

private:
  static const int NB = 256;
  const double budget;
  double nodes = 0;
  bool over = false;
  const int K;
  const std::vector<double>& Pf;
  const std::vector<double>& h;
  std::vector<double> lin, lout, itau2;
  // W column-major with leading dimension K, its row r over the positions 0..r; per depth the
  // subset mean over the coordinates; the coordinate and z = L^-1 h_F at every position
  std::vector<double> W, muStack;
  std::vector<int> F;
  std::vector<double> z;
  // the block of scaled vectors for dsyrk, the second moment (lower triangle), the mean and the
  // inclusion weights, all relative to exp(M), and the subset of largest weight
  std::vector<double> buf, C, mean, incl, best;
  std::vector<int> bestIn, in;
  std::vector<double> l, t, v;
  // per slot A_S^-1 (0 off S), A_S^-1 h_S and h_S' A_S^-1 h_S; slot k + 1 is S after excluding k
  std::vector<double> Vs, ms, qs;
  int nb = 0;
  double M = NEG_INF, S = 0.0, bestLog = NEG_INF;
  bool failed = false, bounded = false;
  static constexpr double PRUNE = 40.0;
  // the probes of size(): one path, the greedy one or a random one, and the running estimate
  bool probing = false, greedy = false;
  double pathWeight = 1.0, estimate = 0.0;
  std::uint64_t rng = 0x9e3779b97f4a7c15ULL;
  bool coin() {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return rng & 1ULL;
  }
  // bounds at the nodes with at least this many coordinates still open
  static const int OPEN = 3;

  // the sum over the open coordinates from j of the larger of exclusion and the bound of inclusion
  double openBound(int from, int slot) const {
    const double* V = Vs.data() + (std::size_t) slot * K * K;
    double b = 0.0;
    for (int j = from; j < K; ++j) b += std::max(lout[j], lin[j] + 0.5 * std::log(V[j + j * K]));
    return b;
  }
  // slot `to` as slot `from` with coordinate k excluded
  void excludeInto(int k, int from, int to) {
    const double* V = Vs.data() + (std::size_t) from * K * K;
    double* V1 = Vs.data() + (std::size_t) to * K * K;
    const double* mu = ms.data() + (std::size_t) from * K;
    double* mu1 = ms.data() + (std::size_t) to * K;
    const double vkk = V[k + k * K], f = mu[k] / vkk;
    for (int c = 0; c < K; ++c) {
      const double vc = V[k + c * K] / vkk;
      for (int r = 0; r < K; ++r) V1[r + c * K] = V[r + c * K] - V[r + k * K] * vc;
    }
    for (int r = 0; r < K; ++r) { mu1[r] = mu[r] - f * V[r + k * K]; V1[r + k * K] = V1[k + r * K] = 0.0; }
    mu1[k] = 0.0;
    qs[to] = qs[from] - mu[k] * f;
  }

  void flush() {
    if (!nb) return;
    const double one = 1.0;
    F77_CALL(dsyrk)("L", "N", &K, &nb, &one, buf.data(), &K, &one, C.data(), &K FCONE FCONE);
    nb = 0;
  }
  // the vector v times sqrt(exp(lw - M)) into the block
  void push(const std::vector<double>& x, double lw) {
    if (nb == NB) flush();
    const double s = std::exp(0.5 * (lw - M));
    double* col = buf.data() + (std::size_t) nb * K;
    for (int k = 0; k < K; ++k) col[k] = s * x[k];
    ++nb;
  }
  void rescale(double newM) {
    flush();
    const double f = std::exp(M - newM);
    if (std::isfinite(M)) {
      S *= f;
      for (double& x : C) x *= f;
      for (double& x : mean) x *= f;
      for (double& x : incl) x *= f;
    }
    M = newM;
  }
  // the log weight of the subsets below the decision on coordinate k, m coordinates included, the
  // bound of S in `slot`
  double node(int k, int m, double lw, double quad, double half, int slot) {
    if (failed) return NEG_INF;
    if (!probing && ++nodes > budget) { over = failed = true; return NEG_INF; }
    if (k == K) {
      if (probing) {
        const double lg = lw + 0.5 * quad - half;
        if (greedy && lg > M) M = lg;
        return lg;
      }
      const double lg = lw + 0.5 * quad - half;
      if (lg > M) rescale(lg);
      const double w = std::exp(lg - M);
      S += w;
      const double* mu = muStack.data() + (std::size_t) m * K;
      for (int c = 0; c < K; ++c) { mean[c] += w * mu[c]; v[c] = mu[c]; }
      push(v, lg);
      if (lg > bestLog) {
        bestLog = lg;
        best.assign(mu, mu + K);
        bestIn = in;
      }
      return lg;
    }
    double bIn = 0.0, bOut = 0.0;
    const bool bound = bounded && K - k >= OPEN;
    if (bound) {
      const double* V = Vs.data() + (std::size_t) slot * K * K;
      bIn = lw + lin[k] + 0.5 * std::log(V[k + k * K]) - half + openBound(k + 1, slot) + 0.5 * qs[slot];
      excludeInto(k, slot, k + 1);
      bOut = lw + lout[k] - half + openBound(k + 1, k + 1) + 0.5 * qs[k + 1];
    }
    const bool skipIn = bound && bIn < M - PRUNE, skipOut = bound && bOut < M - PRUNE;
    double out = NEG_INF, inside = NEG_INF;
    if (probing) {
      const int live = (int) !skipIn + (int) !skipOut;
      if (!live) return NEG_INF;
      pathWeight *= live;
      estimate += pathWeight;
      const bool goIn = greedy ? (!skipIn && (skipOut || bIn > bOut)) : (live == 2 ? coin() : !skipIn);
      return goIn ? include(k, m, lw, quad, half, slot)
                  : node(k + 1, m, lw + lout[k], quad, half, bound ? k + 1 : slot);
    }
    if (bIn > bOut) {
      if (!skipIn) inside = include(k, m, lw, quad, half, slot);
      if (!skipOut && !failed) out = node(k + 1, m, lw + lout[k], quad, half, bound ? k + 1 : slot);
    } else {
      if (!skipOut) out = node(k + 1, m, lw + lout[k], quad, half, bound ? k + 1 : slot);
      if (!skipIn && !failed) inside = include(k, m, lw, quad, half, slot);
    }
    if (failed) return NEG_INF;
    return logSumExp(out, inside);
  }

  // the subsets below with coordinate k included at position m
  double include(int k, int m, double lw, double quad, double half, int slot) {
    // the new row of L is l = W A_F,k, its diagonal d; the new row of W is (-l' W / d, 1 / d)
    const int inc = 1;
    for (int c = 0; c < m; ++c) l[c] = Pf[F[c] + k * K];
    if (m) F77_CALL(dtrmv)("L", "N", "N", &m, W.data(), &K, l.data(), &inc FCONE FCONE FCONE);
    double d2 = Pf[k + k * K] + itau2[k], lz = 0.0;
    for (int c = 0; c < m; ++c) { d2 -= l[c] * l[c]; lz += l[c] * z[c]; }
    if (!(d2 > 0)) { failed = true; return NEG_INF; }
    const double d = std::sqrt(d2), zk = (h[k] - lz) / d;
    std::copy(l.begin(), l.begin() + m, t.begin());
    if (m) F77_CALL(dtrmv)("L", "T", "N", &m, W.data(), &K, t.data(), &inc FCONE FCONE FCONE);
    for (int c = 0; c < m; ++c) W[m + c * K] = -t[c] / d;
    W[m + m * K] = 1.0 / d;
    F[m] = k; z[m] = zk; in[k] = 1;
    double* next = muStack.data() + (std::size_t) (m + 1) * K;
    std::copy(muStack.begin() + (std::size_t) m * K, muStack.begin() + (std::size_t) (m + 1) * K, next);
    for (int c = 0; c <= m; ++c) next[F[c]] += zk * W[m + c * K];
    const double inside = node(k + 1, m + 1, lw + lin[k], quad + zk * zk, half + std::log(d), slot);
    in[k] = 0;
    if (probing || failed || inside == NEG_INF) return inside;
    // row m of W, unchanged below, into the second moment with the weight of the subsets below
    std::fill(v.begin(), v.end(), 0.0);
    for (int c = 0; c <= m; ++c) v[F[c]] = W[m + c * K];
    push(v, inside);
    incl[k] += std::exp(inside - M);
    return inside;
  }
};

Joint slabJoint(const std::vector<double>& Pf, const std::vector<double>& h, const std::vector<double>& p,
                const std::vector<double>& tau, int K, double budget) {
  SlabSum sum(Pf, h, p, tau, K, budget);
  return sum.run(p, tau);
}

// The layout of one call: per group its active stems, per stem its kind and hyperparameters.
struct Layout {
  int G, K;
  std::vector<int> kind;  // 0 spike and slab, 1 Laplace, 2 Dirichlet process
  std::vector<double> prob, tau, lambda, conc;
  std::vector<std::vector<int>> act;  // the active stems of every group
  std::vector<std::vector<double>> P, b;  // per group on its active stems
};

// the prior variance of a coordinate, the start of its site
double priorVar(const Layout& L, int j) {
  if (L.kind[j] == 0) return L.prob[j] * L.tau[j] * L.tau[j];
  if (L.kind[j] == 1) return 2.0 / (L.lambda[j] * L.lambda[j]);
  return L.tau[j] * L.tau[j];
}

// Everything a call returns, on the full G x K layout.
struct Result {
  bool ok = false, converged = false;
  double logI = 0;
  int outer = 0, inner = 0;
  std::vector<double> mean, abs, value, spike, clusters, dprob, dtau, dlambda, dconc, gprob, gtau, glambda,
      free, groups, slab2, Lam, gam;
  std::vector<std::vector<double>> cov;  // per group K x K
};

void initResult(const Layout& L, Result& R) {
  const int GK = L.G * L.K;
  R.mean.assign(GK, 0.0); R.abs.assign(GK, 0.0); R.value.assign(GK, NA_REAL); R.spike.assign(GK, NA_REAL);
  R.clusters.assign(GK, NA_REAL); R.gprob.assign(GK, 0.0); R.gtau.assign(GK, 0.0); R.glambda.assign(GK, 0.0);
  R.Lam.assign(GK, 0.0); R.gam.assign(GK, 0.0);
  R.dprob.assign(L.K, 0.0); R.dtau.assign(L.K, 0.0); R.dlambda.assign(L.K, 0.0); R.dconc.assign(L.K, 0.0);
  R.free.assign(L.K, 0.0); R.groups.assign(L.K, 0.0); R.slab2.assign(L.K, 0.0);
  R.cov.assign(L.G, std::vector<double>(L.K * L.K, 0.0));
  // a coordinate without data keeps its prior
  for (int i = 0; i < L.G; ++i)
    for (int j = 0; j < L.K; ++j) {
      R.cov[i][j + j * L.K] = priorVar(L, j);
      if (L.kind[j] == 0) { R.spike[i + L.G * j] = 1.0 - L.prob[j]; R.clusters[i + L.G * j] = 0; R.value[i + L.G * j] = 0; }
      if (L.kind[j] == 1) R.abs[i + L.G * j] = 1.0 / L.lambda[j];
    }
}

// the reports of a coordinate factor at (group i, stem j)
void reportCoord(const Layout& L, int i, int j, const Coord& c, Result& R) {
  const int at = i + L.G * j;
  R.abs[at] = c.abs;
  if (L.kind[j] == 0) {
    R.value[at] = c.value; R.spike[at] = 1.0 - c.incl; R.clusters[at] = c.selected;
    R.gprob[at] += c.dprob; R.gtau[at] += c.dtau;
    R.dprob[j] += c.dprob; R.dtau[j] += c.dtau;
    R.free[j] += c.incl; R.groups[j] += 1; R.slab2[j] += c.m2;
  } else {
    R.value[at] = NA_REAL; R.spike[at] = NA_REAL; R.clusters[at] = NA_REAL;
    R.glambda[at] += c.dlambda; R.dlambda[j] += c.dlambda;
  }
}

void reportStem(const Layout& L, int j, const std::vector<int>& groups, const Stem& s, Result& R) {
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const int at = groups[g] + L.G * j;
    R.value[at] = s.value[g]; R.spike[at] = s.spike[g]; R.clusters[at] = s.label[g];
  }
  R.dconc[j] += s.dconc; R.dtau[j] += s.dtau;
  R.free[j] += s.free; R.groups[j] += groups.size(); R.slab2[j] += s.slab2;
}

// the cluster numbers of the selected spike-and-slab groups, per stem 1, 2, ... in the order of
// the groups
void numberClusters(const Layout& L, Result& R) {
  for (int j = 0; j < L.K; ++j) {
    if (L.kind[j] != 0) continue;
    int next = 0;
    for (int i = 0; i < L.G; ++i) {
      double& c = R.clusters[i + L.G * j];
      if (!ISNA(c) && c > 0) c = ++next;
    }
  }
}

// ---- expectation consistent ----------------------------------------------------------------

// the coordinates of the call: (group, stem) and its position among the active stems of the group
struct Coords {
  std::vector<int> g, j, pos;
  std::vector<std::vector<int>> at;  // per group the coordinate of each active stem
  std::vector<std::vector<int>> stemGroups, stemCoords;  // per Dirichlet-process stem
};

Coords coordsOf(const Layout& L) {
  Coords C;
  C.at.resize(L.G);
  C.stemGroups.resize(L.K); C.stemCoords.resize(L.K);
  for (int i = 0; i < L.G; ++i)
    for (std::size_t a = 0; a < L.act[i].size(); ++a) {
      const int j = L.act[i][a];
      C.at[i].push_back(C.g.size());
      if (L.kind[j] == 2) { C.stemGroups[j].push_back(i); C.stemCoords[j].push_back(C.g.size()); }
      C.g.push_back(i); C.j.push_back(j); C.pos.push_back(a);
    }
  return C;
}

// The factors against the cavities lam = (lh, lt): log Z summed and the moments per coordinate;
// false outside a domain.
bool factorsAt(const Layout& L, const Coords& C, const std::vector<double>& lh, const std::vector<double>& lt,
               double& logZ, std::vector<double>& m1, std::vector<double>& m2, std::vector<double>& m3,
               std::vector<double>& m4, std::vector<double>& vq, std::vector<Stem>* stems = nullptr) {
  if (stems != nullptr) stems->assign(L.K, Stem());
  const int N = C.g.size();
  m1.assign(N, 0.0); m2.assign(N, 0.0); m3.assign(N, 0.0); m4.assign(N, 0.0); vq.assign(N, 0.0);
  logZ = 0.0;
  for (int k = 0; k < N; ++k) {
    const int j = C.j[k];
    if (L.kind[j] == 2) continue;
    const Coord c = L.kind[j] == 0 ? slab1(lh[k], lt[k], L.prob[j], L.tau[j]) : laplace1(lh[k], lt[k], L.lambda[j]);
    if (!c.ok) return false;
    logZ += c.logZ;
    m1[k] = c.m1; m2[k] = c.m2; m3[k] = c.m3; m4[k] = c.m4; vq[k] = c.var;
  }
  for (int j = 0; j < L.K; ++j) {
    if (L.kind[j] != 2 || C.stemCoords[j].empty()) continue;
    const std::vector<int>& ks = C.stemCoords[j];
    std::vector<double> h(ks.size()), t(ks.size());
    for (std::size_t a = 0; a < ks.size(); ++a) { h[a] = lh[ks[a]]; t[a] = lt[ks[a]]; }
    const Stem s = dpStem(h, t, L.conc[j], L.tau[j], stems != nullptr);
    if (!s.ok) return false;
    logZ += s.logZ;
    if (stems != nullptr) (*stems)[j] = s;
    for (std::size_t a = 0; a < ks.size(); ++a) {
      m1[ks[a]] = s.m1[a]; m2[ks[a]] = s.m2[a]; m3[ks[a]] = s.m3[a]; m4[ks[a]] = s.m4[a];
      vq[ks[a]] = s.m2[a] - s.m1[a] * s.m1[a];
    }
  }
  return true;
}

// The Gaussian part under the sites (gh, gt): A_r summed, the marginal moments per coordinate and
// per group its covariance; false outside its domain.
bool gaussAt(const Layout& L, const Coords& C, const std::vector<double>& gh, const std::vector<double>& gt,
             double& Ar, std::vector<double>& mu, std::vector<double>& var, std::vector<Gauss>& gs) {
  const int N = C.g.size();
  mu.assign(N, 0.0); var.assign(N, 0.0);
  gs.assign(L.G, Gauss());
  Ar = 0.0;
  for (int i = 0; i < L.G; ++i) {
    const int n = L.act[i].size();
    if (!n) continue;
    std::vector<double> lt(n), lh(n);
    for (int a = 0; a < n; ++a) { lt[a] = gt[C.at[i][a]]; lh[a] = gh[C.at[i][a]]; }
    if (!gauss(L.P[i], L.b[i], lt, lh, n, gs[i])) return false;
    double quad = 0.0;
    for (int a = 0; a < n; ++a) quad += (L.b[i][a] + lh[a]) * gs[i].mu[a];
    Ar += 0.5 * (n * LOG_2PI - gs[i].logdet) + 0.5 * quad;
    for (int a = 0; a < n; ++a) { mu[C.at[i][a]] = gs[i].mu[a]; var[C.at[i][a]] = gs[i].C[a + a * n]; }
  }
  return true;
}

double As(const std::vector<double>& nh, const std::vector<double>& nt) {
  double s = 0.0;
  for (std::size_t k = 0; k < nh.size(); ++k) s += 0.5 * (LOG_2PI - std::log(nt[k])) + nh[k] * nh[k] / (2.0 * nt[k]);
  return s;
}

void expectationConsistent(const Layout& L, const Coords& C, const NumericMatrix* sitesLam,
                           const NumericMatrix* sitesGam, Result& R) {
  const int N = C.g.size();
  // the inner loop ends where the moments match to tolIn or, within floor, where a Newton step no
  // longer lowers Psi; the outer one where the marginals match to tolOut or, within floor, where F
  // no longer rises: the rounding of an ill-conditioned Gaussian part
  const double tolIn = 1e-10, tolOut = 1e-10, floor = 1e-4;
  // beyond CROSS_MAX groups the Newton matrix leaves out the Dirichlet-process terms between
  // groups, and the inner loop converges linearly
  const int maxIn = 500, maxOut = 1000;
  std::vector<double> gh(N, 0.0), gt(N), lh(N), lt(N), nh(N), nt(N), mu, var, m1, m2, m3, m4, vq;
  std::vector<Gauss> gs;
  double Ar = 0, Aq = 0;
  bool crossing = false;
  for (int j = 0; j < L.K; ++j)
    if (L.kind[j] == 2 && C.stemCoords[j].size() > 1 && (int) C.stemCoords[j].size() <= CROSS_MAX) crossing = true;
  // sites from the last accepted point where they are valid, else at the prior variances
  bool warm = sitesLam != nullptr && sitesGam != nullptr;
  for (int k = 0; k < N; ++k) {
    gt[k] = warm ? (*sitesLam)(C.g[k], C.j[k]) : 1.0 / priorVar(L, C.j[k]);
    gh[k] = warm ? (*sitesGam)(C.g[k], C.j[k]) : 0.0;
    if (!std::isfinite(gt[k]) || !std::isfinite(gh[k])) warm = false;
  }
  if (!warm || !gaussAt(L, C, gh, gt, Ar, mu, var, gs)) {
    for (int k = 0; k < N; ++k) { gt[k] = 1.0 / priorVar(L, C.j[k]); gh[k] = 0.0; }
    if (!gaussAt(L, C, gh, gt, Ar, mu, var, gs)) return;
  }
  for (int k = 0; k < N; ++k) { nt[k] = 1.0 / var[k]; nh[k] = mu[k] / var[k]; lt[k] = nt[k] - gt[k]; lh[k] = nh[k] - gh[k]; }
  if (!factorsAt(L, C, lh, lt, Aq, m1, m2, m3, m4, vq))
    for (int k = 0; k < N; ++k) { lt[k] = nt[k] / 2.0; lh[k] = nh[k] / 2.0; }
  // Cov_q phi and Cov_r phi on the blocks of the Newton systems, rows (eta, -eta^2 / 2) of the
  // coordinates of a block: one per group, or one of all coordinates where a Dirichlet-process
  // stem couples the groups, with its moments between groups
  struct Block { std::vector<int> ks; std::vector<double> Q, Rr; };
  // the moments of one coordinate in the matrix H of size m at rows x (eta), y (eta^2)
  auto covQ = [&](std::vector<double>& H, int m, int x, int y, int k) {
    H[x + x * m] += vq[k];
    const double c13 = -0.5 * (m3[k] - m1[k] * m2[k]);
    H[x + y * m] += c13;
    H[y + x * m] += c13;
    H[y + y * m] += 0.25 * (m4[k] - m2[k] * m2[k]);
  };
  // Cov_r phi of group i into H of size m, its coordinates at rows pos(a) and pos(a) + off
  auto covR = [&](std::vector<double>& H, int m, int i, int off, const std::vector<int>& pos) {
    const int n = L.act[i].size();
    const std::vector<double>& Cg = gs[i].C;
    for (int a = 0; a < n; ++a)
      for (int c = 0; c < n; ++c) {
        const int ka = C.at[i][a], kc = C.at[i][c], x = pos[a], y = pos[c];
        const double Cac = Cg[a + c * n];
        H[x + y * m] += Cac;
        H[x + (off + y) * m] += -Cac * mu[kc];
        H[(off + y) + x * m] += -Cac * mu[kc];
        H[(off + x) + (off + y) * m] += 0.5 * Cac * Cac + mu[ka] * mu[kc] * Cac;
      }
  };
  auto blocks = [&](std::vector<Block>& B) {
    B.clear();
    if (!crossing) {
      for (int i = 0; i < L.G; ++i) {
        const int n = L.act[i].size();
        if (!n) continue;
        Block b;
        b.ks = C.at[i];
        b.Q.assign(4 * n * n, 0.0); b.Rr.assign(4 * n * n, 0.0);
        std::vector<int> pos(n);
        for (int a = 0; a < n; ++a) pos[a] = a;
        covR(b.Rr, 2 * n, i, n, pos);
        for (int a = 0; a < n; ++a) covQ(b.Q, 2 * n, a, n + a, C.at[i][a]);
        B.push_back(b);
      }
      return true;
    }
    std::vector<Stem> stv;
    double Aq2;
    if (!factorsAt(L, C, lh, lt, Aq2, m1, m2, m3, m4, vq, &stv)) return false;
    const int m = 2 * N;
    Block b;
    b.ks.resize(N);
    for (int k = 0; k < N; ++k) b.ks[k] = k;
    b.Q.assign(m * m, 0.0); b.Rr.assign(m * m, 0.0);
    for (int i = 0; i < L.G; ++i) covR(b.Rr, m, i, N, C.at[i]);
    for (int k = 0; k < N; ++k) if (L.kind[C.j[k]] != 2) covQ(b.Q, m, k, N + k, k);
    for (int j = 0; j < L.K; ++j) {
      const std::vector<int>& ks = C.stemCoords[j];
      if (L.kind[j] != 2 || ks.empty()) continue;
      const Stem& st = stv[j];
      const int Gj = ks.size();
      if (st.e11.empty()) { for (int k : ks) covQ(b.Q, m, k, N + k, k); continue; }
      for (int x = 0; x < Gj; ++x)
        for (int y = 0; y < Gj; ++y) {
          const int k = ks[x], q = ks[y];
          if (x == y) { covQ(b.Q, m, k, N + k, k); continue; }
          b.Q[k + q * m] += st.e11[x + Gj * y] - m1[k] * m1[q];
          b.Q[k + (N + q) * m] += -0.5 * (st.e12[x + Gj * y] - m1[k] * m2[q]);
          b.Q[(N + k) + q * m] += -0.5 * (st.e12[y + Gj * x] - m1[q] * m2[k]);
          b.Q[(N + k) + (N + q) * m] += 0.25 * (st.e22[x + Gj * y] - m2[k] * m2[q]);
        }
    }
    B.push_back(b);
    return true;
  };
  auto solveH = [&](std::vector<double>& H, int m, std::vector<double>& rhs) {
    std::vector<double> Hc = H;
    double ridge = 0.0;
    while (!cholesky(Hc, m)) {
      ridge = ridge ? 10.0 * ridge : 1e-12;
      Hc = H;
      for (int a = 0; a < m; ++a) Hc[a + a * m] += ridge * (1.0 + H[a + a * m]);
      if (ridge > 1e6) return false;
    }
    cholSolve(Hc, m, rhs.data());
    return true;
  };
  // inner: minimise Psi(lam) = A_q(lam) + A_r(nu - lam) by Newton on Cov_q phi + Cov_r phi
  auto psi = [&](const std::vector<double>& h, const std::vector<double>& t, double& value) {
    std::vector<double> sh(N), stt(N);
    for (int k = 0; k < N; ++k) { sh[k] = nh[k] - h[k]; stt[k] = nt[k] - t[k]; }
    if (!factorsAt(L, C, h, t, Aq, m1, m2, m3, m4, vq)) return false;
    if (!gaussAt(L, C, sh, stt, Ar, mu, var, gs)) return false;
    value = Aq + Ar;
    return true;
  };
  auto inner = [&]() {
    double value;
    if (!psi(lh, lt, value)) return false;
    for (int it = 0; it < maxIn; ++it) {
      ++R.inner;
      std::vector<double> g(2 * N);
      double worst = 0.0;
      for (int k = 0; k < N; ++k) {
        // E eta^2 of q minus that of r, as the difference of variances and means
        const double dmean = m1[k] - mu[k];
        g[k] = dmean;
        g[N + k] = -0.5 * (vq[k] - var[k] + dmean * (m1[k] + mu[k]));
        worst = std::max(worst, std::max(std::fabs(dmean) / std::sqrt(var[k]), std::fabs(vq[k] / var[k] - 1.0)));
      }
      if (worst <= tolIn) return true;
      std::vector<double> d(2 * N, 0.0);
      std::vector<Block> B;
      if (!blocks(B)) return false;
      for (const Block& b : B) {
        const int m = b.ks.size(), M = 2 * m;
        std::vector<double> H(M * M), rhs(M);
        for (int a = 0; a < M * M; ++a) H[a] = b.Q[a] + b.Rr[a];
        for (int a = 0; a < m; ++a) { rhs[a] = -g[b.ks[a]]; rhs[m + a] = -g[N + b.ks[a]]; }
        if (!solveH(H, M, rhs)) return false;
        for (int a = 0; a < m; ++a) { d[b.ks[a]] = rhs[a]; d[N + b.ks[a]] = rhs[m + a]; }
      }
      double slope = 0.0;
      for (int k = 0; k < 2 * N; ++k) slope += g[k] * d[k];
      double step = 1.0, trial;
      std::vector<double> th(N), tt(N);
      bool moved = false;
      for (int ls = 0; ls < 60; ++ls, step /= 2.0) {
        for (int k = 0; k < N; ++k) { th[k] = lh[k] + step * d[k]; tt[k] = lt[k] + step * d[N + k]; }
        if (psi(th, tt, trial) && trial <= value + 1e-4 * step * slope) { moved = true; break; }
      }
      // no decrease to resolve: the iterate stays, its moments are those of the last evaluation
      if (!moved) return psi(lh, lt, value) && worst <= floor;
      const bool flat = value - trial <= 1e-14 * (1.0 + std::fabs(value)) && worst <= floor;
      lh = th; lt = tt; value = trial;
      if (flat) return true;
    }
    return false;
  };
  // outer: nu where the marginals of r are the Gaussians of nu; the residual in units of the marginal
  auto residual = [&]() {
    double r = 0.0;
    for (int k = 0; k < N; ++k)
      r = std::max(r, std::max(std::fabs(mu[k] - nh[k] / nt[k]) * std::sqrt(nt[k]), std::fabs(std::log(var[k] * nt[k]))));
    return r;
  };
  // the fixed-point step, nu from the marginals `pm`, `pv` of r with the cavities kept: the CCCP
  // step of F(nu) = min_lam Psi_nu(lam) - A_s(nu), which raises F; where it leaves no valid site,
  // the cavities halve towards nu / 2
  auto plainStep = [&](const std::vector<double>& pm, const std::vector<double>& pv) {
    for (int k = 0; k < N; ++k) { nt[k] = 1.0 / pv[k]; nh[k] = pm[k] / pv[k]; }
    std::vector<double> sh(N), stt(N);
    for (int halve = 0; halve <= 60; ++halve) {
      for (int k = 0; k < N; ++k) { sh[k] = nh[k] - lh[k]; stt[k] = nt[k] - lt[k]; }
      if (gaussAt(L, C, sh, stt, Ar, mu, var, gs)) return true;
      for (int k = 0; k < N; ++k) { lh[k] = (lh[k] + nh[k] / 2.0) / 2.0; lt[k] = (lt[k] + nt[k] / 2.0) / 2.0; }
    }
    return false;
  };
  // the Newton step of mu_q(lam*(nu)) = mu_s(nu) at an inner solution: d mu_q / d nu =
  // Cov_q (Cov_q + Cov_r)^-1 Cov_r, d mu_s / d nu = Cov_s, and the cavities move by
  // d lam* = (Cov_q + Cov_r)^-1 Cov_r d nu
  auto newtonOuter = [&](std::vector<double>& dn, std::vector<double>& dl) {
    dn.assign(2 * N, 0.0); dl.assign(2 * N, 0.0);
    std::vector<Block> B;
    if (!blocks(B)) return false;
    for (const Block& b : B) {
      const int m = b.ks.size(), M = 2 * m;
      std::vector<double> K1(M * M);
      for (int a = 0; a < M * M; ++a) K1[a] = b.Q[a] + b.Rr[a];
      if (!cholesky(K1, M)) return false;
      std::vector<double> X = b.Rr;
      cholSolve(K1, M, X.data(), M);
      std::vector<double> J(M * M, 0.0), rhs(M);
      const double one = 1.0, zero = 0.0;
      F77_CALL(dgemm)("N", "N", &M, &M, &M, &one, b.Q.data(), &M, X.data(), &M, &zero, J.data(), &M FCONE FCONE);
      for (int a = 0; a < m; ++a) {
        const int k = b.ks[a];
        const double v = 1.0 / nt[k], ms = nh[k] / nt[k];
        J[a + a * M] -= v;
        J[a + (m + a) * M] += ms * v;
        J[(m + a) + a * M] += ms * v;
        J[(m + a) + (m + a) * M] -= 0.5 * v * v + ms * ms * v;
        rhs[a] = -(mu[k] - ms);
        rhs[m + a] = 0.5 * ((var[k] + mu[k] * mu[k]) - (v + ms * ms));
      }
      // rows and columns scaled to entries of at most 1: the second moments of a coordinate the
      // data pin are many orders below those of the others
      std::vector<double> dr(M, 0.0), dc(M, 0.0);
      for (int c = 0; c < M; ++c)
        for (int a = 0; a < M; ++a) dr[a] = std::max(dr[a], std::fabs(J[a + c * M]));
      for (int a = 0; a < M; ++a) { if (!(dr[a] > 0)) return false; dr[a] = 1.0 / dr[a]; }
      for (int c = 0; c < M; ++c)
        for (int a = 0; a < M; ++a) dc[c] = std::max(dc[c], std::fabs(dr[a] * J[a + c * M]));
      for (int c = 0; c < M; ++c) { if (!(dc[c] > 0)) return false; dc[c] = 1.0 / dc[c]; }
      for (int c = 0; c < M; ++c)
        for (int a = 0; a < M; ++a) J[a + c * M] *= dr[a] * dc[c];
      for (int a = 0; a < M; ++a) rhs[a] *= dr[a];
      if (!luSolve(J, M, rhs.data())) return false;
      for (int c = 0; c < M; ++c) rhs[c] *= dc[c];
      for (int a = 0; a < m; ++a) {
        dn[b.ks[a]] = rhs[a]; dn[N + b.ks[a]] = rhs[m + a];
        double sh = 0.0, st = 0.0;
        for (int c = 0; c < M; ++c) { sh += X[a + c * M] * rhs[c]; st += X[(m + a) + c * M] * rhs[c]; }
        dl[b.ks[a]] = sh; dl[N + b.ks[a]] = st;
      }
    }
    return true;
  };
  // F at an inner solution, and its gradient in nu, E_r phi - E_s phi
  auto freeEnergy = [&]() { return Aq + Ar - As(nh, nt); };
  auto gradient = [&](std::vector<double>& g) {
    g.assign(2 * N, 0.0);
    for (int k = 0; k < N; ++k) {
      const double v = 1.0 / nt[k], ms = nh[k] / nt[k];
      g[k] = mu[k] - ms;
      g[N + k] = -0.5 * ((var[k] + mu[k] * mu[k]) - (v + ms * ms));
    }
  };
  // F is the merit function of the outer loop: the Newton step where it raises F by the Armijo
  // condition, halved down to a sixteenth, else the CCCP step, which always raises it; the loop ends
  // where the moments match to tolOut, or within floor where F no longer rises
  if (!inner()) return;
  double F = freeEnergy();
  for (int out = 0; out < maxOut; ++out) {
    R.outer = out + 1;
    const double res = residual();
    if (res <= tolOut) { R.converged = true; break; }
    std::vector<double> g, dn, dl;
    gradient(g);
    const std::vector<double> snh = nh, snt = nt, slh = lh, slt = lt, smu = mu, svar = var;
    double Fn = F;
    bool moved = false;
    if (newtonOuter(dn, dl)) {
      double slope = 0.0;
      for (int k = 0; k < 2 * N; ++k) slope += g[k] * dn[k];
      for (double step = 1.0; slope > 0 && step >= 1.0 / 16.0 && !moved; step /= 2.0) {
        bool valid = true;
        for (int k = 0; k < N; ++k) {
          nh[k] = snh[k] + step * dn[k]; nt[k] = snt[k] + step * dn[N + k];
          lh[k] = slh[k] + step * dl[k]; lt[k] = slt[k] + step * dl[N + k];
          valid = valid && nt[k] > 0;
        }
        if (!valid || !inner()) continue;
        Fn = freeEnergy();
        moved = Fn >= F + 1e-4 * step * slope;
      }
      if (!moved) { nh = snh; nt = snt; lh = slh; lt = slt; mu = smu; var = svar; }
    }
    if (!moved) {
      if (!plainStep(smu, svar) || !inner()) return;
      Fn = freeEnergy();
      if (Fn - F <= 1e-14 * (1.0 + std::fabs(F)) && res <= floor) { R.converged = true; break; }
    }
    F = Fn;
  }
  if (!R.converged) return;
  // the value, the moments and the reports at the fixed point
  double value;
  {
    std::vector<double> sh(N), stt(N);
    for (int k = 0; k < N; ++k) { sh[k] = nh[k] - lh[k]; stt[k] = nt[k] - lt[k]; }
    if (!factorsAt(L, C, lh, lt, Aq, m1, m2, m3, m4, vq) || !gaussAt(L, C, sh, stt, Ar, mu, var, gs)) return;
    value = Aq + Ar - As(nh, nt);
    for (int k = 0; k < N; ++k) {
      R.Lam[C.g[k] + L.G * C.j[k]] = stt[k];
      R.gam[C.g[k] + L.G * C.j[k]] = sh[k];
    }
  }
  R.logI = value;
  for (int i = 0; i < L.G; ++i) {
    const int n = L.act[i].size();
    for (int a = 0; a < n; ++a) {
      const int ja = L.act[i][a];
      R.mean[i + L.G * ja] = gs[i].mu[a];
      for (int c = 0; c < n; ++c) R.cov[i][ja + L.K * L.act[i][c]] = gs[i].C[a + c * n];
    }
  }
  for (int k = 0; k < N; ++k) {
    const int j = C.j[k];
    if (L.kind[j] == 2) continue;
    const Coord c = L.kind[j] == 0 ? slab1(lh[k], lt[k], L.prob[j], L.tau[j]) : laplace1(lh[k], lt[k], L.lambda[j]);
    reportCoord(L, C.g[k], j, c, R);
  }
  for (int j = 0; j < L.K; ++j) {
    if (L.kind[j] != 2 || C.stemCoords[j].empty()) continue;
    const std::vector<int>& ks = C.stemCoords[j];
    std::vector<double> h(ks.size()), t(ks.size());
    for (std::size_t a = 0; a < ks.size(); ++a) { h[a] = lh[ks[a]]; t[a] = lt[ks[a]]; }
    reportStem(L, j, C.stemGroups[j], dpStem(h, t, L.conc[j], L.tau[j]), R);
  }
  numberClusters(L, R);
  R.ok = true;
}

// ---- exact factors ---------------------------------------------------------------------------

// The integral where the precision couples none of the factors: per group its spike-and-slab stems
// jointly (SlabSum), every other stem alone, every Dirichlet-process stem over its groups, each
// against the linear term and the precision of its own coordinates. `over` where a joint sum left
// its budget of nodes, R then unchanged.
bool exactFactors(const Layout& L, double budget, Result& R, bool& over) {
  const int G = L.G, K = L.K;
  over = false;
  double logI = 0.0;
  Result T = R;
  for (int i = 0; i < G; ++i) {
    const int n = L.act[i].size();
    std::vector<int> slabs;
    for (int a = 0; a < n; ++a) if (L.kind[L.act[i][a]] == 0) slabs.push_back(a);
    if (!slabs.empty()) {
      const std::size_t m = slabs.size();
      std::vector<double> Pf(m * m), h(m), p(m), tau(m);
      for (std::size_t x = 0; x < m; ++x) {
        const int jx = L.act[i][slabs[x]];
        p[x] = L.prob[jx]; tau[x] = L.tau[jx]; h[x] = L.b[i][slabs[x]];
        for (std::size_t y = 0; y < m; ++y) Pf[x + y * m] = L.P[i][slabs[x] + slabs[y] * n];
      }
      const Joint J = slabJoint(Pf, h, p, tau, m, budget);
      if (J.over) { over = true; return false; }
      if (!J.ok) return false;
      logI += J.logZ;
      for (std::size_t x = 0; x < m; ++x) {
        const int jx = L.act[i][slabs[x]], at = i + G * jx;
        T.mean[at] = J.mean[x];
        for (std::size_t y = 0; y < m; ++y) T.cov[i][jx + K * L.act[i][slabs[y]]] = J.cov[x + y * m];
        T.value[at] = J.value[x]; T.spike[at] = 1.0 - J.incl[x]; T.clusters[at] = J.selected[x];
        T.gprob[at] += J.dprob[x]; T.gtau[at] += J.dtau[x]; T.dprob[jx] += J.dprob[x]; T.dtau[jx] += J.dtau[x];
        T.free[jx] += J.incl[x]; T.groups[jx] += 1; T.slab2[jx] += J.slab2[x];
      }
    }
    for (int a = 0; a < n; ++a) {
      const int j = L.act[i][a];
      if (L.kind[j] != 1) continue;
      const Coord c = laplace1(L.b[i][a], L.P[i][a + a * n], L.lambda[j]);
      if (!c.ok) return false;
      logI += c.logZ;
      T.mean[i + G * j] = c.m1;
      T.cov[i][j + K * j] = c.m2 - c.m1 * c.m1;
      reportCoord(L, i, j, c, T);
    }
  }
  for (int j = 0; j < K; ++j) {
    if (L.kind[j] != 2) continue;
    std::vector<int> gr;
    std::vector<double> h, t;
    for (int i = 0; i < G; ++i)
      for (std::size_t a = 0; a < L.act[i].size(); ++a)
        if (L.act[i][a] == j) {
          gr.push_back(i); h.push_back(L.b[i][a]); t.push_back(L.P[i][a + a * L.act[i].size()]);
        }
    if (gr.empty()) continue;
    const Stem st = dpStem(h, t, L.conc[j], L.tau[j]);
    if (!st.ok) return false;
    logI += st.logZ;
    for (std::size_t x = 0; x < gr.size(); ++x) {
      T.mean[gr[x] + G * j] = st.m1[x];
      T.cov[gr[x]][j + K * j] = st.m2[x] - st.m1[x] * st.m1[x];
    }
    reportStem(L, j, gr, st, T);
  }
  T.logI = logI;
  numberClusters(L, T);
  T.ok = T.converged = true;
  R = T;
  return true;
}

// Whether the precision of group i links two of its factors: its spike-and-slab coordinates
// jointly, every other coordinate alone.
bool groupCoupled(const Layout& L, int i) {
  const int n = L.act[i].size();
  for (int a = 0; a < n; ++a)
    for (int c = 0; c < n; ++c) {
      if (a == c || L.P[i][a + c * n] == 0.0) continue;
      if (!(L.kind[L.act[i][a]] == 0 && L.kind[L.act[i][c]] == 0)) return true;
    }
  return false;
}

// The groups one by one, where no Dirichlet-process stem joins two of them: exactly where the
// precision couples none of the factors of a group and its joint sum stays within `budget` nodes,
// else by expectation consistent, every coordinate a factor; `ec` marks the groups of the second.
void byGroup(const Layout& L, const NumericMatrix* sitesLam, const NumericMatrix* sitesGam, double budget,
             std::vector<int>& ec, Result& R) {
  R.logI = 0.0;
  ec.assign(L.G, 0);
  for (int i = 0; i < L.G; ++i) {
    if (L.act[i].empty()) continue;
    Layout S;
    S.G = 1; S.K = L.K; S.kind = L.kind; S.prob = L.prob; S.tau = L.tau; S.lambda = L.lambda; S.conc = L.conc;
    S.act.assign(1, L.act[i]); S.P.assign(1, L.P[i]); S.b.assign(1, L.b[i]);
    Result Ri;
    initResult(S, Ri);
    bool over = false;
    if (groupCoupled(L, i) || !exactFactors(S, budget, Ri, over)) {
      if (!over && !groupCoupled(L, i)) { R.ok = R.converged = false; return; }
      NumericMatrix lam1, gam1;
      const NumericMatrix *l1 = nullptr, *g1 = nullptr;
      if (sitesLam != nullptr && sitesGam != nullptr) {
        lam1 = NumericMatrix(1, L.K); gam1 = NumericMatrix(1, L.K);
        for (int j = 0; j < L.K; ++j) { lam1(0, j) = (*sitesLam)(i, j); gam1(0, j) = (*sitesGam)(i, j); }
        l1 = &lam1; g1 = &gam1;
      }
      initResult(S, Ri);
      expectationConsistent(S, coordsOf(S), l1, g1, Ri);
      ec[i] = 1;
    }
    R.outer = std::max(R.outer, Ri.outer);
    R.inner += Ri.inner;
    if (!Ri.ok) { R.ok = R.converged = false; return; }
    R.logI += Ri.logI;
    R.cov[i] = Ri.cov[0];
    for (int j = 0; j < L.K; ++j) {
      const int at = i + L.G * j;
      R.mean[at] = Ri.mean[j]; R.abs[at] = Ri.abs[j]; R.value[at] = Ri.value[j]; R.spike[at] = Ri.spike[j];
      R.clusters[at] = Ri.clusters[j]; R.gprob[at] = Ri.gprob[j]; R.gtau[at] = Ri.gtau[j];
      R.glambda[at] = Ri.glambda[j]; R.Lam[at] = Ri.Lam[j]; R.gam[at] = Ri.gam[j];
      R.dprob[j] += Ri.dprob[j]; R.dtau[j] += Ri.dtau[j]; R.dlambda[j] += Ri.dlambda[j]; R.dconc[j] += Ri.dconc[j];
      R.free[j] += Ri.free[j]; R.groups[j] += Ri.groups[j]; R.slab2[j] += Ri.slab2[j];
    }
  }
  numberClusters(L, R);
  R.ok = R.converged = true;
}

}  // namespace

// The integral over the sparse coordinates: P a list of the G precisions on all stems, b the
// G x K linear terms, per stem its kind (0 spike and slab, 1 Laplace, 2 Dirichlet process) and
// hyperparameters, `active` the G x K coordinates the data reach, of which those with a row of P
// and an entry of b exactly 0 keep their prior as well. Exact where P couples no two factors and
// the joint spike-and-slab sum of a group stays within `budget` nodes, else by expectation
// consistent from the sites `sitesLam`, `sitesGam` where given.
// [[Rcpp::export]]
List sparse_integral(List P, NumericMatrix b, IntegerVector kind, NumericVector prob, NumericVector tau,
                     NumericVector lambda, NumericVector conc, LogicalMatrix active,
                     Nullable<NumericMatrix> sitesLam, Nullable<NumericMatrix> sitesGam, double budget) {
  Layout L;
  L.G = b.nrow(); L.K = b.ncol();
  L.kind = as<std::vector<int>>(kind);
  L.prob = as<std::vector<double>>(prob); L.tau = as<std::vector<double>>(tau);
  L.lambda = as<std::vector<double>>(lambda); L.conc = as<std::vector<double>>(conc);
  L.act.resize(L.G); L.P.resize(L.G); L.b.resize(L.G);
  for (int i = 0; i < L.G; ++i) {
    NumericMatrix Pi = P[i];
    // a coordinate the data do not see, its row of P and its b exactly 0, keeps its prior
    for (int j = 0; j < L.K; ++j) {
      if (!active(i, j)) continue;
      bool flat = b(i, j) == 0.0;
      for (int c = 0; c < L.K && flat; ++c) flat = Pi(j, c) == 0.0 && Pi(c, j) == 0.0;
      if (!flat) L.act[i].push_back(j);
    }
    const int n = L.act[i].size();
    L.P[i].assign(n * n, 0.0); L.b[i].assign(n, 0.0);
    for (int a = 0; a < n; ++a) {
      L.b[i][a] = b(i, L.act[i][a]);
      for (int c = 0; c < n; ++c) L.P[i][a + c * n] = Pi(L.act[i][a], L.act[i][c]);
    }
  }
  Result R;
  initResult(L, R);
  NumericMatrix lv, gv;
  const NumericMatrix *l0 = nullptr, *g0 = nullptr;
  if (sitesLam.isNotNull() && sitesGam.isNotNull()) {
    lv = NumericMatrix(sitesLam); gv = NumericMatrix(sitesGam);
    if (lv.nrow() == L.G && lv.ncol() == L.K && gv.nrow() == L.G && gv.ncol() == L.K) { l0 = &lv; g0 = &gv; }
  }
  // all groups at once where a Dirichlet-process stem joins two of them, else one by one
  const Coords C = coordsOf(L);
  bool across = false;
  for (int j = 0; j < L.K; ++j) across = across || (L.kind[j] == 2 && C.stemGroups[j].size() > 1);
  std::vector<int> ec(L.G, 0);
  if (across) {
    bool coupled = false, over = false;
    for (int i = 0; i < L.G; ++i) coupled = coupled || groupCoupled(L, i);
    if (coupled || !exactFactors(L, budget, R, over)) {
      if (coupled || over) {
        initResult(L, R);
        expectationConsistent(L, C, l0, g0, R);
        ec.assign(L.G, 1);
      }
    }
  } else {
    byGroup(L, l0, g0, budget, ec, R);
  }
  bool anyEc = false;
  for (int e : ec) anyEc = anyEc || e;
  auto mat = [&](const std::vector<double>& v, int nr, int nc) {
    NumericMatrix m(nr, nc);
    std::copy(v.begin(), v.end(), m.begin());
    return m;
  };
  List cov(L.G);
  for (int i = 0; i < L.G; ++i) cov[i] = mat(R.cov[i], L.K, L.K);
  LogicalVector ecGroups(L.G);
  for (int i = 0; i < L.G; ++i) ecGroups[i] = ec[i];
  return List::create(
      _["ok"] = R.ok, _["converged"] = R.converged, _["coupling"] = anyEc ? "expectationConsistent" : "none",
      _["ec"] = ecGroups, _["logI"] = R.logI,
      _["mean"] = mat(R.mean, L.G, L.K), _["cov"] = cov, _["abs"] = mat(R.abs, L.G, L.K),
      _["values"] = mat(R.value, L.G, L.K), _["spike"] = mat(R.spike, L.G, L.K),
      _["clusters"] = mat(R.clusters, L.G, L.K), _["dprob"] = R.dprob, _["dtau"] = R.dtau,
      _["dlambda"] = R.dlambda, _["dconc"] = R.dconc, _["gprob"] = mat(R.gprob, L.G, L.K),
      _["gtau"] = mat(R.gtau, L.G, L.K), _["glambda"] = mat(R.glambda, L.G, L.K),
      _["stats"] = List::create(_["free"] = R.free, _["groups"] = R.groups, _["slab2"] = R.slab2),
      _["Lam"] = mat(R.Lam, L.G, L.K), _["gam"] = mat(R.gam, L.G, L.K), _["outer"] = R.outer,
      _["inner"] = R.inner);
}

// The factor of one spike-and-slab, Laplace or Dirichlet-process coordinate against the tilts
// (h, t), for the tests of the integrals.
// [[Rcpp::export]]
List sparse_factor(int kind, NumericVector h, NumericVector t, double prob, double tau, double lambda,
                   double conc) {
  if (kind == 2) {
    const Stem s = dpStem(as<std::vector<double>>(h), as<std::vector<double>>(t), conc, tau);
    return List::create(_["ok"] = s.ok, _["logZ"] = s.logZ, _["m1"] = s.m1, _["m2"] = s.m2, _["m3"] = s.m3,
                        _["m4"] = s.m4, _["dconc"] = s.dconc, _["dtau"] = s.dtau, _["spike"] = s.spike,
                        _["labels"] = s.label);
  }
  const Coord c = kind == 0 ? slab1(h[0], t[0], prob, tau) : laplace1(h[0], t[0], lambda);
  return List::create(_["ok"] = c.ok, _["logZ"] = c.logZ, _["m1"] = c.m1, _["m2"] = c.m2, _["m3"] = c.m3,
                      _["m4"] = c.m4, _["abs"] = c.abs, _["dprob"] = c.dprob, _["dtau"] = c.dtau,
                      _["dlambda"] = c.dlambda);
}
